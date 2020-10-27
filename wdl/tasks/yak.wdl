version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "sum.wdl" as sum_t

workflow runYakAssemblyStats {

    input {
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        Array[File] sampleReadsILM
        File assemblyFastaPat
        File assemblyFastaMat
        File? referenceFasta
        Int kmerSize = 21
        Int shardLinesPerFile = 256000000
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_yak:latest"
    }

    # extract reads
    scatter (readFile in maternalReadsILM) {
        call extractReads_t.extractReads as maternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in paternalReadsILM) {
        call extractReads_t.extractReads as paternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results (for yak counting)
    call sum_t.sum as maternalReadSize {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call sum_t.sum as paternalReadSize {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }
    call sum_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    # do counting
    call yakCount as yakCountMat {
        input:
            readFiles=maternalReadsExtracted.extractedRead,
            sampleName="mat",
            diskSizeGB=maternalReadSize.value * 2,
            dockerImage=dockerImage
    }
    call yakCount as yakCountPat {
        input:
            readFiles=paternalReadsExtracted.extractedRead,
            sampleName="pat",
            diskSizeGB=paternalReadSize.value * 2,
            dockerImage=dockerImage
    }
    call yakCount as yakCountSample {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            sampleName="sample",
            diskSizeGB=sampleReadSize.value * 2,
            dockerImage=dockerImage
    }

    # get stats
    call yakAssemblyStats {
        input:
            assemblyFastaPat=assemblyFastaPat,
            assemblyFastaMat=assemblyFastaMat,
            patYak=yakCountPat.outputYak,
            matYak=yakCountMat.outputYak,
            sampleYak=yakCountSample.outputYak,
            dockerImage=dockerImage
    }

	output {
		File assemblyStats = yakAssemblyStats.outputTarball
		File maternalYak = yakCountMat.outputYak
		File paternalYak = yakCountPat.outputYak
	}

}


task yakCount {
    input{
        Array[File] readFiles
        String sampleName
        Int bloomSize=37
        # runtime configurations
        Int memSizeGB=64
        Int threadCount=32
        Int diskSizeGB=256
        String dockerImage="tpesout/hpp_yak:latest"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # Kmer counting with https://github.com/lh3/yak.
        yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(cat ~{sep=" " readFiles}) <(cat ~{sep=" " readFiles})
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
    }

    output {
        File outputYak = "~{sampleName}.yak"
    }
}


task yakAssemblyStats {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File patYak
        File matYak
        File sampleYak
        String genomeSize = "3.2g"
        String minSequenceLength = "100k"
        # runtime configurations
        Int memSizeGB = 64
        Int threadCount = 32
        Int diskSizeGB = 256
        String dockerImage = "tpesout/hpp_yak:latest"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # name
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/.[pm]at$//')

        # Computing error rates
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{assemblyFastaPat} > $PREFIX.pat.yak_error_stats.txt
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{assemblyFastaMat} > $PREFIX.mat.yak_error_stats.txt

        # QV
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaPat} > $PREFIX.pat.yak_asm-sr_qv.txt
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaMat} > $PREFIX.mat.yak_asm-sr_qv.txt

        # condense
        tar czvf $PREFIX.yak_stats.tar.gz *txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
    }

    output {
        File outputTarball = glob("*.yak_stats.tar.gz")[0]
    }
}