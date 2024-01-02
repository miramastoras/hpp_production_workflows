version 1.0

import "../tasks/extract_reads.wdl" as extractReads_t
import "../tasks/shard_reads.wdl" as shardReads_t
import "../tasks/arithmetic.wdl" as arithmetic_t

workflow runYakCount {

    input {
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        Array[File] sampleReadsILM
        File? referenceFasta
        Int shardLinesPerFile = 256000000
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
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
                dockerImage="mobinasri/bio_base:v0.2"
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
                dockerImage="mobinasri/bio_base:v0.2"
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
                dockerImage="mobinasri/bio_base:v0.2"
        }
    }

    # get file size of results (for yak counting)
    call arithmetic_t.sum as maternalReadSize {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as paternalReadSize {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }
    call arithmetic_t.sum as sampleReadSize {
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
    output {
      File sampleYak = yakCountSample.outputYak
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
        Int memSizeGB=512
        Int threadCount=32
        Int diskSizeGB=512
        String dockerImage="juklucas/hpp_yak:latest"
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
        preemptible: 1
    }

    output {
        File outputYak = "~{sampleName}.yak"
    }
}
