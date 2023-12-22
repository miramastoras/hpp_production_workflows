version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runMeryl {

    input {
        Array[File] sampleReadsILM
        File? referenceFasta
        Int kmerSize = 21
        Boolean compress = false
        Int merylCountMemSizeGB = 42
        Int merylCountThreadCount = 64
        Int merylUnionSumMemSizeGB = 32
        Int merylUnionSumThreadCount = 32
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"

    # extract reads
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

    # get file size of results (for union sum)
    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    # get max file size (for meryl count sum)
    call arithmetic_t.max as sampleReadSizeMax {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    # do the meryl counting
    scatter (readFile in sampleReadsExtracted.extractedRead) {
        call merylCount as sampleMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                compress=compress,
                threadCount=merylCountThreadCount,
                memSizeGB=merylCountMemSizeGB,
                diskSizeGB=sampleReadSizeMax.value * 4,
                dockerImage=dockerImage
        }
    }

    # do the meryl merging
    call merylUnionSum as sampleMerylUnionSum {
        input:
            merylCountFiles=sampleMerylCount.merylDb,
            identifier="sample",
            threadCount=merylUnionSumThreadCount,
            memSizeGB=merylUnionSumMemSizeGB,
            diskSizeGB=sampleReadSize.value * 4,
            dockerImage=dockerImage
    }

	output {
		File sampleMerylDB = sampleMerylUnionSum.merylDb
	}
}


task merylCount {
    input {
        File readFile
        Int kmerSize=21
        Boolean compress = false
        Int memSizeGB = 42
        Int threadCount = 64
        Int diskSizeGB = 64
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

    String compress_arg = if compress then "compress" else ""

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
        OMP_NUM_THREADS=~{threadCount}

        # generate meryl db for each read
        ID=`basename ~{readFile} | sed 's/.gz$//' | sed 's/.f[aq]\(st[aq]\)*$//'`
        meryl ~{compress_arg} k=~{kmerSize} threads=~{threadCount} memory=$((~{memSizeGB}-10)) count output $ID.meryl ~{readFile}

        # package
        tar cvf $ID.meryl.tar $ID.meryl

        # cleanup
        rm -rf $ID.meryl
	>>>
	output {
		File merylDb = glob("*.meryl.tar")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task merylUnionSum {
    input {
        Array[File] merylCountFiles
        String identifier
        Int memSizeGB = 32
        Int threadCount = 32
        Int diskSizeGB = 64
        String dockerImage = "juklucas/hpp_merqury:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # extract meryl dbs
        mkdir extracted/
        cd extracted/
        for m in ~{sep=" " merylCountFiles} ; do
            tar xf $m &
        done
        wait
        cd ../

        # merge meryl dbs
        meryl union-sum output ~{identifier}.meryl extracted/*

        # package
        tar cvf ~{identifier}.meryl.tar ~{identifier}.meryl

        # cleanup
        rm -rf extracted/
	>>>
	output {
		File merylDb = identifier + ".meryl.tar"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
