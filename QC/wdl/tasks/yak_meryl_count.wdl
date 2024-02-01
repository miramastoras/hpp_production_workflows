version 1.0

import "../tasks/yak_count.wdl" as yak_count_t
import "../tasks/extract_reads.wdl" as extractReads_t
import "../tasks/arithmetic.wdl" as arithmetic_t


workflow yakMerylCount {

    input {
        Array[File] sampleReadsIlm
        File? referenceFasta
        Int kmerSize
        Int threadCount=32
        String sampleID
    }

    # extract reads
    scatter (readFile in sampleReadsIlm) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=256,
                dockerImage="mobinasri/bio_base:v0.2"
        }
    }
    call arithmetic_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }

    call arithmetic_t.max as sampleReadSizeMax {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }
    # do the meryl counting
    call merylCount as sampleMerylCount {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            kmerSize=kmerSize,
            threadCount=threadCount,
            diskSizeGB=sampleReadSizeMax.value * 4
        }

    call yak_count_t.yakCount as yakCountSample {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            sampleName="sample",
            diskSizeGB=sampleReadSize.value * 2,
            kmerSize=kmerSize
    }
}

task merylCount {
    input {
        Array[File] readFiles
        Int kmerSize=21
        Boolean compress = false
        Int memSizeGB = 128
        Int threadCount = 32
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"
        String sampleID
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
        meryl ~{compress_arg} k=~{kmerSize} threads=~{threadCount} memory=$((~{memSizeGB}-10)) count output ~{sampleID}.meryl ~{sep " " readFile}

        # package
        tar cvf ~{sampleID}.meryl.tar ~{sampleID}.meryl

        # cleanup
        rm -rf ~{sampleID}.meryl
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
