version 1.0

import "extract_reads.wdl" as extractReads_t

workflow runMeryl {

    input {
        Array[File] sampleReadsILM
        Array[File] sampleReadsHiFi
        File? referenceFasta
        String identifier
        Int kmerSize = 21

        String dockerImage = "juklucas/hpp_merqury:latest"
        Int threadCount = 32
    }

    # extract reads
    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as ilmReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=256,
                dockerImage="mobinasri/bio_base:v0.2"
        }
    }
    scatter (readFile in sampleReadsHiFi) {
        call extractReads_t.extractReads as hifiReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=256,
                dockerImage="mobinasri/bio_base:v0.2"
        }
    }

    call merylHybrid as makeHybridDB {
        input:
            ilmReads=ilmReadsExtracted.extractedRead,
            hifiReads=hifiReadsExtracted.extractedRead,
            threadCount=threadCount,
            kmerSize=kmerSize,
            identifier=identifier
    }

	output {
		File hybridMerylDB = makeHybridDB.merylDb
	}
}


task merylHybrid {
    input {
        Array[File] ilmReads
        Array[File] hifiReads

        String identifier
        Int memSizeGB = 256
        Int threadCount = 32
        Int diskSizeGB = 256
        Int kmerSize = 21
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        meryl count threads=~{threadCount} k=~{kmerSize} ~{sep=" " ilmReads} output ilm.meryl
        meryl count threads=~{threadCount} k=~{kmerSize} ~{sep=" " hifiReads} output hifi.meryl

        meryl greater-than 1 threads=~{threadCount} ilm.meryl output ilm.gt1.meryl
        meryl greater-than 1 threads=~{threadCount} hifi.meryl output hifi.gt1.meryl

        meryl union-sum threads=~{threadCount} ilm.gt1.meryl hifi.gt1.meryl output ~{identifier}.hybrid.meryl

        tar zcvf ~{identifier}.hybrid.meryl.tar.gz ~{identifier}.hybrid.meryl
	>>>
	output {
		File merylDb= identifier + ".hybrid.meryl.tar.gz"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
