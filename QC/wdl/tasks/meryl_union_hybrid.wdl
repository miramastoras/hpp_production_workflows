version 1.0

# Takes two DBs and makes a hybrid meryl DB of them

workflow runMerylUnionHybrid {

    input {
        File TarDB1
        File TarDB2
        String sampleID
        Int kmerSize

        String dockerImage = "juklucas/hpp_merqury:latest"
        Int threadCount = 32
    }

    call merylHybrid as makeHybridDB {
        input:
            TarDB1=TarDB1,
            TarDB2=TarDB2,
            threadCount=threadCount,
            kmerSize=kmerSize,
            sampleID=sampleID,
            dockerImage=dockerImage
    }

	output {
		File hybridMerylDB = makeHybridDB.merylDb
	}
}


task merylHybrid {
    input {
        File TarDB1
        File TarDB2

        String sampleID
        Int memSizeGB = 400
        Int threadCount = 32
        Int diskSizeGB = 400
        Int kmerSize
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir extracted

        tar xvf ~{TarDB1} --no-same-owner -C .
        tar xvf ~{TarDB2} --no-same-owner -C .

        BASE1=`basename ~{TarDB1}`
        BASE2=`basename ~{TarDB2}`

        DB_1_BASE="${BASE1%.tar}"
        DB_2_BASE="${BASE2%.tar}"

        meryl greater-than 1 threads=~{threadCount} $DB_1_BASE output DB1.gt1.meryl
        meryl greater-than 1 threads=~{threadCount} $DB_2_BASE output DB2.gt1.meryl

        meryl union-sum threads=~{threadCount} DB1.gt1.meryl DB2.gt1.meryl output ~{sampleID}.k~{kmerSize}.hybrid.meryl

        tar zcvf ~{sampleID}.k~{kmerSize}.hybrid.meryl.tar.gz ~{sampleID}.k~{kmerSize}.hybrid.meryl
	>>>
	output {
		File merylDb= glob("*.hybrid.meryl.tar.gz")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
