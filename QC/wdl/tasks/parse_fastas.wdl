  version 1.0

workflow runParseFastas {
    call parseFastas
    output {
      File Hap1RawFasta = parseFastas.hap1RawFasta
      File Hap1RawFastaIndex = parseFastas.hap1RawFastaIndex
      File Hap2RawFasta = parseFastas.hap2RawFasta
      File Hap2RawFastaIndex = parseFastas.hap2RawFastaIndex
      File dipRawFastaGz = parseFastas.dipRawFastaGz
      File dipRawFastaGzIndex = parseFastas.dipRawFastaGzIndex
    }
}

task parseFastas {
    input {
        File hap1Fasta
        File hap2Fasta
        String sampleName
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 48
        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
    }


    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        HAP1_FILENAME=$(basename -- "~{hap1Fasta}")

        HAP2_FILENAME=$(basename -- "~{hap2Fasta}")

        if [[ $HAP1_FILENAME =~ \.gz$ ]]; then
            cp ~{hap1Fasta} hap1.fasta.gz
            gunzip -f hap1.fasta.gz

        else
            ln -s ~{hap1Fasta} hap1.fasta
        fi

        samtools faidx hap1.fasta

        if [[ $HAP2_FILENAME =~ \.gz$ ]]; then
            cp ~{hap2Fasta} hap2.fasta.gz
            gunzip -f hap2.fasta.gz

        else
            ln -s ~{hap2Fasta} hap2.fasta
        fi

        samtools faidx hap2.fasta

        cat hap1.fasta hap2.fasta > diploid.fasta

        bgzip diploid.fasta

        samtools faidx diploid.fasta.gz

    >>>

    output {
        File hap1RawFasta = "hap1.fasta"
        File hap1RawFastaIndex = "hap1.fasta.fai"
        File hap2RawFasta = "hap2.fasta"
        File hap2RawFastaIndex = "hap2.fasta.fai"
        File dipRawFastaGz = "diploid.fasta.gz"
        File dipRawFastaGzIndex = "diploid.fasta.gz.fai"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
