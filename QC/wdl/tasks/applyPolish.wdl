version 1.0

# This is a task level wdl workflow to apply a set of polishing variants to a diploid assembly using bcftools consensus

workflow runApplyPolish {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Polish assembly with input vcf"
    }
    call applyPolish
    output {
        File hap1AsmPolished = applyPolish.hap1Polished
        File hap2AsmPolished = applyPolish.hap2Polished
    }
}

task applyPolish{
    input {
        File hap1PolishingVcf
        File hap2PolishingVcf
        File hap1AsmRaw
        File hap2AsmRaw
        String outPrefix

        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        H1_VCF_FILENAME="~{hap1PolishingVcf}"
        H2_VCF_FILENAME="~{hap2PolishingVcf}"

        H1_FILENAME=$(basename -- "~{hap1PolishingVcf}")
        H2_FILENAME=$(basename -- "~{hap2PolishingVcf}")

        H1_SUFFIX="${H1_FILENAME##*.}"
        H2_SUFFIX="${H2_FILENAME##*.}"

        if [[ "$H1_SUFFIX" != "gz" ]] ; then
            bcftools view -Oz ~{hap1PolishingVcf} > "~{hap1PolishingVcf}".gz
            H1_VCF_FILENAME="~{hap1PolishingVcf}".gz
        fi

        if [[ "$H2_SUFFIX" != "gz" ]] ; then
            bcftools view -Oz ~{hap2PolishingVcf} > "~{hap2PolishingVcf}".gz
            H2_VCF_FILENAME="~{hap2PolishingVcf}".gz
        fi

        bcftools index $H1_VCF_FILENAME
        bcftools index $H2_VCF_FILENAME

        bcftools consensus -f ~{hap1AsmRaw} -H 2 $H1_VCF_FILENAME > ~{outPrefix}.hap1.polished.fasta
        bcftools consensus -f ~{hap2AsmRaw} -H 2 $H2_VCF_FILENAME > ~{outPrefix}.hap2.polished.fasta
    >>>
    output {
        File hap1Polished = "~{outPrefix}.hap1.polished.fasta"
        File hap2Polished = "~{outPrefix}.hap2.polished.fasta"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
