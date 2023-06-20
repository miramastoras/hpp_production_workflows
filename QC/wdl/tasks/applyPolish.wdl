version 1.0

# This is a task level wdl workflow to apply a set of variants to an assembly for polishing using bcftools consensus

workflow runApplyPolish {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Polish assembly with input vcf"
    }
    call applyPolish
    output {
        File asmPolished = applyPolish.asmPolished
    }
}

task applyPolish{
    input {
        File polishingVcf
        File asmRaw
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

        VCF_FILENAME="~{polishingVcf}"

        FILENAME=$(basename -- "~{polishingVcf}")
        SUFFIX="${FILENAME##*.}"

        if [[ "$SUFFIX" != "gz" ]] ; then
            bcftools view -Oz ~{polishingVcf} > "~{polishingVcf}".gz
            VCF_FILENAME="~{polishingVcf}".gz
        fi

        bcftools index $VCF_FILENAME

        bcftools consensus -f ~{asmRaw} -H 2 $VCF_FILENAME > ~{outPrefix}.polished.fasta
    >>>
    output {
        File asmPolished = "~{outPrefix}.polished.fasta"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
