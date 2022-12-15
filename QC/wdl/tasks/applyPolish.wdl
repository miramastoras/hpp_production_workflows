version 1.0

# This is a task level wdl workflow to apply a set of variants to an assembly for polishing using bcftools consensus

workflow runApplyPolish {

    call applyPolish
    output {
        File hap1AsmPolished = applyPolish.hap1Polished
        File hap2AsmPolished = applyPolish.hap2Polished
    }
}

task applyPolish{
    input {
        File polishingVcf
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

        bcftools consensus -f ~{hap1AsmRaw} -H 1 ~{polishingVcf} > ~{outPrefix}.hap1.polished.fasta
        bcftools consensus -f ~{hap2AsmRaw} -H 1 ~{polishingVcf} > ~{outPrefix}.hap2.polished.fasta
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
