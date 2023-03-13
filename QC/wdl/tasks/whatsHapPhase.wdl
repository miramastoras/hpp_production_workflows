version 1.0

# This is a task level wdl workflow to apply a set of variants to an assembly for polishing using bcftools consensus

workflow runWhatsHapPhase {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Phase variants in a vcf using a bamfile"
    }
    call WhatsHapPhase
    output {
        File phasedVcf = WhatsHapPhase.phasedVcf
    }
}

task WhatsHapPhase {
    input {
        File vcfFile
        File vcfFileIdx
        File refFile
        File refFileIdx
        File bamFile
        String outPrefix

        String dockerImage = "tpesout/whatshap:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        whatshap phase -o ~{outPrefix}.vcf.gz --indels -r ~{refFile} ~{vcfFile} ~{bamFile} --ignore-read-groups
    >>>
    output {
        File phasedVcf = "~{outPrefix}.vcf.gz"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
