version 1.0

workflow runWhatsHapPhase {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Phase variants in a vcf using a bamfile with whatshap"
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
        File bamFileIdx
        String outPrefix

        String dockerImage = "tpesout/whatshap:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 256
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # soft link data and indexes so they are in same place
        REF=$(basename ~{refFile})
        REF_IDX=$(basename ~{refFileIdx})

        cp ~{refFile} ./$REF
        cp ~{refFileIdx} ./$REF_IDX

        VCF=$(basename ~{vcfFile})
        VCF_IDX=$(basename ~{vcfFileIdx})

        cp ~{vcfFile} ./$VCF
        cp ~{vcfFileIdx} ./$VCF_IDX

        BAM=$(basename ~{bamFile})
        BAM_IDX=$(basename ~{bamFileIdx})

        cp ~{bamFile} ./$BAM
        cp ~{bamFileIdx} ./$BAM_IDX

        whatshap phase -o ~{outPrefix}.vcf.gz --indels -r $REF $VCF $BAM --ignore-read-groups
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
