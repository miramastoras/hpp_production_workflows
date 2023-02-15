version 1.0

workflow runMarginPhase {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "phase variants with margin phase"
    }

    call marginPhase

    output {
        File out_margin_phase_svs = marginPhase.phasedVcf
    }
}

task marginPhase {
    input {
        File vcfFile
        File refFile
        File bamFile
        String sampleName
        String HifiOrONT

        String dockerImage = "miramastoras/marginphase_sv:latest"
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # Set param file based on input hifi or ont read alignments
        if [[ ~{HifiOrONT} =~ Hifi ]]; then
            PARAMS=/opt/margin/params/phase/allParams.phase_vcf.ont.json
        else
            PARAMS=/opt/margin/params/phase/allParams.phase_vcf.pb-hifi.json
        fi

        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/
        margin phase ~{bamFile} ~{refFile} ~{vcfFile} $PARAMS -t ~{threads} -o output/~{sampleName} -M
    >>>
    output {
    File phasedVcf = "output/~{sampleName}.margin_phased.vcf"
    }

    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
