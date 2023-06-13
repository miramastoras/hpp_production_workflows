version 1.0

workflow runMarginPhaseSmall {

    call marginPhaseSmall

    output {
        File out_margin_phase_small = marginPhaseSmall.phasedVcf
    }
}


task marginPhaseSmall {
    input {
        File vcfFile
        File refFile
        File bamFile
        String sampleName
        String dockerImage
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools index -@ ~{threads} ~{bamFile}
        samtools faidx ~{refFile}
        mkdir output/
        margin phase ~{bamFile} ~{refFile} ~{vcfFile} /opt/margin/params/phase/allParams.phase_vcf.ont.json -t ~{threads} -o output/~{sampleName} -M
        bgzip output/~{sampleName}.phased.vcf
    >>>
    output {
    File phasedVcf = "output/~{sampleName}.phased.vcf.gz"
    }

    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
