version 1.0

# This is a task level wdl workflow to run Merfin for filtering variants for polishing

workflow runMerfin {
    call genomeScope
    call Merfin
    output {
        File hap1Bam = Separate.hap1Bam
        File hap2Bam = Separate.hap2Bam
    }
}

task genomeScope{
    input {
        File 
}

task Merfin{
    input {
        File dipBam
        File hap1Fai
        File hap2Fai

        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        FILENAME=$(basename ~{dipBam})
        PREFIX=${FILENAME%.bam}

        cut -f1-2 ~{hap1Fai} | awk '{print $1"\t""0""\t"$2-1}' > ~{hap1Fai}.bed
        cut -f1-2 ~{hap2Fai} | awk '{print $1"\t""0""\t"$2-1}' > ~{hap2Fai}.bed

        samtools view -@ ~{threadCount} -bh -L ~{hap1Fai}.bed ~{dipBam} > ${PREFIX}.hap1.bam
        samtools view -@ ~{threadCount} -bh -L ~{hap2Fai}.bed ~{dipBam} > ${PREFIX}.hap2.bam

    >>>
    output {
        File hap1Bam = glob("output/*hap1.bam")[0]
        File hap2Bam = glob("output/*hap2.bam")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
