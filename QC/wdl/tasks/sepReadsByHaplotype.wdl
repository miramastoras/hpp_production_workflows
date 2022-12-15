version 1.0

# This is a task level wdl workflow to separate reads from a diploid bamfile by haplotype

workflow runSepReadsByHap {

    call Separate
    output {
        File outputFile = Separate.hap1Bam
        File outputFile = Separate.hap2Bam
    }
}

task Separate{
    input {
        File dipBam
        File hap1Fasta
        File hap2Fasta
        String SampleName

        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    parameter_meta{
     inputVcf: "Reads aligned to assembly. Must be in BAM format."
     SampleName: "Sample name. Will be used in output VCF file."
     outputFileTag: "Output file tag to tag files by type of read data (HiFi/Ont)."
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        FILENAME=$(basename ~{dipBam})
        PREFIX=${FILENAME%.bam}

        samtools index ~{hap1Fasta}
        samtools index ~{hap2Fasta}

        cut -f1-2 ~{hap1Fasta}.fai | awk '{print $1"\t""0""\t"$2-1}' > ~{hap1Fasta}.bed
        cut -f1-2 ~{hap2Fasta}.fai | awk '{print $1"\t""0""\t"$2-1}' > ~{hap2Fasta}.bed

        samtools view -@ ~{threadCount} -bh -L ~{hap1Fasta}.bed ~{dipBam} > ${PREFIX}.hap1.bam
        samtools view -@ ~{threadCount} -bh -L ~{hap2Fasta}.bed ~{dipBam} > ${PREFIX}.hap2.bam

    >>>
    output {
        File hap1Bam = "${PREFIX}.hap1.bam"
        File hap2Bam = "${PREFIX}.hap2.bam"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
