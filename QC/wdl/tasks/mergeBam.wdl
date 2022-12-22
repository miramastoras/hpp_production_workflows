version 1.0

workflow mergeBam {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "merge two bamfiles "
    }
    call Merge
    output {
        File mergedBam = Merge.mergedBam
        File mergedBai = Merge.mergedBai
    }
}

task Merge{
    input {
        File bam1
        File bam2

        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        BAM1ID=`basename ~{bam1} | sed 's/.bam$//'`
        BAM2ID=`basename ~{bam2} | sed 's/.bam$//'`

        samtools merge -@ ~{threadCount} ${BAM1ID}.${BAM2ID}.bam ~{bam1} ~{bam2}
        samtools index ${BAM1ID}.${BAM2ID}.bam

    >>>
    output {
        File mergedBam =  glob("*bam")[0]
        File mergedBai = glob("*.bai")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
