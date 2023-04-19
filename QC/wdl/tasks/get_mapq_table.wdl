version 1.0

workflow getMapQTable {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "make table of mapq values for correctBam"
    }

    call getMapQTable

    output {
        File mapqTable = getMapQTable.mapqTable
    }
}

task getMapQTable {
    input {
        File allHifiToMatBam
        File allHifiToMatBai
        File allHifiToPatBam
        File allHifiToPatBai

        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # softlink bam and index so they are in same directory

        ln -s ~{allHifiToMatBam} mat.bam
        ln -s ~{allHifiToMatBai} mat.bai

        ln -s ~{allHifiToPatBam} pat.bam
        ln -s ~{allHifiToPatBai} pat.bai

        # Set param file based on input hifi or ont read alignments
        samtools view pat.bam | awk '{print $1"\t"$3"\t"$4"\t"$5}' > mapq_table_pat.tsv
        samtools view mat.bam | awk '{print $1"\t"$3"\t"$4"\t"$5}' > mapq_table_mat.tsv
        cat mapq_table_pat.tsv mapq_table_mat.tsv > mapq_table_all.tsv

    >>>
    output {
        File mapqTable = "mapq_table_all.tsv"
    }

    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
