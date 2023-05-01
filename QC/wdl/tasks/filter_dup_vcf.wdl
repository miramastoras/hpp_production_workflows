version 1.0

# This is a task level wdl workflow to apply a set of variants to an assembly for polishing using bcftools consensus

workflow runFilterDupsVCF {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Remove duplicate vcf entries"
    }
    call filterDups
    output {
        File
    }
}

task filterDups{
    input {
        File VCF

        String dockerImage = "miramastoras/polishing:latest"
        Int memSizeGB = 64
        Int threadCount = 1
        Int diskSizeGB = 64
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        python3 /opt/useful_scripts/rmdups_vcf.py -i ~{VCF}
    >>>
    output {
        File filtVCF=glob("*filt.vcf.gz")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
