version 1.0


workflow concatVCF {

    call bcftoolsConcat

    output {
        File filtVcfOut=filter_t.vcfOut
        File filtVcfOutIdx=filter_t.vcfOutIdx
    }
}

task bcftoolsConcat {
    input {
        File vcf1
        File vcf2

        Int memSizeGB = 8
        Int threadCount = 4
        Int diskSizeGB = 50
        String dockerImage = "kishwars/pepper_deepvariant:r0.8"

    }

    parameter_meta {
        vcf1: "VCF #1 to combine. Must be bgzipped."
        vcf2: "VCF #2 to combine. Must be bgzipped."
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace


        bcftools concat -a ~{}

        tabix -p vcf ~{outputFile}
    >>>

    output {
        File vcfOut = outputFile
        File vcfOutIdx = outputFileIdx
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}
