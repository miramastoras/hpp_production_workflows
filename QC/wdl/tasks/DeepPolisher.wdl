version 1.0

# This is a task level wdl workflow to run DeepPolisher

workflow runDeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Run DeepPolisher on phased reads aligned to one haplotype"
    }
    input {
        File Bam
        File Bai
        File Fasta
        File ModelFilesTarGZ
        String sampleName
        String dockerImage
    }
    call DeepPolisher {
        input:
            Bam=Bam,
            Bai=Bai,
            Fasta=Fasta,
            ModelFilesTarGZ=ModelFilesTarGZ,
            sampleName=sampleName,
            dockerImage=dockerImage
    }
    call DPPostProcess {
        input:
            VCFsTarGz=DeepPolisher.VCFsTarGz
    }
    output {
        File PolisherVcf = DPPostProcess.vcfFile
        File PolisherVcfTbi = DPPostProcess.vcfFileTbi
    }
}

task DeepPolisher{
    input {
        File Bam
        File Bai
        File Fasta
        File ModelFilesTarGZ
        String sampleName

        String dockerImage
        Int memSizeGB = 128
        Int threadCount = 32
        Int diskSizeGB = 256
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        mkdir images
        mkdir vcf

        # softlink bam and index so they are in same directory
        BAMFILE=$(basename ~{Bam})
        BAIFILE=$(basename ~{Bai})

        ln -s ~{Bam} ./${BAMFILE}
        ln -s ~{Bai} ./${BAIFILE}

        # untar model files
        # they need to be tar'd in a folder called "checkpoint" and the name of the tar file needs to match
        # the prefix for the .autofdo_profile, .index and .data-00000-of-00001 files example:
        # tar -zcvf checkpoint-657.tar.gz checkpoint/

        CHECKPOINT=( $(basename ~{ModelFilesTarGZ} | sed 's/.gz$//' | sed 's/.tar$//') )

        echo $CHECKPOINT

        echo checkpoint/${CHECKPOINT}

        tar xvf ~{ModelFilesTarGZ} --no-same-owner

        # Make images
        time polisher make_images --bam ${BAMFILE} --fasta ~{Fasta} --output images/images --cpus ~{threadCount}

        # Inference on images to generate VCFs
        time polisher inference --input_dir images --out_dir vcf/ --checkpoint checkpoint/${CHECKPOINT} --reference_file ~{Fasta} --sample_name ~{sampleName} --cpus ~{threadCount}

        tar -zcvf vcf.tar.gz vcf/

    >>>
    output {
        File VCFsTarGz="vcf.tar.gz"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task DPPostProcess{
    input{
        File VCFsTarGz

        String dockerImage = "miramastoras/polishing:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # untar vcf files
        tar xvf ~{VCFsTarGz} --no-same-owner

        # get vcf names
        VCF_NAMES=$(ls vcf/*vcf.gz | tr '\n' ' ')

        # combine vcfs
        vcf-concat $VCF_NAMES > polisher_output.unsorted.vcf

        # Sort the calls
        bcftools view polisher_output.unsorted.vcf --no-header | vcf-sort > polisher_output.sorted.calls_only

        # Get the header
        bcftools view polisher_output.unsorted.vcf --header > polisher_output.header

        # remove stars for now until bug is fixed
        grep -v "*" polisher_output.sorted.calls_only > tmp; mv tmp polisher_output.sorted.calls_only

        cat polisher_output.header polisher_output.sorted.calls_only > polisher_output.vcf

        bgzip polisher_output.vcf
        tabix -p vcf polisher_output.vcf.gz

        >>>

        output {
            File vcfFile="polisher_output.vcf.gz"
            File vcfFileTbi="polisher_output.vcf.gz.tbi"
        }
        runtime {
            memory: memSizeGB + " GB"
            cpu: threadCount
            disks: "local-disk " + diskSizeGB + " SSD"
            docker: dockerImage
            preemptible: 1
        }
}
