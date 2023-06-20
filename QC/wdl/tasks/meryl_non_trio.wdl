version 1.0

# This is a task level wdl workflow to run Meryl for generating either a hybrid hifi + illumina or just illumina readmer database
# from reads extracted from a diploid bam file

workflow runMeryl {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Generate hybrid (hifi + illumina) or just illumina meryl database"
    }
    call Meryl
    output {
        File merylDb = Meryl.merylDb
    }
}

task Meryl{
    input {
        File ilmBam
        File? hifiBam
        Int kmerSize = 21

        String dockerImage = "juklucas/hpp_merqury:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        mkdir output
        ILM_ID=`basename ~{ilmBam} | sed 's/.bam$//'`

        # make illumina meryls
        samtools fastq -@~{threadCount} ~{ilmBam} > output/${ILM_ID}.fq
        meryl count threads=~{threadCount} k=~{kmerSize} output/${ILM_ID}.fq output output/ilm.k~{kmerSize}.meryl

        # make hybrid db if hifi supplied
        if [[ ! -z "~{hifiBam}" ]]; then

          HIFI_ID=`basename ~{hifiBam} | sed 's/.bam$//'`
          samtools fastq -@~{threadCount} ~{hifiBam} > output/${HIFI_ID}.fq
          meryl count threads=~{threadCount} k=~{kmerSize} output/${HIFI_ID}.fq output output/hifi.k~{kmerSize}.meryl

          # remove low freq kmers to avoid overestimating QV 
          meryl greater-than 1 output/hifi.k~{kmerSize}.meryl output output/hifi.k~{kmerSize}.gt1.meryl
          meryl greater-than 1 output/ilm.k~{kmerSize}.meryl output output/ilm.k~{kmerSize}.gt1.meryl
          rm -rf output/hifi.k~{kmerSize}.meryl

          # merge hifi and ilm
          meryl union-max output/ilm.k~{kmerSize}.gt1.meryl output/hifi.k~{kmerSize}.gt1.meryl output hybrid.k~{kmerSize}.gt1.meryl

          # tarball
          tar -zcvf hybrid.k~{kmerSize}.gt1.meryl.tar.gz hybrid.k~{kmerSize}.gt1.meryl

        else # if not hybrid, return just illumina meryl db
          mv output/ilm.k~{kmerSize}.meryl ilm.k~{kmerSize}.meryl
          tar -zcvf ilm.k~{kmerSize}.meryl.tar.gz ilm.k~{kmerSize}.meryl

        fi
    >>>
    output {
        File merylDb=glob("*.meryl.tar.gz")[0]

    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
