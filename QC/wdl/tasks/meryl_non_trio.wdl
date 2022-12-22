version 1.0

# This is a task level wdl workflow to run Meryl for generating either a hybrid hifi + illumina or just illumina readmer database
# from reads extracted from a diploid bam file

workflow runMeryl {
    call Meryl
    output {
        File merylDb = Meryl.merylDb
    }
}

task Meryl{
    input {
        File ilmBam
        File? hifiBam

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

        # make illumina meryls
        samtools fastq -@~{threadCount} ~{ilmBam} > output/${ilmBam}.fq
        meryl count threads=~{threadCount} k=21 output/${ilmBam}.fq output output/ilm.k21.meryl
        meryl greater-than 1 output/ilm.k21.meryl output output/ilm.k21.gt1.meryl
        rm -rf output/ilm.k21.meryl

        # make hybrid db if hifi supplied
        if [[ ! -z "~{hifiReads}" ]]; then
          samtools fastq -@~{threadCount} ~{hifiBam} > output/${hifiBam}.fq
          meryl count threads=~{threadCount} k=21 output/${hifiBam}.fq output output/hifi.k21.meryl
          meryl greater-than 1 output/hifi.k21.meryl output output/hifi.k21.gt1.meryl
          rm -rf output/hifi.k21.meryl

          # merge with ilm
          meryl union-max output/ilm.k21.gt1.meryl output/hifi.k21.gt1.meryl output hybrid.k21.gt1.meryl

          # tarball
          tar -zcvf hybrid.k21.gt1.meryl.tar.gz hybrid.k21.gt1.meryl

        else: # if not hybrid, return just illumina meryl db
          mv output/ilm.k21.gt1.meryl ilm.k21.gt1.meryl
          tar -zcvf ilm.k21.gt1.meryl.tar.gz ilm.k21.gt1.meryl

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
