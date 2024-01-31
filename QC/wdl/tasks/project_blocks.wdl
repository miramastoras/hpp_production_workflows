version 1.0

workflow projectBlocks {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Project bed coordinates from one genome assembly to another"
      }

    call project_blocks
    output {
        File projectionBedFile = project_blocks.projectionBedFile
        File projectableBedFile = project_blocks.projectableBedFile
    }
}

task project_blocks{
    input {
        File pafFile
        File bedFile
        String mode='ref2asm'
        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "mobinasri/flagger:latest"
      }

      command <<<
            set -o pipefail
            set -e
            set -u
            set -o xtrace

            pafBasename=$(basename ~{pafFile} | sed 's/.gz$//' | sed 's/.paf$//')
            bedBasename=$(basename ~{bedFile} | sed 's/.gz$//' | sed 's/.bed$//')

            bedtools merge -c 1 -o count -i ~{bedFile} > ${bedBasename}.mrg.bed
            bedtools sort -i ${bedBasename}.mrg.bed > ${bedBasename}.mrg.srt.bed

            python3 /home/programs/src/project_blocks_multi_thread.py \
            --threads ~{threadCount} \
            --mode ~{mode} \
            --paf ~{pafFile} \
            --blocks ${bedBasename}.mrg.srt.bed \
            --outputProjectable ${pafBasename}_${bedBasename}.projectable.bed \
            --outputProjection ${pafBasename}_${bedBasename}.projection.bed
    	>>>

      output {
          File projectionBedFile = glob("*.projection.bed")[0]
          File projectableBedFile = glob("*.projectable.bed")[0]
      }

      runtime{
          memory: memSizeGB + " GB"
          cpu: threadCount
          disks: "local-disk " + diskSizeGB + " SSD"
          docker: dockerImage
      }
}
