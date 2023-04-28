
version 1.0

workflow findHomozygousRegions {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "detect homozygous regions from paf alignment of two assemblies"
}
    call FindHomozygousRegions

    output{
        File bed=FindHomozygousRegions.bed
        File extendedBed=FindHomozygousRegions.extendedBed
    }
}

task FindHomozygousRegions{
    input {
        File pafFile
        String minWindowSizeBp
        String extendBp
        String outPrefix

        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        String dockerImage = "mobinasri/secphase:dev-v0.2.0"
    }
    command <<<
        set -eux -o pipefail
        set -o xtrace

        python3 /home/programs/src/find_homozygous_regions.py -p ~{pafFile} -m ~{minWindowSizeBp} -e ~{extendBp} -o ~{outPrefix}
  	>>>
  	output {
  		  File bed = "~{outPrefix}.bed"
  		  File extendedBed = "~{outPrefix}*flanking*.bed"
  	}
      runtime {
          memory: memSizeGB + " GB"
          cpu: threadCount
          disks: "local-disk " + diskSizeGB + " SSD"
          docker: dockerImage
          preemptible: 1
      }
}
