version 1.0

workflow subFastaByBed {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "subset fasta file by coordinates in a bed file "
}
    call SubFastaByBed
    output {
        File subFasta = SubFastaByBed.subFasta
    }
}

task SubFastaByBed {
    input {
        File Fasta
        File Bed

        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        String dockerImage = "pegi3s/bedtools"
    }

	command <<<
      set -eux -o pipefail
      set -o xtrace

      BEDID=`basename ~{Bed} | sed 's/.bed$//'`
      FASTAID=`basename ~{Fasta} | sed 's/.bed$//'`

      bedtools getfasta -fi ~{Fasta} -bed ~{Bed} -fo ${FASTAID}_sub_${BEDID}.fasta

	>>>
	output {
		  File subFasta = glob("*sub*.fasta")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
