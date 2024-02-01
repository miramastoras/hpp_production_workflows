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
        String outputLabel
        String sampleID

        Int memSizeGB = 128
        Int threadCount = 4
        Int diskSizeGB = 128
        String dockerImage = "pegi3s/bedtools"
    }

	command <<<
      set -eux -o pipefail
      set -o xtrace

      FASTA_FILENAME=$(basename -- "~{Fasta}")

      if [[ $FASTA_FILENAME =~ \.gz$ ]]; then
          cp ~{Fasta} file.fasta.gz
          gunzip -f file.fasta.gz
      else
          cp ~{Fasta} file.fasta
      fi

      bedtools getfasta -fi file.fasta -bed ~{Bed} -fo ~{sampleID}.~{outputLabel}.subBed.fasta

	>>>
	output {
		  File subFasta = glob("*subBed.fasta")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
