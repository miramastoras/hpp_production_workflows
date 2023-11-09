version 1.0

workflow subFastaByBed {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "subset bamfile by coordinates in a bed file "
}
    call SubFastaByBed
    output {
        File subFasta = SubFastaByBed.subFasta
    }
}

task SubBamByBed {
    input {
        File Fasta
        File Bed

        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
        String dockerImage = "mobinasri/flagger:latest"
    }

	command <<<
      set -eux -o pipefail
      set -o xtrace

      BAMID=`basename ~{Bam} | sed 's/.bam$//'`
      BEDID=`basename ~{Bed} | sed 's/.bed$//'`

      # softlink bam and index so they are in same directory
      ln -s ~{Bam} ./bamfile.bam
      ln -s ~{Bai} ./baifile.bai

      samtools view -@ ~{threadCount} -b -h -L ~{Bed} ./bamfile.bam > ${BAMID}_sub_${BEDID}.bam
      samtools index ${BAMID}_sub_${BEDID}.bam

	>>>
	output {
		  File subBam = glob("*sub*.bam")[0]
		  File subBai = glob("*sub*.bam.bai")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
