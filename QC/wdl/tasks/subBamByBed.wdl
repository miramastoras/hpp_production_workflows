version 1.0

workflow subBamByBed {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "subset bamfile by coordinates in a bed file "
}
    call SubBamByBed
    output {
        File subBam = SubBamByBed.subBam
        File subBamBai = SubBamByBed.subBai
    }
}

task SubBamByBed {
    input {
        File Bam
        File Bai
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
      BAMFILE=$(basename ~{BAM})
      BAIFILE=$(basename ~{BAI})

      ln -s ~{Bam} ./$BAMFILE
      ln -s ~{Bai} ./$BAIFILE

      samtools view -@ ~{threadCount} -b -h -L ~{Bed} ~{Bai} > ${BAMID}_sub_${BEDID}.bam
      samtools index ${BAMID}_sub_${BEDID}.bam

	>>>
	output {
		  File subBam = glob("*sub*.bam")[0]
		  File subBai = glob("*sub*.bai")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
