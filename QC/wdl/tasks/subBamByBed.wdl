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
        String Prefix

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
      BAMFILE=$(basename ~{Bam})
      BAIFILE=$(basename ~{Bai})

      ln -s ~{Bam} ./$BAMFILE
      ln -s ~{Bai} ./$BAIFILE

      samtools view -@ ~{threadCount} -b -h -L ~{Bed} ./$BAMFILE > ~{Prefix}.sub.bam
      samtools index ~{Prefix}.sub.bam

	>>>
	output {
		  File subBam = "~{Prefix}.sub.bam"
		  File subBai = "~{Prefix}.sub.bai"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
