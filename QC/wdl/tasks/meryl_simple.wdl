version 1.0

workflow runMeryl {

    input {
        Array[File] readFiles
        Int kmerSize
        Boolean compress = false
        Int memSizeGB = 500
        Int threadCount = 32
        Int diskSize =500
        String dockerImage = "juklucas/hpp_merqury:latest"
        String sampleID
    }

    call merylCount as sampleMerylCount {
        input:
            readFiles=readFiles,
            kmerSize=kmerSize,
            compress=compress,
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=diskSize,
            dockerImage=dockerImage,
            sampleID=sampleID
      }

	  output {
		  File sampleMerylDB = sampleMerylCount.merylDb
	     }
}


task merylCount {
    input {
        Array[File] readFiles
        String sampleID
        Int kmerSize=21
        Boolean compress = false
        Int memSizeGB = 42
        Int threadCount = 64
        Int diskSizeGB = 64
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

    String compress_arg = if compress then "compress" else ""

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ID=~{sampleID}

        # generate meryl db
        meryl \
          ~{compress_arg} \
          k=~{kmerSize} \
          threads=~{threadCount} \
          memory=~{memSizeGB} \
          count \
          output $ID.meryl \
          ~{sep=' ' readFiles}

        # package
        tar cvf $ID.meryl.tar $ID.meryl

        # cleanup
        rm -rf $ID.meryl
	>>>
	output {
		File merylDb = glob("*.meryl.tar")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
