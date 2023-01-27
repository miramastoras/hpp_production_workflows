version 1.0

# This is a task level wdl workflow to use Merqury to get switch error statistics and BED file of FP kmers in the assembly

workflow runMerquryEval {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Eval switch errors and FP kmers with Merqury"
    }
    input {
        File hap1Fasta
        File hap2Fasta
        File kmerTarball
    }
    call combineFA {
        input:
            hap1Fasta=hap1Fasta,
            hap2Fasta=hap2Fasta
    }
    call merqurySwitch {
        input:
            assemblyFasta = assemblyFasta,
            kmerTarball=kmerTarball
    }
    call merqurySpectraCN {
        input:
            assemblyFasta = assemblyFasta,
            kmerTarball=kmerTarball
    }
    output {
        File merqurySwitchTarball=merqurySwitch.outputTarball
        File merqurySpectraCNTarBall=merqurySpectraCN.outputTarball
    }
}

task combineFA {
    input {
        File hap1Fasta
        File hap2Fasta

        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"
    }
    command <<<
          set -o pipefail
          set -e
          set -u
          set -o xtrace

          cat {hap1Fasta} {hap2Fasta} > diploidFasta.fa
    >>>
    output {
      File outputFasta = "diploidFasta.fa"
    }
      runtime {
          memory: memSizeGB + " GB"
          cpu: threadCount
          disks: "local-disk " + diskSizeGB + " SSD"
          docker: dockerImage
          preemptible: 1
      }
}

task merqurySwitch {
    input {
        File assemblyFasta
        File kmerTarball
        File matKmerTarball
        File patKmerTarball
        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        OMP_NUM_THREADS=~{threadCount}

        # get filename
        ASM_ID=$(basename ~{assemblyFasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # initilize command
        cmd=( /opt/merqury/trio/phase_block.sh )

        # link primary asm file
        FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $FILENAME
            mv ${FILENAME%\.gz} asm.fasta
        else
            ln -s ~{assemblyFasta} asm.fasta
        fi
        cmd+=( asm.fasta )

        # extract kmers
        tar xvf ~{kmerTarball} &
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            tar xvf ~{matKmerTarball} &
            tar xvf ~{patKmerTarball} &
        fi
        wait

        cmd+=( $(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            cmd+=( $(basename ~{matKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
            cmd+=( $(basename ~{patKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        fi

        # prep output
        cmd+=( $ASM_ID.merqury_switch )

        # run command
        ${cmd[@]}

        # get output
        tar czvf $ASM_ID.merqury.tar.gz $ASM_ID.merqury*

        # /opt/merqury/trio/phase_block.sh <asm.fasta> <hap1.meryl> <hap2.meryl> <out>


	>>>
	output {
		File outputTarball = glob("*.merqury_switch.tar.gz")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task merqurySpectraCN {
    input {
        File assemblyFasta
        File kmerTarball
        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        OMP_NUM_THREADS=~{threadCount}

        # get filename
        ASM_ID=$(basename ~{assemblyFasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # initilize command
        cmd=( /opt/merqury/eval/spectra-cn.sh )

        # extract kmers
        tar xvf ~{kmerTarball}

        cmd+=( $(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )

        # link primary asm file
        FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $FILENAME
            mv ${FILENAME%\.gz} asm.fasta
        else
            ln -s ~{assemblyFasta} asm.fasta
        fi
        cmd+=( asm.fasta )

        # prep output
        cmd+=( $ASM_ID.merqury_spectracn )

        # run command
        ${cmd[@]}

        # get output
        tar czvf $ASM_ID.merqury_spectracn.tar.gz $ASM_ID.merqury_spectracn*

        # spectra-cn.sh <read.meryl> <asm1.fasta> [asm2.fasta] out-prefix


	>>>
	output {
		File outputTarball = glob("*.merqury_spectracn.tar.gz")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
