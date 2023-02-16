version 1.0

# This is a task level wdl workflow to use Merqury to get switch error statistics and BED file of FP kmers in the assembly

workflow runMerqurySwitch {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Switch switch errors  Merqury"
    }

    call merqurySwitch
    output {
        File merqurySwitchTarball=merqurySwitch.outputSwitchTarball
        File merqurySpectraCNTarball=merqurySwitch.outputSpectraCNTarball
    }
}

task merqurySwitch {
    input {
        File hap1Fasta
        File hap2Fasta
        File kmerTarball
        File matKmerTarball
        File patKmerTarball
        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "miramastoras/merqury:latest"
    }

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        OMP_NUM_THREADS=~{threadCount}

        # initialize commands
        cmdPhase=( /opt/merqury/trio/phase_block.sh )
        cmdSpectra=( /opt/merqury/Switch/spectra-cn.sh  )
        ASM_ID=$(basename ~{hap1Fasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        ##### spectra-cn.sh <read.meryl> <asm1.fasta> [asm2.fasta] out-prefix

        # extract kmers
        tar xvf ~{kmerTarball}

        cmdSpectra+=( $(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )

        # link primary asm hap1 file
        FILENAME=$(basename -- "~{hap1Fasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{hap1Fasta} .
            gunzip $FILENAME
            mv ${FILENAME%\.gz} asmhap1.fasta
        else
            ln -s ~{hap1Fasta} asmhap1.fasta
        fi
        cmdSpectra+=( asmhap1.fasta )

        # link primary asm hap2 file
        FILENAME=$(basename -- "~{hap2Fasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{hap2Fasta} .
            gunzip $FILENAME
            mv ${FILENAME%\.gz} asmhap2.fasta
        else
            ln -s ~{hap1Fasta} asmhap2.fasta
        fi
        cmdSpectra+=( asmhap2.fasta )

        # prep output
        cmdSpectra+=( $ASM_ID.merqury_spectracn )

        #### phase_block.sh <asm.fasta> <hap1.meryl> <hap2.meryl> <out>

        cat asmhap1.fasta asmhap2.fasta > asmdip.fasta
        cmdPhase+=( asmdip.fasta )

        # extract kmers
        tar xvf ~{patKmerTarball}
        tar xvf ~{matKmerTarball}

        cmdPhase+=( $(basename ~{patKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        cmdPhase+=( $(basename ~{matKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )

        # prep output
        cmdPhase+=( $ASM_ID.merqury_switch )

        # run commands
        ${cmdPhase[@]}
        ${cmdSpectra[@]}

        # get output
        tar czvf $ASM_ID.merqury_switch.tar.gz $ASM_ID.merqury_switch*
        tar czvf $ASM_ID.merqury_spectracn.tar.gz $ASM_ID.merqury_spectracn*

	>>>
	output {
		File outputSwitchTarball = glob("*.merqury_switch.tar.gz")[0]
    File outputSpectraCNTarball = glob("*.merqury_spectracn.tar.gz")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
