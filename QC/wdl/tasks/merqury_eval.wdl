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
        ASM_ID=$(basename ~{hap1Fasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        #### phase_block.sh <asm.fasta> <hap1.meryl> <hap2.meryl> <out>
        # docker run -it -v /scratch/mira:/scratch/mira "miramastoras/merqury:latest" /opt/merqury/trio/phase_block.sh /scratch/mira/test_polishing/merqury_eval/HG002.f1_assembly_v2_genbank.dip.fasta /scratch/mira/test_polishing/merqury_eval/paternal.hapmers.meryl /scratch/mira/test_polishing/merqury_eval/maternal.hapmers.meryl merqury_switch

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

        # get output
        tar czvf $ASM_ID.merqury_switch.tar.gz $ASM_ID.merqury_switch*

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
