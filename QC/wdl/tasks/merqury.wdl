version 1.0

workflow runMerqury {
    call merqury
    output {
        File QV = merqury.QV
        File outputTarball = merqury.outputTarball
        File FPkmers = merqury.FPkmers
    }
}

task merqury {
    input {
        File assemblyFasta
        File? altHapFasta
        File kmerTarball
        File? matKmerTarball
        File? patKmerTarball
        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_merqury:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
        OMP_NUM_THREADS=~{threadCount}

        # get filename
        ASM_ID=$(basename ~{assemblyFasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # extract kmers
        tar xvf ~{kmerTarball} --no-same-owner &
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            tar xvf ~{matKmerTarball} --no-same-owner &
            tar xvf ~{patKmerTarball} --no-same-owner &
        fi
        wait

        # initilize command
        cmd=( merqury.sh )
        cmd+=( $(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            cmd+=( $(basename ~{matKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
            cmd+=( $(basename ~{patKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//') )
        fi

        # link primary asm file
        FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $FILENAME
            mv ${FILENAME%\.gz} asm.fasta
        else
            cp ~{assemblyFasta} asm.fasta
        fi
        cmd+=( asm.fasta )

        # link secondary asm file
        if [[ -f "~{altHapFasta}" ]]; then
            FILENAME=$(basename -- "~{altHapFasta}")
            if [[ $FILENAME =~ \.gz$ ]]; then
                cp ~{altHapFasta} .
                gunzip $FILENAME
                mv ${FILENAME%\.gz} altHap.fasta
            else
                cp ~{altHapFasta} altHap.fasta
            fi
            cmd+=( altHap.fasta )
        fi

        # prep output
        cmd+=( $ASM_ID.merqury )

        # run command
        ${cmd[@]}

        # get qv output
        tar czvf $ASM_ID.merqury.tar.gz $ASM_ID.merqury*

        # get FP kmers
        tar czvf $ASM_ID.FPkmers.tar.gz *_only.*

	>>>
	output {
		File QV = glob("*.merqury.qv")[0]
		File outputTarball = glob("*.merqury.tar.gz")[0]
    File FPkmers = glob("*.FPkmers.tar.gz")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
