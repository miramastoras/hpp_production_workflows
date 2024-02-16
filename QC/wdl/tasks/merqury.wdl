version 1.0

workflow runMerqury {
    call merqury
    output {
        File QV = merqury.QV
        File outputTarball = merqury.outputTarball
        File altHapFPkmers = merqury.altHapFPkmers
        File asmFPkmers = merqury.asmFPkmers
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
        ASM_ID=$(basename ~{assemblyFasta} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//' | sed 's/[._]Hap1//' | sed 's/[._]Hap2//' | sed 's/[._]hap1//' | sed 's/[._]hap2//' )
        KMER_ID=$(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar//' )

        echo `printing kmer ID `
        echo `$KMER_ID`
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            MAT_KMER_ID=$(basename ~{matKmerTarball} | sed 's/.gz$//' | sed 's/.tar//' )
            PAT_KMER_ID=$(basename ~{patKmerTarball} | sed 's/.gz$//' | sed 's/.tar//' )
        fi

        # extract kmers
        tar xvf ~{kmerTarball} --no-same-owner -C . &
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            tar xvf ~{matKmerTarball} --no-same-owner -C . &
            tar xvf ~{patKmerTarball} --no-same-owner -C .
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

        echo "printing cmd"
        echo $cmd

        echo "printing kmer id"
        echo $KMER_ID

        echo "printing directory"
        echo `ls .`

        # run command
        ${cmd[@]}

        mkdir output

        cp altHap_only.bed $ASM_ID.merqury.altHap_only.bed
        cp asm_only.bed $ASM_ID.merqury.asm_only.bed

        # get qv output
        tar czvf $ASM_ID.merqury.tar.gz $ASM_ID.merqury*

	>>>
	output {
		File QV = glob("*.merqury.qv")[0]
		File outputTarball = glob("*.merqury.tar.gz")[0]
    File altHapFPkmers = "altHap_only.bed"
    File asmFPkmers = "asm_only.bed"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
