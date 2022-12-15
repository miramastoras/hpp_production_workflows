version 1.0

# This is a task level wdl workflow to run Merfin for filtering variants for polishing

workflow runMerfin {
    call GenomeScope
    call Merfin {
        input:
            genomeScopeStdOut = GenomeScope.genomeScopeStdOut,
            lookupTable       = GenomeScope.lookupTable
    }
    output {
        File merfinFilteredVcf = Merfin.filteredVCF
    }
}

task GenomeScope{
    input {
        File merylHist

        String dockerImage = "dmolik/genomescope2"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # run genomescope
        xvfb-run \
        genomescope2.0/genomescope.R \
        -i ~{merylHist} \
        -k 21 \
        -o genomescope_outfiles \
        -p 1 \
        --fitted_hist &> genomescope.stdout
    >>>
    output {
        File genomeScopeStdOut = "genomescope.stdout"
        File lookupTable = "genomescope_outfiles/lookup_table.txt"

    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task Merfin{
    input {
        File genomeScopeStdout
        File readmerDBTarball
        File lookupTable
        File vcfFile
        File refFasta

        String dockerImage = "miramastoras/merfin:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        FILENAME=$(basename ~{vcfFile})
        PREFIX=${FILENAME%.bam}

        # untar readmer dbs
        tar xvf ~{readmerDBTarball}
        READMER_DIR=$(basename ~{readmerDBTarball} | sed 's/.gz$//' | sed 's/.tar$//')

        # Pull out peak value from genomescope output
        KCOV=$(grep "kcov" ~{genomeScopeStdout} | cut -d" " -f 4 | cut -d":" -f 2)

        merfin -polish \
                -vcf ~{vcfFile} \
                -threads ~{threadCount} \
                -sequence ~{refFasta} \
                -readmers $READMER_DIR \
                -prob ~{lookupTable} \
                -peak $KCOV \
                -output ${PREFIX}.merfin
    >>>
    output {
        File filteredVCF=glob("*.merfin.polish.vcf")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
