version 1.0

# This is a task level wdl workflow to run Merfin for filtering variants for polishing

workflow runMerfin {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Run merfin to filter variants for polishing"
    }
    input {
        File readmerDBTarball
        String runID
    }
    call CombineVCF {
        input:
            runID = runID
    }
    call MerylHist {
        input:
            readmerDBTarball  = readmerDBTarball
    }
    call GenomeScope {
        input:
            merylHist = MerylHist.hist
    }
    call Merfin {
        input:
            genomeScopeStdout = GenomeScope.genomeScopeStdOut,
            lookupTable       = GenomeScope.lookupTable,
            readmerDBTarball  = readmerDBTarball,
            dipVCF = CombineVCF.dipVCF,
            runID = runID
    }
    output {
        File merfinFilteredVcf = Merfin.filteredVCF
    }
}

task CombineVCF {
    input {
        File hap1VcfFile
        File hap2VcfFile
        String runID

        String dockerImage = "biocontainers/bcftools:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # combine haplotype variant calls into one vcf 
        bcftools concat -a ~{hap1VcfFile} ~{hap2VcfFile} -o ~{runID}.merged_variants.diploid.vcf
        bcftools view -Oz -f "PASS" ~{runID}.merged_variants.diploid.vcf > ~{runID}.merged_variants.diploid.PASS.vcf.gz
    >>>
    output {
        File dipVCF=glob("*merged_variants.diploid.PASS.vcf.gz")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task MerylHist {
    input {
        File readmerDBTarball

        String dockerImage = "juklucas/hpp_merqury:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # untar readmer dbs
        tar xvf ~{readmerDBTarball}
        READMER_DIR=$(basename ~{readmerDBTarball} | sed 's/.gz$//' | sed 's/.tar$//')

        # run genomescope
        meryl histogram ${READMER_DIR} > meryl.hist
    >>>
    output {
        File hist = "meryl.hist"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task GenomeScope{
    input {
        File merylHist

        String dockerImage = "dmolik/genomescope2:latest"
        Int memSizeGB = 128
        Int threadCount = 64
        Int diskSizeGB = 128
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # run genomescope
        xvfb-run /home/genomics/genomescope2.0/genomescope.R -i ~{merylHist} -k 21 -o genomescope_outfiles -p 1 --fitted_hist &> genomescope.stdout
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
        String? extraArgs

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
        PREFIX=${FILENAME%.vcf.gz}

        # untar readmer dbs
        tar xvf ~{readmerDBTarball}
        READMER_DIR=$(basename ~{readmerDBTarball} | sed 's/.gz$//' | sed 's/.tar$//')

        # Pull out peak value from genomescope output
        KCOV=$(grep "kcov" ~{genomeScopeStdout} | cut -d" " -f 4 | cut -d":" -f 2)

        ## Pass optional arguments if extraArgs is set, if not just pass empty string
        if [ -z "~{extraArgs}" ]
        then
            EXTRA_ARGS=""
        else
            EXTRA_ARGS="~{extraArgs}"
        fi

        merfin -polish ${EXTRA_ARGS} -vcf ~{vcfFile} -threads ~{threadCount} -sequence ~{refFasta} -readmers $READMER_DIR -prob ~{lookupTable} -peak $KCOV -output ${PREFIX}.merfin
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
