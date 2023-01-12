version 1.0

import "../tasks/dipcall.wdl" as runDipcall
import "../tasks/hapDotPy.wdl" as runHappy

workflow dipcall_happy {

    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "This is a wdl workflow to use dipcall + hap.py to benchmark quality of known variants (ie GIAB) in a diploid assembly"
    }
    input {
        File referenceFasta
        File referenceFastaFai
        File confidenceBedFile
    }

    ## Run dipcall
    call runDipcall.dipcall as dipcall_t {
        input:
            referenceFasta=referenceFasta
    }

    ## Call bedtoolsIntersect to subset happy -f bedfile by dipcall bedfile. Intersections written relative to BED1
    call bedtoolsIntersect {
        input:
            BED1=confidenceBedFile,
            BED2=dipcall_t.outputBED
    }

    ## run happy to evaluate assembly variants against truthset variants from reference
    call runHappy.hapDotPy as happy_t {
        input:
            queryVCF      = dipcall_t.outputVCF,
            assembly      = referenceFasta,
            assemblyIndex = referenceFastaFai,
            bedRegions    = bedtoolsIntersect.outputBED
    }

    output{
        File dipCallhapDotPyVCF    = happy_t.vcfOut
        File dipCallhapDotPyVCFIdx = happy_t.vcfIdxOut
        File dipCallhapDotPyTar    = happy_t.happyTar
    }
}

task bedtoolsIntersect {
    input{
        File BED1
        File BED2

        Int memSizeGB = 8
        Int threadCount = 2
        Int diskSizeGB = 128
        String dockerImage = "mobinasri/flagger"
    }

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        BED1_ID=`basename ~{BED1} | sed 's/.bed$//'`
        BED2_ID=`basename ~{BED2} | sed 's/.bed$//'`

        bedtools intersect  -a ~{BED1} -b ~{BED2} > ${BED1_ID}_intersect_${BED2_ID}.bed

    >>>
    output{
        File outputBED = glob("*intersect*.bed")[0]
    }

    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
