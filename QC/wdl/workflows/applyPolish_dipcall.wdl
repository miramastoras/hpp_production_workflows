version 1.0

import "../tasks/applyPolish.wdl" as apply_polish_wf
import "../tasks/dipcall.wdl" as dipcall_wf

workflow applyPolish_dipcall {

    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Apply polishing edits, run dipcall, get intersected bedfile for happy"
    }
    input {
        File referenceFasta
        File referenceFastaFai
        File confidenceBedFile
        File hap1PolishingVcf
        File hap2PolishingVcf
        File hap1Fasta
        File hap2Fasta
        String sampleID
        String GenotypeToPolish="2"
    }

    call apply_polish_wf.applyPolish as applyPolishHap1 {
        input:
            polishingVcf=hap1PolishingVcf,
            asmRaw=hap1Fasta,
            outPrefix=sampleID,
            HaplotypeLabel="hap1",
            GenotypeToPolish=GenotypeToPolish
    }
    call apply_polish_wf.applyPolish as applyPolishHap2 {
        input:
            polishingVcf=hap2PolishingVcf,
            asmRaw=hap2Fasta,
            outPrefix=sampleID,
            HaplotypeLabel="hap2",
            GenotypeToPolish=GenotypeToPolish
    }
    ## Run dipcall
    call dipcall_wf.dipcall as dipcall_t {
        input:
            referenceFasta=referenceFasta,
            referenceFai=referenceFastaFai,
            assemblyFastaPat=applyPolishHap1.asmPolished,
            assemblyFastaMat=applyPolishHap2.asmPolished
    }

    ## Call bedtoolsIntersect to subset happy -f bedfile by dipcall bedfile. Intersections written relative to BED1
    call bedtoolsIntersect as intersectBeds {
        input:
            BED1=confidenceBedFile,
            BED2=dipcall_t.outputBED
    }
    output{
        File dipCallVCF = dipcall_t.outputVCF
        File happyBed = intersectBeds.outputBED
        File dipCallTar = dipcall_t.outputTarball
        File polishedHap1= applyPolishHap1.asmPolished
        File polishedHap2= applyPolishHap2.asmPolished
    }
}

task bedtoolsIntersect {
    input{
        File BED1
        File BED2
        File? BED3

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

        ## if EXTRA BED is set, do another subset
        if [[ ! -z "~{BED3}" ]]
        then
            BED3_ID=`basename ~{BED3} | sed 's/.bed$//'`
            bedtools intersect -a ${BED1_ID}_intersect_${BED2_ID}.bed -b ~{BED3} > ${BED1_ID}_intersect_${BED2_ID}_intersect_${BED3_ID}.bed
            rm ${BED1_ID}_intersect_${BED2_ID}.bed
        fi

        echo "Size of final bedfile used in hap.py"
        awk '{sum += $3-$2}END{print sum}' ${BED1_ID}_intersect_${BED2_ID}.bed
    >>>
    output{
        File outputBED = glob("*intersect*.bed")[0]
    }

    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
