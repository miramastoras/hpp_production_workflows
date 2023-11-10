version 1.0

import "../tasks/applyPolish.wdl" as apply_polish_wf
import "../tasks/dipcall.wdl" as dipcall_wf
import "../tasks/dipcall_happy_eval.wdl" as dipcall_happy_wf

workflow dipcall_happy {

    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Apply polishing edits, run dipcall, get intersected bedfile for happy"
    }
    input {
        File referenceFasta
        File referenceFastaFai
        File confidenceBedFile
        File polishingVcf
        File hap1Fasta
        File hap2Fasta
        String sampleID
    }

    call apply_polish_wf.applyPolish as applyPolishHap1 {
        input:
            polishingVcf=polishingVcf,
            asmRaw=hap1Fasta,
            outPrefix=sampleID,
            HaplotypeLabel="hap1"
    }
    call apply_polish_wf.applyPolish as applyPolishHap2 {
        input:
            polishingVcf=polishingVcf,
            asmRaw=hap2Fasta,
            outPrefix=sampleID,
            HaplotypeLabel="hap2"
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
    call dipcall_happy_wf.bedtoolsIntersect as intersectBeds {
        input:
            BED1=confidenceBedFile,
            BED2=dipcall_t.outputBED
    }
    output{
        File dipCallVCF = dipcall_t.outputVCF
        File happyBed = intersectBeds.outputBED
        File dipCallTar = dipcall_t.outputTarball
    }
}