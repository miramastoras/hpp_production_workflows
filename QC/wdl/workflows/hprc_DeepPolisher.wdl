version 1.0

import "../tasks/long_read_aligner_scattered_PhaseHom.wdl" as long_read_aligner_scattered_t
import "./phasing_homozygous_v4.wdl" as phasing_homozygous_t
import "./diploid_DeepPolisher.wdl" as diploid_deepPolisher_t
import "../tasks/applyPolish.wdl" as applyPolish_t

workflow hprc_DeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Full HPRC polishing pipeline: aligning Hifi reads to raw diploid assembly,
        phasing homozygous regions with ONT UL, running DeepPolisher and applying polish to raw assemblies"
    }

    input {
        File paternalRawFasta
        File paternalRawFastaIdx
        File maternalRawFasta
        File maternalRawFastaIdx
        File diploidRawFasta
        File diploidRawFastaIdx

        String deepPolisherDocker
        String sampleName
    }

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignHifiToDiploid {
        assembly=diploidRawFasta,
        aligner="winnowmap",
        preset="map-pb",
        sampleName=sampleName,
        options="--cs --eqx -L -Y -I8g",
        dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToPat {
        assembly=paternalRawFasta,
        aligner="winnowmap",
        preset="map-ont",
        sampleName=sampleName,
        options="--cs --eqx -L -Y",
        dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToMat {


    }

    call phasing_homozygous_t as phaseHomozygousRegions {
        paternalFasta=paternalRawFasta,
        paternalFastaIndex=paternalRawFastaIdx,
        maternalFasta=maternalRawFasta,
        maternalFastaIndex=maternalRawFastaIdx,
        diploidFaGz=diploidRawFasta,
        allHifiToDiploidBam=alignHifiToDiploid.bamFile,
        allHifiToDiploidBai=alignHifiToDiploid.baiFile
        allONTToMatBam
    File allONTToPatBam
    File allONTToMatBai
    File allONTToPatBai

    String sampleName
    }

    call deepPolisher_t

    call applyPolish_t

    output {
        File phasedVcfMat=marginPhaseMat.phasedVcf

    }
}
