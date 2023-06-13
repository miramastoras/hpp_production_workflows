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
        File paternalRawFastaIndex
        File maternalRawFasta
        File maternalRawFastaIndex
        File diploidRawFasta
        File diploidRawFastaIndex

        Array[File] ONTReadsUL
        Array[File] HifiReads

        String DeepPolisherDocker
        String sampleName
        String DeepPolisherModelFilesTarGZ
    }

    ## Align all hifi reads to diploid assembly
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignHifiToDiploid {
        assembly=diploidRawFasta,
        readFiles=HifiReads,
        aligner="winnowmap",
        preset="map-pb",
        sampleName=sampleName,
        options="--cs --eqx -L -Y -I8g",
        dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    ## Align all ONT UL reads to paternal haplotype
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToPat {
        assembly=paternalRawFasta,
        readFiles=ONTReadsUL,
        aligner="winnowmap",
        preset="map-ont",
        sampleName=sampleName,
        options="--cs --eqx -L -Y",
        dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    ## Align all ONT UL reads to maternal haplotype
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToMat {
        assembly=maternalRawFasta,
        readFiles=ONTReadsUL,
        aligner="winnowmap",
        preset="map-ont",
        sampleName=sampleName,
        options="--cs --eqx -L -Y",
        dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    ## Phase reads in homozygous regions with UL, secphase marker mode in non-homoyzgous regions
    call phasing_homozygous_t as phaseHomozygousRegions {
        paternalFasta=paternalRawFasta,
        paternalFastaIndex=paternalRawFastaIndex,
        maternalFasta=maternalRawFasta,
        maternalFastaIndex=maternalRawFastaIndex,
        diploidFaGz=diploidRawFasta,
        allHifiToDiploidBam=alignHifiToDiploid.bamFile,
        allHifiToDiploidBai=alignHifiToDiploid.baiFile
        allONTToMatBam=alignONTToMat.bamFile,
        allONTToPatBam=alignONTToPat.bamFile,
        allONTToMatBai=alignONTToMat.baiFile,
        allONTToPatBai=alignONTToPat.baiFile,
        sampleName=sampleName
    }

    ## Pass final phased hifi alignments to deepPolisher to produce polishing variants
    call diploid_deepPolisher_t.diploid_DeepPolisher as DeepPolisher {
        hifiPhasedDipBam=phaseHomozygousRegions.finalPhasedDipBam,
        hifiPhasedDipBai=phaseHomozygousRegions.finalPhasedDipBai,
        paternalRawFasta=paternalRawFasta,
        paternalRawFastaIndex=paternalRawFastaIndex,
        maternalRawFasta=maternalRawFasta,
        maternalRawFastaIndex=maternalRawFastaIndex,
        ModelFilesTarGZ=DeepPolisherModelFilesTarGZ,
        DeepPolisherDocker=DeepPolisherDocker,
        sampleName=sampleName
    }

    ## Apply polishing variants to assemblies 
    call applyPolish_t.applyPolish as applyDPPolish {
        hap1PolishingVcf=DeepPolisher.Hap1_DeepPolisherVcf,
        hap2PolishingVcf=DeepPolisher.Hap2_DeepPolisherVcf,
        hap1AsmRaw=paternalRawFasta,
        hap2AsmRaw=maternalRawFasta,
        outPrefix=sampleName
    }

    output {
        File polishedAsmHap1=applyDPPolish.hap1Polished
        File polishedAsmHap2=applyDPPolish.hap2Polished

    }
}
