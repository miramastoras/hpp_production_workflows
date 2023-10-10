version 1.0

import "../tasks/long_read_aligner_scattered_PhaseHom.wdl" as long_read_aligner_scattered_t
import "./PHARAOH.wdl" as PHARAOH_t
import "../tasks/DeepPolisher.wdl" as deepPolisher_t
import "../tasks/applyPolish.wdl" as applyPolish_t

workflow hprc_DeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Full HPRC polishing pipeline: aligning Hifi reads to raw diploid assembly, phasing homozygous regions with ONT UL, running DeepPolisher and applying polish to raw assemblies"
    }

    input {
        File paternalRawFasta
        File paternalRawFastaIndex
        File maternalRawFasta
        File maternalRawFastaIndex
        File dipRawFastaGz
        File dipRawFastaGzIndex
        File DeepPolisherModelFilesTarGZ

        Array[File] ONTReadsUL
        Array[File] HifiReads

        String DeepPolisherDocker
        String sampleName
        Boolean useMargin=false

        # for minimap2, use k=19 for "map-hifi" and k=15 for "map-ont"
        # for winnowmap, use k=15 for preset "map-pb" and "map-ont"
        # default is minimap2

        String hifiAlignerToUse="winnowmap"
        String ONTAlignerToUse="minimap2"
        String alignerHiFiPreset="map-pb"
        String alignerONTPreset="map-ont"
        String alignerHiFiKmerSize="15"
        String alignerONTKmerSize="15"
    }

    ## Align all hifi reads to diploid assembly
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignHifiToDiploid {
        input:
          assembly=dipRawFastaGz,
          readFiles=HifiReads,
          aligner=hifiAlignerToUse,
          preset=alignerHiFiPreset,
          kmerSize=alignerHiFiKmerSize,
          sampleName=sampleName,
          options="--cs --eqx -L -Y -I8g",
          dockerImage="mobinasri/long_read_aligner:v0.3.3",
          sampleSuffix="hifi.to.diploid.asm"
    }

    ## Align all ONT UL reads to paternal haplotype
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToPat {
        input:
          assembly=paternalRawFasta,
          readFiles=ONTReadsUL,
          aligner=ONTAlignerToUse,
          preset=alignerONTPreset,
          kmerSize=alignerONTKmerSize,
          sampleName=sampleName,
          options="--cs --eqx -L -Y",
          dockerImage="mobinasri/long_read_aligner:v0.3.3",
          sampleSuffix="ONT.to.pat.asm"
    }

    ## Align all ONT UL reads to maternal haplotype
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToMat {
        input:
          assembly=maternalRawFasta,
          readFiles=ONTReadsUL,
          aligner=ONTAlignerToUse,
          preset=alignerONTPreset,
          kmerSize=alignerONTKmerSize,
          sampleName=sampleName,
          options="--cs --eqx -L -Y",
          dockerImage="mobinasri/long_read_aligner:v0.3.3",
          sampleSuffix="ONT.to.mat.asm"
    }

    ## Phase reads in homozygous regions with UL, secphase marker mode in non-homoyzgous regions
    call PHARAOH_t.PHARAOH as PHARAOH {
        input:
          paternalFasta=paternalRawFasta,
          paternalFastaIndex=paternalRawFastaIndex,
          maternalFasta=maternalRawFasta,
          maternalFastaIndex=maternalRawFastaIndex,
          diploidFaGz=dipRawFastaGz,
          allHifiToDiploidBam=alignHifiToDiploid.bamFile,
          allHifiToDiploidBai=alignHifiToDiploid.baiFile,
          allONTToMatBam=alignONTToMat.bamFile,
          allONTToPatBam=alignONTToPat.bamFile,
          allONTToMatBai=alignONTToMat.baiFile,
          allONTToPatBai=alignONTToPat.baiFile,
          sampleName=sampleName,
          useMargin=useMargin,
          PharaohAligner=hifiAlignerToUse,
          PharaohHiFiPreset=alignerHiFiPreset,
          PharaohKmerSize=alignerHiFiKmerSize

    }

    ## Pass final phased hifi alignments to deepPolisher to produce polishing variants
    call deepPolisher_t.runDeepPolisher as DeepPolisher {
        input:
          Bam=PHARAOH.finalPhasedDipBam,
          Bai=PHARAOH.finalPhasedDipBai,
          Fasta=dipRawFastaGz,
          ModelFilesTarGZ=DeepPolisherModelFilesTarGZ,
          dockerImage=DeepPolisherDocker,
          sampleName=sampleName
    }

    ## Apply polishing variants to assemblies
    call applyPolish_t.applyPolish as applyDPPolish {
        input:
          polishingVcf=DeepPolisher.PolisherVcf,
          asmRaw=dipRawFastaGz,
          outPrefix=sampleName
    }

    output {
        File polishedAsm=applyDPPolish.asmPolished
        File finalPhasedDipBam=PHARAOH.finalPhasedDipBam
        File finalPhasedDipBai=PHARAOH.finalPhasedDipBai
        File DeepPolisherVcf=DeepPolisher.PolisherVcf
        File allONTToMatBam=alignONTToMat.bamFile
        File allONTToPatBam=alignONTToPat.bamFile
        File allONTToMatBai=alignONTToMat.baiFile
        File allONTToPatBai=alignONTToPat.baiFile
        File asm2asmPaf=PHARAOH.asm2asmPaf
        File allHomToPatBamMaxDiv=PHARAOH.allHomToPatBamMaxDiv
        File allHomToMatBamMaxDiv=PHARAOH.allHomToMatBamMaxDiv
        File secphaseVariantBlocks=PHARAOH.secphaseVariantBlocks
        File deepVariantPatFilt=PHARAOH.deepVariantPatFilt
        File deepVariantMatFilt=PHARAOH.deepVariantMatFilt
        File phasedVcf=PHARAOH.phasedVcf
        File homExtendedbed=PHARAOH.homExtendedbed
        File homBed=PHARAOH.homBed
    }
}
