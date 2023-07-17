version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/find_homozygous_regions.wdl" as findHomozygousRegions_t
import "../tasks/subBamByBed.wdl" as subBamByBed_t
import "../tasks/extract_reads.wdl" as extract_reads_t
import "../tasks/correct_bam.wdl" as correct_bam_t
import "../tasks/deepvariant.wdl" as deepvariant_t
import "../tasks/pepperMarginDeepVariant.wdl" as pmdv_t
import "../tasks/whatsHapPhase.wdl" as whatshap_phase_t
import "../tasks/long_read_aligner_scattered_PhaseHom.wdl" as long_read_aligner_scattered_t
import "../tasks/secphase.wdl" as secphase_t
import "../tasks/concatVcf.wdl" as concatVcf_t
import "../tasks/get_mapq_table.wdl" as get_mapq_t


workflow phasingHomozygous{
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "workflow to phase hifi reads in homozygous regions using ONT UL data"
    }
    input {
        File paternalFasta
        File paternalFastaIndex
        File maternalFasta
        File maternalFastaIndex

        File diploidFaGz

        File allHifiToDiploidBam
        File allHifiToDiploidBai

        File allONTToMatBam
        File allONTToPatBam
        File allONTToMatBai
        File allONTToPatBai

        String sampleName
    }

    ## Align maternal to paternal assembly
    call long_read_aligner_t.alignmentPaf as alignmentPaf{
        input:
            aligner="winnowmap",
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=maternalFasta,
            refAssembly=paternalFasta,
            suffix="mat2pat",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    ## Get Homozygous regions
    call findHomozygousRegions_t.FindHomozygousRegions as findHomozygousRegions{
        input:
            pafFile=alignmentPaf.pafFile,
            minWindowSizeBp=20000,
            extendBp=50000,
            outPrefix=sampleName
    }
    ## subset diploid bamfile to homozygous regions
    call subBamByBed_t.SubBamByBed as subDipBamByHomozygous{
        input:
            Bam=allHifiToDiploidBam,
            Bai=allHifiToDiploidBai,
            Bed=findHomozygousRegions.extendedBed
    }

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignAllToPatScattered{
        input:
            readFiles=[subDipBamByHomozygous.subBam],
            assembly=paternalFasta,
            aligner="winnowmap",
            preset="map-pb",
            sampleName=sampleName,
            sampleSuffix="all2pat.winnowmap",
            options="--cs --eqx -Y -L",
            dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignAllToMatScattered{
        input:
            readFiles=[subDipBamByHomozygous.subBam],
            assembly=maternalFasta,
            aligner="winnowmap",
            preset="map-pb",
            sampleName=sampleName,
            sampleSuffix="all2mat.winnowmap",
            options="--cs --eqx -Y -L",
            dockerImage="mobinasri/long_read_aligner:v0.2"
    }


    ## correct bams for maxDivergence of reads
    call correct_bam_t.correctBam as correctBamMaxDivergencePat {
        input:
            Bam=alignAllToPatScattered.bamFile,
            options="--maxDiv 0.02",
            suffix="maxDiv.02",
            dockerImage="mobinasri/secphase:dev-v0.2.0-hom"

    }

    call correct_bam_t.correctBam as correctBamMaxDivergenceMat {
        input:
            Bam=alignAllToMatScattered.bamFile,
            options="--maxDiv 0.02",
            suffix="maxDiv.02",
            dockerImage="mobinasri/secphase:dev-v0.2.0-hom"

    }

    ## call variants on each bam
    call deepvariant_t.DeepVariant as DeepVariantPat{
        input:
            inputReads=correctBamMaxDivergencePat.correctedBam,
            inputReadsIdx=correctBamMaxDivergencePat.correctedBamIndex,
            assembly=paternalFasta,
            assemblyIndex=paternalFastaIndex,
            sample=sampleName,
            modelType = "PACBIO"
    }
    call deepvariant_t.DeepVariant as DeepVariantMat{
        input:
            inputReads=correctBamMaxDivergenceMat.correctedBam,
            inputReadsIdx=correctBamMaxDivergenceMat.correctedBamIndex,
            assembly=maternalFasta,
            assemblyIndex=maternalFastaIndex,
            sample=sampleName,
            modelType = "PACBIO"
    }

    ## filter variants by GQ
    call pmdv_t.bcftoolsFilter as FilterDVPat{
        input:
          inputVCF=DeepVariantPat.vcfOut,
          excludeExpr="'FORMAT/GQ<=10'",
          applyFilters=""
    }
    call pmdv_t.bcftoolsFilter as FilterDVMat{
        input:
          inputVCF=DeepVariantMat.vcfOut,
          excludeExpr="'FORMAT/GQ<=10'",
          applyFilters=""
    }

    ## Phase variants with UL reads
    call whatshap_phase_t.WhatsHapPhase as WhatsHapPhasePat {
        input:
          vcfFile=FilterDVPat.vcfOut,
          vcfFileIdx=FilterDVPat.vcfOutIdx,
          refFile=paternalFasta,
          refFileIdx=paternalFastaIndex,
          bamFile=allONTToPatBam,
          bamFileIdx=allONTToPatBai,
          outPrefix="phased_Vcf_UL_Pat"
    }
    call whatshap_phase_t.WhatsHapPhase as WhatsHapPhaseMat {
        input:
          vcfFile=FilterDVMat.vcfOut,
          vcfFileIdx=FilterDVMat.vcfOutIdx,
          refFile=maternalFasta,
          refFileIdx=maternalFastaIndex,
          bamFile=allONTToMatBam,
          bamFileIdx=allONTToMatBai,
          outPrefix="phased_Vcf_UL_Mat"
    }

    call concatVcf_t.bcftoolsConcat as bcftoolsConcat {
        input:
          vcf1=WhatsHapPhasePat.phasedVcf,
          vcf2=WhatsHapPhaseMat.phasedVcf
    }

    call secphase_t.runSecPhase as runSecPhase {
        input:
          inputBam=allHifiToDiploidBam,
          diploidAssemblyFastaGz=diploidFaGz,
          phasedVcf=bcftoolsConcat.vcfOut,
          variantBed=findHomozygousRegions.bed
    }

    call get_mapq_t.getMapQTable as getMapQTable {
        input:
          allHifiToMatBam=correctBamMaxDivergenceMat.correctedBam,
          allHifiToMatBai=correctBamMaxDivergenceMat.correctedBamIndex,
          allHifiToPatBam=correctBamMaxDivergencePat.correctedBam,
          allHifiToPatBai=correctBamMaxDivergencePat.correctedBamIndex,
          secPhaseBed=runSecPhase.variantBlocksBed
    }

    call correct_bam_t.correctBam as correctBamSecPhase {
        input:
          Bam=allHifiToDiploidBam,
          mapqTableText=getMapQTable.mapqTable,
          phasingLogText=runSecPhase.outLog,
          suffix="UL_phased",
          options="--maxDiv 0.002"
    }

    output {
        File phasedVcfMat=WhatsHapPhaseMat.phasedVcf
        File phasedVcfPat=WhatsHapPhasePat.phasedVcf

        File asm2asmPaf=alignmentPaf.pafFile

        File homExtendedbed=findHomozygousRegions.extendedBed
        File homBed=findHomozygousRegions.bed

        File dipBamHomozygous=subDipBamByHomozygous.subBam

        File allHomToPatBam=alignAllToPatScattered.bamFile
        File allHomToMatBam=alignAllToMatScattered.bamFile

        File allHomToPatBamMaxDiv=correctBamMaxDivergencePat.correctedBam
        File allHomToMatBamMaxDiv=correctBamMaxDivergenceMat.correctedBam

        File deepVariantPat=DeepVariantPat.vcfOut
        File deepVariantMat=DeepVariantMat.vcfOut

        File deepVariantPatFilt=FilterDVPat.vcfOut
        File deepVariantMatFilt=FilterDVMat.vcfOut

        File secphaseOutLog=runSecPhase.outLog
        File secphaseVariantBlocks=runSecPhase.variantBlocksBed

        File finalPhasedDipBam=correctBamSecPhase.correctedBam
        File finalPhasedDipBai=correctBamSecPhase.correctedBamIndex
    }
}
