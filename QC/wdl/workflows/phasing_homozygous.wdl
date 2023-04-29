
version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/find_homozygous_regions.wdl" as findHomozygousRegions_t
import "../tasks/subBamByBed.wdl" as subBamByBed_t
import "../tasks/extract_reads.wdl" as extract_reads_t


workflow phasingHomozygous{

    input {
        File paternalFasta
        File maternalFasta
        File allReadsToDiploidBam
        File allReadsToDiploidBai
        String sampleName
    }

    ## Align maternal to paternal assembly
    call long_read_aligner_t.alignmentPaf as alignmentPaf{
        input:
            aligner="winnowmap",
            preset="asm5",
            options="-Y -L --eqx --cs -c",
            readFastq_or_queryAssembly=maternalFasta,
            refAssembly=paternalFasta,
            suffix="mat2pat",
            diskSize=256,
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
            Bam=allReadsToDiploidBam,
            Bai=allReadsToDiploidBai,
            Bed=findHomozygousRegions.extendedBed
    }

    ## convert homozygous reads to fastq
    call extract_reads_t.extractReads as extractReads {
        input:
            readFile=subDipBamByHomozygous.subBam
    }
    ## align all homozygous reads to pat
    call long_read_aligner_t.alignmentBam as alignAllToPat {
        input:
            readFastq_or_queryAssembly=extractReads.extractedRead,
            refAssembly=paternalFasta,
            aligner="winnowmap",
            preset="map-pb",
            suffix="all2pat.winnowmap",
            options="--cs --eqx -Y -L",
            dockerImage="mobinasri/long_read_aligner:v0.2",
            diskSize=256
    }

    ## align all homozygous reads to mat
    call long_read_aligner_t.alignmentBam as alignAllToMat {
        input:
            readFastq_or_queryAssembly=extractReads.extractedRead,
            refAssembly=maternalFasta,
            aligner="winnowmap",
            preset="map-pb",
            suffix="all2mat.winnowmap",
            options="--cs --eqx -Y -L",
            dockerImage="mobinasri/long_read_aligner:v0.2",
            diskSize=256
    }

    output {
        File allHifiToMat=alignAllToMat.sortedBamFile,
        File allHifiToPat=alignAllToPat.sortedBamFile
    }

}

# bcftools view -e 'FORMAT/GQ<=10' -Oz ~{inputVCF} > ~{outputFile}
