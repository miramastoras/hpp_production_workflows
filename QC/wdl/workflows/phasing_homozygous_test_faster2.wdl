version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/find_homozygous_regions.wdl" as findHomozygousRegions_t
import "../tasks/subBamByBed.wdl" as subBamByBed_t
import "../tasks/extract_reads.wdl" as extract_reads_t
import "../tasks/correct_bam.wdl" as correct_bam_t
import "../tasks/deepvariant.wdl" as deepvariant_t
import "../tasks/pepperMarginDeepVariant.wdl" as pmdv_t
import "../tasks/whatsHapPhase.wdl" as whatshap_phase_t
import "../tasks/long_read_aligner_scattered.wdl" as long_read_aligner_scattered_t


workflow phasingHomozygous{

    input {
        File paternalFasta
        File paternalFastaIndex
        File maternalFasta
        File maternalFastaIndex

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

    call splitbamContigWise as splitbamContigWisePat{
        input:
            assemblyFasta = paternalFasta,
            bam = allONTToPatBam,
            bamIndex = allONTToPatBai,
            splitNumber = 16,
            threadCount = 16,
            diskSize = 2 * ceil(size(allONTToPatBam, "GB")) + 64
    }
    scatter (part in zip(splitbamContigWisePat.splitBams, splitbamContigWisePat.splitBeds)) {
        call whatshap_phase_t.WhatsHapPhase as WhatsHapPhasePat {
            input:
              vcfFile=FilterDVPat.vcfOut,
              vcfFileIdx=FilterDVPat.vcfOutIdx,
              refFile=paternalFasta,
              refFileIdx=paternalFastaIndex,
              bamFile=allONTToPatBam,
              bamFileIdx=allONTToPatBai,
              outPrefix="phased_Vcf_UL_Pat",
              diskSizeGB = 2 * ceil(size(part.left, "GB")) + 64,
              }
    }
    call mergeVcf as mergeVcfPat{
        input:
            vcfGzFiles = WhatsHapPhasePat.phasedVcf,
            outputName = basename("${allONTToPatBam}", ".bam")
    }

    call splitbamContigWise as splitbamContigWiseMat{
        input:
            assemblyFasta = maternalFasta,
            bam = allONTToMatBam,
            bamIndex = allONTToMatBai,
            splitNumber = 16,
            threadCount = 16,
            diskSize = 2 * ceil(size(allONTToMatBam, "GB")) + 64
    }
    scatter (part in zip(splitbamContigWiseMat.splitBams, splitbamContigWiseMat.splitBeds)) {
        call whatshap_phase_t.WhatsHapPhase as WhatsHapPhaseMat {
            input:
              vcfFile=FilterDVMat.vcfOut,
              vcfFileIdx=FilterDVMat.vcfOutIdx,
              refFile=maternalFasta,
              refFileIdx=maternalFastaIndex,
              bamFile=allONTToMatBam,
              bamFileIdx=allONTToMatBai,
              outPrefix="phased_Vcf_UL_Mat",
              diskSizeGB = 2 * ceil(size(part.left, "GB")) + 64,
              }
    }
    call mergeVcf as mergeVcfMat{
        input:
            vcfGzFiles = WhatsHapPhaseMat.phasedVcf,
            outputName = basename("${allONTToMatBam}",".bam")
    }
    output {
        File phasedVcfMat=mergeVcfMat.vcfGz
        File phasedVcfPat=mergeVcfPat.vcfGz
    }
}

task splitbamContigWise{
    input{
        File assemblyFasta
        File bam
        File bamIndex
        Int splitNumber
        Int memSize=32
        Int threadCount
        Int diskSize=512
        String dockerImage="mobinasri/flagger:v0.2"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        ## unzip fasta file and produce its index file
        ASSEMBLY_NAME=$(basename ~{assemblyFasta})
        ASSEMBLY_PREFIX=${ASSEMBLY_NAME%%.fa}
        ln -f ~{assemblyFasta} > ${ASSEMBLY_PREFIX}.fa
        samtools faidx ${ASSEMBLY_PREFIX}.fa

        ## hard link the bam and bai files to the working directory
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln -f ~{bam} > ${BAM_PREFIX}.bam
        ln -f ~{bamIndex} > ${BAM_PREFIX}.bam.bai

        ## make a bed file that covers the whole assembly
        cat ${ASSEMBLY_PREFIX}.fa.fai | awk '{print $1"\t"0"\t"$2}' > ${ASSEMBLY_PREFIX}.bed

        ## split the bed file of the whole assembly into multiple bed files
        mkdir split_beds split_bams
        python3 ${SPLIT_BED_CONTIG_WISE_PY} --bed ${ASSEMBLY_PREFIX}.bed --n ~{splitNumber} --dir split_beds --prefix ${ASSEMBLY_PREFIX}

        ## make a separate bam for each bed file
        n=$(ls split_beds/ | wc -l)
        seq 1 ${n} | xargs -I {} -n 1 -P ~{threadCount} sh -c "samtools view -h -b -L split_beds/${ASSEMBLY_PREFIX}_{}.bed ${BAM_PREFIX}.bam > split_bams/${BAM_PREFIX}_{}.bam"
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        Array[File] splitBams = glob("split_bams/*.bam")
        Array[File] splitBeds = glob("split_beds/*.bed")
    }
}

task mergeVcf{
    input{
        Array[File] vcfGzFiles
        String outputName
        # runtime configurations
        Int memSize=32
        Int threadCount=16
        Int diskSize=512
        String dockerImage="mobinasri/bio_base:v0.1"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        mkdir vcf_files
        ln ~{sep=" " vcfGzFiles} vcf_files
        files=(vcf_files/*)

        ## make a header for the merged vcf file
        zcat ${files[0]} | awk 'substr($0,1,1) == "#"' | awk -v colname="~{outputName}" '{if ($1 == "#CHROM"){ for(i =1; i < 10; i++){printf("%s\t",$i)}; printf("%s\n",colname)} else {print $0}}' > merged.vcf
        zcat vcf_files/*.vcf.gz | awk 'substr($0,1,1) != "#"' >> merged.vcf

        ## sort the merged vcf file and produce the final gzipped vcf file
        mkdir final_vcf
        bcftools sort -o final_vcf/~{outputName}.vcf.gz -Oz merged.vcf
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File vcfGz = "final_vcf/~{outputName}.vcf.gz"
    }

}
