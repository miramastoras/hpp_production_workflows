version 1.0

import "extract_reads.wdl" as extractReads_t

workflow longReadAlignmentNoScatter {
    input {
        String aligner
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        File assembly
        File? referenceFasta
        Int preemptible=2
        Int extractReadsDiskSize=256
        String zones
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="tpesout/hpp_base:latest"
        }
    }
    ## align reads to the assembly
    call alignmentBam as alignment{
        input:
            aligner =  aligner,
            preset = preset,
            refAssembly=assembly,
            readFastq_or_queryAssembly = extractReads.extractedRead,
            diskSize = extractReads.fileSizeGB * 3,
            preemptible = preemptible,
            zones = zones
   }
    output {
        File sortedBamFile = mergeBams.mergedBam
        File baiFile = mergeBams.mergedBai
    }

}

task alignmentBam{
    input{
        String aligner
        String preset
        String sampleID
        String suffix=""
        String options=""
        Array[File] readFastq_or_queryAssembly
        File refAssembly
        Int kmerSize=15
        # runtime configurations
        Int memSize=512
        Int threadCount=128
        Int diskSize
        String dockerImage="mobinasri/long_read_aligner:v0.3.3"
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


        if [[ ~{aligner} == "winnowmap" ]]; then
            meryl count k=~{kmerSize} output merylDB ~{refAssembly}
            meryl print greater-than distinct=0.9998 merylDB > repetitive_k~{kmerSize}.txt
            ALIGNER_CMD="winnowmap -W repetitive_k~{kmerSize}.txt"
        elif [[ ~{aligner} == "minimap2" ]] ; then
            ALIGNER_CMD="minimap2 -k ~{kmerSize}"
        else
             echo "UNSUPPORTED ALIGNER (expect minimap2 or winnowmap): ~{aligner}"
             exit 1
        fi

        fileBasename=~{sampleID}
        echo '${ALIGNER_CMD} -a -x ~{preset} ~{options} -t~{threadCount} ~{refAssembly} ~{sep=" " readFastq_or_queryAssembly} | samtools view -h -b > ${fileBasename%.*.*}.bam'

        ${ALIGNER_CMD} -a -x ~{preset} ~{options} -t~{threadCount} ~{refAssembly} ~{sep=" " readFastq_or_queryAssembly} | samtools view -h -b > ${fileBasename%.*.*}.bam

        if [ -z ~{suffix} ]; then
            OUTPUT_FILE=${fileBasename%.*.*}.sorted.bam
        else
            OUTPUT_FILE=${fileBasename%.*.*}.~{suffix}.sorted.bam
        fi
        samtools sort -@~{threadCount} -o ${OUTPUT_FILE} ${fileBasename%.*.*}.bam
        du -s -BG ${OUTPUT_FILE} | sed 's/G.*//' > outputsize.txt
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
        File sortedBamFile = glob("*.sorted.bam")[0]
        Int fileSizeGB = read_int("outputsize.txt")
    }
}
