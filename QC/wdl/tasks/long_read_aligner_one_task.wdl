version 1.0

workflow longReadAlignmentOneTask {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Aligns reads to assembly with either minimap2 or winnowmap in a single wdl task"
    }
    call alignmentBam as alignment

    output {
        File sortedBamFile = alignment.sortedBamFile
        File baiFile = alignment.sortedBamFileBai
    }

}

task alignmentBam{
    input{
        String aligner
        String preset
        String sampleID
        String suffix=""
        String options=""
        Array[File] readFiles
        File fastaForAlignment
        File? cramExtractionFasta
        Int minReadLength=0
        Int kmerSize
        # runtime configurations
        Int memSize=512
        Int threadCount=128
        Int diskSize=1024
        String dockerImage="mobinasri/long_read_aligner:v0.3.3"
        Int preemptible=2
        String zones="us-west2-a"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Extract reads ##
        mkdir output
        FASTQ_FOLDER="output"

        for file in ~{sep=' ' readFiles}
        do
            FILENAME=$(basename -- $file)
            PREFIX="${FILENAME%.*}"
            SUFFIX="${FILENAME##*.}"

            if [[ "$SUFFIX" == "bam" ]] ; then
                samtools fastq -@~{threadCount} ${file} > output/${PREFIX}.fq
            elif [[ "$SUFFIX" == "cram" ]] ; then
                if [[ ! -f "~{cramExtractionFasta}" ]] ; then
                    echo "Could not extract $FILENAME, reference file not supplied"
                    exit 1
                fi
                ln -s ~{cramExtractionFasta}
                samtools fastq -@~{threadCount} --reference `basename ~{cramExtractionFasta}` ${file} > output/${PREFIX}.fq
            elif [[ "$SUFFIX" == "gz" ]] ; then
                gunzip -c ${file} > output/${PREFIX}
                ls output
            elif [[ "$SUFFIX" == "fastq" ]] || [[ "$SUFFIX" == "fq" ]] ; then
                ln -s ${file} output/${PREFIX}.fq .
            elif [[ "$SUFFIX" != "fastq" ]] && [[ "$SUFFIX" != "fq" ]] && [[ "$SUFFIX" != "fasta" ]] && [[ "$SUFFIX" != "fa" ]] ; then
                echo "Unsupported file type: ${SUFFIX}"
                exit 1
            fi
        done
        echo "Here"
        ## Filter short reads if param specified
        if minReadLength > 0:
          mkdir output_filtered
          FASTQ_FOLDER="output_filtered"
          minLenKb=$(echo ~{minReadLength} | awk '{printf "%.0f",$1/1e3}')

          for fastq in output/*
          do
              FILENAME=$(basename -- ${fastq})
              PREFIX="${FILENAME%.*}"

              # filter reads shorter than minReadLength
              awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>~{minReadLength}{print a"\n"b"\n"c"\n"$0;}' ${fastq} > output_filtered/${PREFIX}.gt_${minLenKb}kb.fastq
          done
        echo "Here 2"
        ## Run alignment

        if [[ ~{aligner} == "winnowmap" ]]; then
            meryl count k=~{kmerSize} output merylDB ~{fastaForAlignment}
            meryl print greater-than distinct=0.9998 merylDB > repetitive_k~{kmerSize}.txt
            ALIGNER_CMD="winnowmap -W repetitive_k~{kmerSize}.txt"
        elif [[ ~{aligner} == "minimap2" ]] ; then
            ALIGNER_CMD="minimap2 -k ~{kmerSize}"
        else
             echo "UNSUPPORTED ALIGNER (expect minimap2 or winnowmap): ~{aligner}"
             exit 1
        fi

        fileBasename=~{sampleID}
        echo '${ALIGNER_CMD} -a -x ~{preset} ~{options} -t~{threadCount} ~{fastaForAlignment} ${FASTQ_FOLDER}/*.fastq | samtools view -h -b > ${fileBasename%.*.*}.bam'

        ${ALIGNER_CMD} -a -x ~{preset} ~{options} -t~{threadCount} ~{fastaForAlignment} ${FASTQ_FOLDER}/*.fastq | samtools view -h -b > ${fileBasename%.*.*}.bam

        if [ -z ~{suffix} ]; then
            OUTPUT_FILE=${fileBasename%.*.*}.sorted.bam
        else
            OUTPUT_FILE=${fileBasename%.*.*}.~{suffix}.sorted.bam
        fi

        samtools sort -@~{threadCount} -o ${OUTPUT_FILE} ${fileBasename%.*.*}.bam
        samtools index ${OUTPUT_FILE}
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
        File sortedBamFileBai = glob("*.sorted.bam.bai")[0]
        Int fileSizeGB = read_int("outputsize.txt")
    }
}
