version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "filter_hifi_adapter.wdl" as adapter_t
import "tar.wdl" as tar_t

workflow runTrioHifiasm{
    input {
        File paternalYak
        File maternalYak
        Array[File] childReadsHiFi
        Array[File] childReadsONT=[]
        Int? homCov
        Int minOntReadLength=100000
        Int minHiFiReadLength=1
        String childID
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        Boolean filterAdapters
        Int removeLastFastqLines
        Int memSizeGB
        Int threadCount
        Int preemptible
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "quay.io/masri2019/hpp_hifiasm:latest"
        String zones = "us-west2-a"
    }

    scatter (readFile in childReadsHiFi) {
        call extractReads_t.extractReads as childReadsHiFiExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
        # filter short HiFi reads
        call filterShortReads as extractLongHiFiReads{
            input:
                readFastq = childReadsHiFiExtracted.extractedRead,
                diskSizeGB = fileExtractionDiskSizeGB,
                minReadLength = minHiFiReadLength
        }
        if (filterAdapters){
            call adapter_t.cutadapt as filterAdapterHiFi {
                input:
                    readFastq = extractLongHiFiReads.longReadFastq,
                    removeLastLines = removeLastFastqLines,
                    diskSizeGB = fileExtractionDiskSizeGB
            } 
        }
        File hifiProcessed = select_first([filterAdapterHiFi.filteredReadFastq, extractLongHiFiReads.longReadFastq])
    }

    # if ONT reads are provided
    if ("${sep="" childReadsONT}" != ""){
        scatter (readFile in childReadsONT) {
            call extractReads_t.extractReads as childReadsOntExtracted {
                input:
                    readFile=readFile,
                    referenceFasta=referenceFasta,
                    memSizeGB=4,
                    threadCount=4,
                    diskSizeGB=fileExtractionDiskSizeGB,
                    dockerImage=dockerImage
            }
            # filter ONT reads to get UL reads
            call filterShortReads as extractUltraLongReads{
                input:
                    readFastq = childReadsOntExtracted.extractedRead,
                    diskSizeGB = fileExtractionDiskSizeGB,
                    minReadLength = minOntReadLength
            }
         }
         call arithmetic_t.sum as childReadULSize {
             input:
                 integers=extractUltraLongReads.fileSizeGB
         }
    }

    call arithmetic_t.sum as childReadHiFiSize {
        input:
            integers=extractLongHiFiReads.fileSizeGB
    }

    # if no ONT data is provided then it would be zero
    Int readULSize = select_first([childReadULSize.value, 0])

    call trioHifiasm {
        input:
            paternalYak=paternalYak,
            maternalYak=maternalYak,
            childReadsHiFi=hifiProcessed,
            childReadsUL=extractUltraLongReads.longReadFastq, # optional argument
            homCov = homCov,
            childID=childID,
            extraOptions=hifiasmExtraOptions,
            inputBinFilesTarGz=inputBinFilesTarGz,
            memSizeGB=memSizeGB,
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize.value + readULSize) * 2.5) + 64,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    output {
        File outputPaternalGfa = trioHifiasm.outputPaternalGfa
        File outputMaternalGfa = trioHifiasm.outputMaternalGfa
        File outputPaternalContigGfa = trioHifiasm.outputPaternalContigGfa
        File outputMaternalContigGfa = trioHifiasm.outputMaternalContigGfa
        File outputRawUnitigGfa = trioHifiasm.outputRawUnitigGfa
        File outputBinFiles = trioHifiasm.outputBinFiles
    }
}

task trioHifiasm {
    input{
        File paternalYak
        File maternalYak
        Array[File] childReadsHiFi
        Array[File]? childReadsUL
        Int? homCov
        String childID
	String? extraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        # runtime configurations
        Int memSizeGB
        Int threadCount
        Int diskSizeGB
        Int preemptible
        String dockerImage
        String zones
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

        ## if bin files are given we have to extract them in the directory where hifiasm is being run
        ## this enables hifiasm to skip the time-consuming process of finding read overlaps
        if [ ! -v ~{inputBinFilesTarGz} ]; then
            tar -xzf ~{inputBinFilesTarGz} --strip-components 1
            rm -rf ~{inputBinFilesTarGz}
        fi

        ## it is not possible to pipe multiple fastq files to hifiasm
        cat ~{sep=" " childReadsHiFi} > ~{childID}.fastq
        rm -rf ~{sep=" " childReadsHiFi}

        ## run trio hifiasm https://github.com/chhylp123/hifiasm
        # If ONT ultra long reads are provided
        if [[ -n "~{sep="" childReadsUL}" ]];
        then
            hifiasm ~{extraOptions} -o ~{childID} --ul ~{sep="," childReadsUL} --hom-cov ~{homCov} -t~{threadCount} -1 ~{paternalYak} -2 ~{maternalYak} ~{childID}.fastq 
        else 
            hifiasm ~{extraOptions} -o ~{childID} -t~{threadCount} -1 ~{paternalYak} -2 ~{maternalYak} ~{childID}.fastq
        fi

        #Move bin and gfa files to saparate folders and compress them 
        mkdir ~{childID}.raw_unitig_gfa
        mkdir ~{childID}.pat.contig_gfa
        mkdir ~{childID}.mat.contig_gfa
        mkdir ~{childID}.binFiles
        
        ln ~{childID}.dip.r_utg.* ~{childID}.raw_unitig_gfa
        ln *.hap1.p_ctg.* ~{childID}.pat.contig_gfa
        ln *.hap2.p_ctg.* ~{childID}.mat.contig_gfa
        ln *.bin ~{childID}.binFiles
        
        
        # make archives
        tar -cf ~{childID}.raw_unitig_gfa.tar ~{childID}.raw_unitig_gfa
        tar -cf ~{childID}.pat.contig_gfa.tar ~{childID}.pat.contig_gfa
        tar -cf ~{childID}.mat.contig_gfa.tar ~{childID}.mat.contig_gfa
        tar -cf ~{childID}.binFiles.tar ~{childID}.binFiles
        
        # compress
        pigz -p~{threadCount} ~{childID}.raw_unitig_gfa.tar
        pigz -p~{threadCount} ~{childID}.pat.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.mat.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.binFiles.tar
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        zones : zones
    }

    output {
        File outputPaternalGfa = "~{childID}.dip.hap1.p_ctg.gfa"
        File outputMaternalGfa = "~{childID}.dip.hap2.p_ctg.gfa"
        File outputPaternalContigGfa = "~{childID}.pat.contig_gfa.tar.gz"
        File outputMaternalContigGfa = "~{childID}.mat.contig_gfa.tar.gz"
        File outputRawUnitigGfa = "~{childID}.raw_unitig_gfa.tar.gz"
        File outputBinFiles = "~{childID}.binFiles.tar.gz"
    }
}

# This task does not work properly if sequences have methylation tags
task filterShortReads {
    input{
        File readFastq
        Int minReadLength
        # runtime configurations
        Int memSizeGB=8
        Int threadCount=4
        Int diskSizeGB=512
        Int preemptible=1
        String dockerImage="quay.io/masri2019/hpp_hifi_adapter_filt:latest"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir data
        cd data
        FILENAME=$(basename -- "~{readFastq}")
        PREFIX="${FILENAME%.*}"
        # filter reads shorter than minReadLength
        awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>~{minReadLength}{print a"\n"b"\n"c"\n"$0;}' ~{readFastq} > ${PREFIX}.long.fastq
        OUTPUTSIZE=`du -s -BG *.long.fastq | sed 's/G.*//'` 
        echo $OUTPUTSIZE > outputsize
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File longReadFastq = glob("data/*.long.fastq")[0]
        Int fileSizeGB = read_int("data/outputsize")
    }
}

