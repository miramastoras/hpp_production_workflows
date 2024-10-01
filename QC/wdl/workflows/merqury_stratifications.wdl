version 1.0

import "../tasks/subFastaByBed.wdl" as subset_fasta_t
import "../tasks/merqury.wdl" as merqury_t

workflow merqury_stratifications {

    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Subsets an assembly and runs merqury whole genome, in and out of provided bed regions."
      }

    input {
      File Hap1Fasta
      File Hap2Fasta
      File dipFai
      File bedFile
      File ilmMerylDBTarGz
      String sampleID

    }

    call getInverseBed {
        input:
            dipFai=dipFai
            bedFile=bedFile
    }
    call subset_fasta_t.SubFastaByBed as subHap1InsideBed {
        input:
            Fasta=Hap1Fasta,
            Bed=bedFile,
            outputLabel="hap1.insideBed",
            sampleID=sampleID
    }
    call subset_fasta_t.SubFastaByBed as subHap2InsideBed {
        input:
            Fasta=Hap2Fasta,
            Bed=bedFile,
            outputLabel="hap2.insideBed",
            sampleID=sampleID
    }
    call subset_fasta_t.SubFastaByBed as subHap1OutsideBed {
        input:
            Fasta=Hap1Fasta,
            Bed=getInverseBed.inverseBed,
            outputLabel="hap1.outsideBed",
            sampleID=sampleID
    }
    call subset_fasta_t.SubFastaByBed as subHap2OutsideBed {
        input:
            Fasta=Hap2Fasta,
            Bed=getInverseBed.inverseBed,
            outputLabel="hap2.outsideBed",
            sampleID=sampleID
    }
    call merqury_t.merqury as merquryInsideBed {
        input:
            assemblyFasta=subHap1InsideBed.subFasta,
            altHapFasta=subHap2InsideBed.subFasta,
            kmerTarball=ilmMerylDBTarGz
    }
    call merqury_t.merqury as merquryOutsideBed {
        input:
            assemblyFasta=subHap1OutsideBed.subFasta,
            altHapFasta=subHap2OutsideBed.subFasta,
            kmerTarball=ilmMerylDBTarGz
    }
    output {
        File insideQV = merquryInsideBed.QV
        File insideMerquryTarball = merquryInsideBed.outputTarball
        File outsideQV = merquryOutsideBed.QV
        File outsideMerquryTarball = merquryOutsideBed.outputTarball
    }
}

task getInverseBed {
  input {
      File dipFai
      File bedFile

      Int memSizeGB = 12
      Int threadCount = 4
      Int diskSizeGB = 64
      String dockerImage = "pegi3s/bedtools"
  }

  command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        awk '{print $1"\t""0""\t"$2}' ~{dipFai} > dip.fai.bed

        FILENAME=$(basename -- "~{bedFile}")

        bedtools subtract -a dip.fai.bed -b ~{bedFile} > ${FILENAME}_inverse.bed
  >>>

  output {
      File inverseBed = glob("*inverse.bed")[0]
  }

  runtime{
      memory: memSizeGB + " GB"
      cpu: threadCount
      disks: "local-disk " + diskSizeGB + " SSD"
      docker: dockerImage
  }

}
