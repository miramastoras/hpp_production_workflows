version 1.0

import "../tasks/merqury.wdl" as merqury_t
import "../tasks/yak.wdl" as yak_t
import "../tasks/project_blocks.wdl" as project_blocks_t

workflow kmerPolishingEval {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Evaluate assemblies before and after polishing with kmer-based metrics (QV,switch error,hamming error) "
      }

    input {
      File rawHap1Fasta
      File rawHap2Fasta
      File polishedHap1Fasta
      File polishedHap2Fasta

      File grch38Fasta
      File grch38InsideConfRegions
      File grch38OutsideConfRegions

      File ilmMerylDBTarGz
      }

    call project_blocks_t as projectInsideConfRawHap1 {
        input:
            assemblyFasta=rawHap1Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t as projectInsideConfRawHap2 {
        input:
            assemblyFasta=rawHap2Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t as projectOutsideConfRawHap1 {
        input:
            assemblyFasta=rawHap1Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38OutsideConfRegions
    }
    call project_blocks_t as projectOutsideConfRawHap2 {
        input:
            assemblyFasta=rawHap2Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38OutsideConfRegions
    }
    call project_blocks_t as projectInsideConfPolishHap1 {
        input:
            assemblyFasta=polishedHap1Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t as projectInsideConfPolishHap2 {
        input:
            assemblyFasta=polishedHap2Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t as projectOutsideConfPolishHap1 {
        input:
            assemblyFasta=polishedHap1Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38OutsideConfRegions
    }
    call project_blocks_t as projectOutsideConfPolishHap2 {
        input:
            assemblyFasta=polishedHap2Fasta,
            refFasta=grch38Fasta,
            bedFile=grch38OutsideConfRegions
    }
    call subset_fastas {

    }
    call

    output {
        File QV = merqury.QV
        File outputTarball = merqury.outputTarball
        File FPkmers = merqury.FPkmers
    }
}

# subset fasta to inside and outside conf regions
task subset_fastas {

}



# project inside and outside GIAB conf regions to raw and polished assembly
# subset raw and polished fasta to inside and outside conf regions
# run merqury and yak Qv on all 6 fasta files
# run yak switch and hamming on whole genome fasta files
# combine results into one csv file
