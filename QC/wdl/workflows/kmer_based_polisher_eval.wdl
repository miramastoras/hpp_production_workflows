version 1.0

import "../tasks/merqury.wdl" as merqury_t
import "../tasks/yak.wdl" as yak_t

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


    }
    call runMeryl
    output {
        File QV = merqury.QV
        File outputTarball = merqury.outputTarball
        File FPkmers = merqury.FPkmers
    }
}

# project inside and outside GIAB conf regions to raw and polished assembly
# subset raw and polished fasta to inside and outside conf regions
# run merqury and yak Qv on all 6 fasta files
# run yak switch and hamming on whole genome fasta files
# combine results into one csv file
