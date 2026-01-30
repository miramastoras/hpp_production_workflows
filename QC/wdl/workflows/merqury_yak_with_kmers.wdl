version 1.0

import "../tasks/merqury.wdl" as merqury_t
import "../tasks/yak.wdl" as yak_t
import "../tasks/yak_non_trio.wdl" as yak_non_trio_t
import "../tasks/yak_meryl_count.wdl" as yak_meryl_count_t

workflow runMerquryAndYak {

    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Runs meryl, yak count, merqury, and yak for a set kmer size, illumina reads only"
      }

    input {
      File Hap1Fasta
      File Hap2Fasta

      File maternalYak
      File paternalYak

      Array[File] sampleReadsIlm

      String sampleID

      Int yakMerylKmerSize=31
    }

    call yak_meryl_count_t.runYakMerylCount as countYakMerylKmers {
        input:
            sampleReadsIlm=sampleReadsIlm,
            kmerSize=yakMerylKmerSize,
            threadCount=32,
            sampleID=sampleID
    }

    # Run merqury QV whole genome
    call merqury_t.merqury as merquryWholeGenome {
        input:
            assemblyFasta=Hap1Fasta,
            altHapFasta=Hap2Fasta,
            kmerTarball=countYakMerylKmers.merylDbTarGz
    }

    call yak_t.yakAssemblyStats as yakQCWholeGenome {
        input:
            matYak=maternalYak,
            patYak=paternalYak,
            sampleYak=countYakMerylKmers.sampleYak,
            assemblyFastaPat=Hap1Fasta,
            assemblyFastaMat=Hap2Fasta,
            minSequenceLength="0",
            dockerImage="miramastoras/hpp_yak:latest"
    }

    output {
      File QV_whole_genome = merquryWholeGenome.QV
      File merquryWGTarBall = merquryWholeGenome.outputTarball
      File yakSummary=yakQCWholeGenome.outputSummary
      File yakTarBallWG = yakQCWholeGenome.outputTarball
    }
}
