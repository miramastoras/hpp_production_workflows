version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/merqury.wdl" as merqury_t
import "../tasks/yak.wdl" as yak_t
import "../tasks/project_blocks.wdl" as project_blocks_t

workflow kmerPolishingEval {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Evaluate diploid assembly with kmer-based metrics (QV,switch error,hamming error)"
      }

    input {
      File hap1Fasta
      File hap2Fasta

      File grch38Fasta
      File grch38InsideConfRegions
      File grch38OutsideConfRegions

      String sampleID
      String pafAligner="minimap2"
      String mode="ref2asm"

      File ilmMerylDBTarGz
      }

    call long_read_aligner_t.alignmentPaf as alignHap1ToRef{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=hap1Fasta,
            refAssembly=grch38Fasta,
            suffix="asmToRef",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }
    call long_read_aligner_t.alignmentPaf as alignHap1ToRef{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=hap2Fasta,
            refAssembly=grch38Fasta,
            suffix="asmToRef",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    call project_blocks_t.project_blocks as projectInsideConfHap1 {
        input:
            pafFile=alignHap1ToRef.pafFile,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t.project_blocks as projectInsideConfHap2 {
        input:
            pafFile=alignHap2ToRef.pafFile,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t.project_blocks as projectOutsideConfHap1 {
        input:
            pafFile=alignHap1ToRef.pafFile,
            bedFile=grch38OutsideConfRegions
    }
    call project_blocks_t.project_blocks as projectOutsideConfHap2 {
        input:
            pafFile=alignHap2ToRef.pafFile,
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
