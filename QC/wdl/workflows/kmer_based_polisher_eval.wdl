version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/merqury.wdl" as merqury_t
import "../tasks/yak.wdl" as yak_t
import "../tasks/project_blocks.wdl" as project_blocks_t
import "../tasks/subFastaByBed.wdl" as subset_fasta_t

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
      Array[File] maternalReadsILM
      Array[File] paternalReadsILM
      Array[File] sampleReadsILM
      }

    # Align hap1 and hap2 to grch38, in paf format
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
    call long_read_aligner_t.alignmentPaf as alignHap2ToRef{
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
    # project inside and outside confidence regions to Hap1 and Hap2 coordinates
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

    call subset_fasta_t.SubFastaByBed as subHap1InsideConf {
        input:
            Fasta=hap1Fasta,
            Bed=projectInsideConfHap1.projectionBedFile
    }
    call subset_fasta_t.SubFastaByBed as subHap2InsideConf {
        input:
            Fasta=hap2Fasta,
            Bed=projectInsideConfHap2.projectionBedFile
    }
    call subset_fasta_t.SubFastaByBed as subHap1OutsideConf {
        input:
            Fasta=hap1Fasta,
            Bed=projectOutsideConfHap1.projectionBedFile
    }
    call subset_fasta_t.SubFastaByBed as subHap2OutsideConf {
        input:
            Fasta=hap2Fasta,
            Bed=projectOutsideConfHap2.projectionBedFile
    }

    # Run merqury QV whole genome
    call merqury_t.merqury as merquryWholeGenome {
        input:
            assemblyFasta=hap1Fasta,
            altHapFasta=hap2Fasta,
            kmerTarball=ilmMerylDBTarGz
    }
    call merqury_t.merqury as merquryInsideConf {
        input:
            assemblyFasta=subHap1InsideConf.subFasta,
            altHapFasta=subHap2InsideConf.subFasta,
            kmerTarball=ilmMerylDBTarGz
    }
    call merqury_t.merqury as merquryOutsideConf {
        input:
            assemblyFasta=subHap1OutsideConf.subFasta,
            altHapFasta=subHap2OutsideConf.subFasta,
            kmerTarball=ilmMerylDBTarGz
    }

    # Run Yak QV whole genome, inside and outside conf
    call yak_t.runYakAssemblyStats as yakQCWholeGenome {
        input:
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            sampleReadsILM=sampleReadsILM,
            assemblyFastaPat=hap1Fasta,
            assemblyFastaMat=hap2Fasta
    }

    call yak_t.runYakAssemblyStats as yakQCInsideConf {
        input:
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            sampleReadsILM=sampleReadsILM,
            assemblyFastaPat=subHap1InsideConf.subFasta,
            assemblyFastaMat=subHap2InsideConf.subFasta
    }
    call yak_t.runYakAssemblyStats as yakQCOutsideConf {
        input:
            maternalReadsILM=maternalReadsILM,
            paternalReadsILM=paternalReadsILM,
            sampleReadsILM=sampleReadsILM,
            assemblyFastaPat=subHap1OutsideConf.subFasta,
            assemblyFastaMat=subHap2OutsideConf.subFasta
    }

    output {
        File QV_whole_genome = merquryWholeGenome.QV
        File QV_inside_conf = merquryInsideConf.QV
        File QV_outside_conf = merquryOutsideConf.QV
        File yakTarBallWG=yakQCWholeGenome.outputTarball
        File yakTarBallInsideConf=yakQCInsideConf.outputTarball
        File yakTarBallOutsideConf=yakQCOutsideConf.outputTarball
    }
}
