version 1.0

import "./kmer_based_polisher_eval.wdl" as kmer_based_polisher_eval_wf
import "../tasks/project_blocks.wdl" as project_blocks_t
import "../tasks/meryl_non_trio.wdl" as meryl_t
import "../tasks/long_read_aligner.wdl" as long_read_aligner_t

workflow hprc_polishing_QC {

    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Runs kmer based QC for hprc assemblies, before and after polishing"
      }

    input {
      File rawHap1Fasta
      File rawHap2Fasta
      File polishedHap1Fasta
      File polishedHap2Fasta
      File polishingVcf

      File sampleYak
      File maternalYak
      File paternalYak

      Array[File] sampleIlmFastq

      File grch38Fasta
      File grch38InsideConfRegions
      File grch38OutsideConfRegions

      String polSampleID
      String rawSampleID
      String pafAligner="minimap2"
    }
    call meryl_t.runMeryl as makeMerylDB {
        input:
          sampleReadsILM=sampleIlmFastq
    }

    call kmer_based_polisher_eval_wf.kmerPolishingEval as kmerPolishingEvalRaw {
        input:
          hap1Fasta=rawHap1Fasta,
          hap2Fasta=rawHap2Fasta,
          grch38Fasta=grch38Fasta,
          grch38InsideConfRegions=grch38InsideConfRegions,
          grch38OutsideConfRegions=grch38OutsideConfRegions,
          sampleID=rawSampleID,
          ilmMerylDBTarGz=makeMerylDB.sampleMerylDB,
          sampleYak=sampleYak,
          paternalYak=paternalYak,
          maternalYak=maternalYak
    }

    call kmer_based_polisher_eval_wf.kmerPolishingEval as kmerPolishingEvalPolished {
        input:
          hap1Fasta=polishedHap1Fasta,
          hap2Fasta=polishedHap2Fasta,
          grch38Fasta=grch38Fasta,
          grch38InsideConfRegions=grch38InsideConfRegions,
          grch38OutsideConfRegions=grch38OutsideConfRegions,
          sampleID=polSampleID,
          ilmMerylDBTarGz=makeMerylDB.sampleMerylDB,
          sampleYak=sampleYak,
          paternalYak=paternalYak,
          maternalYak=maternalYak
    }

    # Align hap1 and hap2 polished to hap1 and hap2 raw
    call long_read_aligner_t.alignmentPaf as alignHap1ToRaw{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=polishedHap1Fasta,
            refAssembly=rawHap1Fasta,
            suffix="hap1PolToRaw",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }
    call long_read_aligner_t.alignmentPaf as alignHap2ToRaw{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=polishedHap2Fasta,
            refAssembly=rawHap1Fasta,
            suffix="hap2PolToRaw",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    call project_blocks_t.project_blocks as projectFPKmersToRawHap1 {
        input:
            pafFile=alignHap1ToRaw.pafFile,
            bedFile=kmerPolishingEvalPolished.merquryAsmFPkmers
    }

    call project_blocks_t.project_blocks as projectFPKmersToRawHap2 {
        input:
            pafFile=alignHap2ToRaw.pafFile,
            bedFile=kmerPolishingEvalPolished.merquryAltHapFPkmers
    }

    call countEditsOverlappingFPKmers {
        input:
            polishingVcf=polishingVcf,
            hap1FPKmersProjectedBed=projectFPKmersToRawHap1.projectionBedFile,
            hap2FPKmersProjectedBed=projectFPKmersToRawHap2.projectionBedFile
    }

    output {
      File editsIntersectingFPKmersTxt=countEditsOverlappingFPKmers.countsFile
      File totalEditsTxt=countEditsOverlappingFPKmers.totalEdits
      
      File hap1FPKmersProjectedBed=projectFPKmersToRawHap1.projectionBedFile
      File hap2FPKmersProjectedBed=projectFPKmersToRawHap2.projectionBedFile

      File wholeGenomeQVRawMerq=kmerPolishingEvalRaw.QV_whole_genome
      File insideConfQVRawMerq=kmerPolishingEvalRaw.QV_inside_conf
      File outsideConfRawMerq=kmerPolishingEvalRaw.QV_outside_conf

      File wholeGenomeQVPolMerq=kmerPolishingEvalPolished.QV_whole_genome
      File insideConfQVPolMerq=kmerPolishingEvalPolished.QV_inside_conf
      File outsideConfPolMerq=kmerPolishingEvalPolished.QV_outside_conf

      File yakTarBallWGRaw=kmerPolishingEvalRaw.yakTarBallWG
      File yakTarBallInsideConfRaw=kmerPolishingEvalRaw.yakTarBallInsideConf
      File yakTarBallOutsideConfRaw=kmerPolishingEvalRaw.yakTarBallOutsideConf

      File yakTarBallWGPol=kmerPolishingEvalPolished.yakTarBallWG
      File yakTarBallInsideConfPol=kmerPolishingEvalPolished.yakTarBallInsideConf
      File yakTarBallOutsideConfPol=kmerPolishingEvalPolished.yakTarBallOutsideConf
    }
}

task countEditsOverlappingFPKmers {
  input {
      File polishingVcf
      File hap1FPKmersProjectedBed
      File hap2FPKmersProjectedBed
      Int memSizeGB = 12
      Int threadCount = 4
      Int diskSizeGB = 64
      String dockerImage = "mobinasri/flagger:latest"
  }

  command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cat ~{hap1FPKmersProjectedBed} ~{hap2FPKmersProjectedBed} > dip.FPkmers.projected
        bedtools intersect -a ~{polishingVcf} -b dip.FPkmers.projected | sort | uniq | wc -l > edits_intersecting_FPkmers.txt
        gunzip ~{polishingVcf} | grep -v "^#" | sort | uniq | wc -l > total_edits.txt
  >>>

  output {
      File countsFile = "edits_intersecting_FPkmers.txt"
      File totalEdits="total_edits.txt"
  }

  runtime{
      memory: memSizeGB + " GB"
      cpu: threadCount
      disks: "local-disk " + diskSizeGB + " SSD"
      docker: dockerImage
  }
}

task collateResults {
  input {
      File editsIntersectingFPKmersTxt

      File wholeGenomeQVRawMerq
      File insideConfQVRawMerq
      File outsideConfRawMerq

      File wholeGenomeQVPolMerq
      File insideConfQVPolMerq
      File outsideConfPolMerq

      File yakTarBallWGRaw
      File yakTarBallInsideConfRaw
      File yakTarBallOutsideConfRaw

      File yakTarBallWGPol
      File yakTarBallInsideConfPol
      File yakTarBallOutsideConfPol

      String sampleID

      Int memSizeGB = 12
      Int threadCount = 4
      Int diskSizeGB = 64
      String dockerImage = "mobinasri/flagger:latest"
  }

  command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        echo "sampleID,edits_overlapFPkmers,total_edits,WholeGenomeQV_Merqury_Hap1,WholeGenomeQV_," > ~{sampleID}.polishing.QC.tsv
  >>>

  output {
      File QC_stats = glob("*polishing.QC.tsv")[0]
  }

  runtime{
      memory: memSizeGB + " GB"
      cpu: threadCount
      disks: "local-disk " + diskSizeGB + " SSD"
      docker: dockerImage
  }
}
