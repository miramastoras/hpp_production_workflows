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

      File toilRunLog
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
        gunzip -c ~{polishingVcf} | grep -v "^#" | sort | uniq | wc -l > total_edits.txt
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
      File totalEditsTxt

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

      File toilRunLog

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

        # define header
        echo "sampleID,Runtime,total_edits,edits_overlapping_FPkmers,WholeGenomeQV_Merqury_Hap1,WholeGenomeQV_Merqury_Hap2,WholeGenomeQV_Merqury_Dip,WholeGenomeQV_Yak_Hap1,WholeGenomeQV_Yak_Hap2,WholeGenomeQV_Yak_Dip,InsideConfQV_Merqury_Hap1,InsideConfQV_Merqury_Hap2,InsideConfQV_Merqury_Dip,InsideConfQV_Yak_Hap1,InsideConfQV_Yak_Hap2,InsideConfQV_Yak_Dip,OutsideConfQV_Merqury_Hap1,OutsideConfQV_Merqury_Hap2,OutsideConfQV_Merqury_Dip,OutsideConfQV_Yak_Hap1,OutsideConfQV_Yak_Hap2,OutsideConfQV_Yak_Dip" > header.csv

        # add sample ID
        echo ~{sampleID} >> sample.csv

        # add total edits and total edits overlapping FP kmers
        paste -d "," sample.csv ~{totalEditsTxt} ~{editsIntersectingFPKmersTxt} > tmp ; mv tmp sample.csv

        # add runtime
        

        # combine results with header
        cat header.csv sample.csv > ~{sampleID}.polishing.QC.csv

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
