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

      Array[File] sampleReadsIlm

      File grch38Fasta
      File grch38InsideConfRegions
      File grch38OutsideConfRegions

      String sampleID
      String pafAligner="minimap2"

      File toilRunLog
    }

    call meryl_t.runMeryl as makeMerylDB {
        input:
          sampleReadsILM=sampleReadsIlm,
          merylCountMemSizeGB=138,
          merylUnionSumMemSizeGB=138,
          fileExtractionDiskSizeGB=512
    }

    call kmer_based_polisher_eval_wf.kmerPolishingEval as kmerPolishingEvalRaw {
        input:
          hap1Fasta=rawHap1Fasta,
          hap2Fasta=rawHap2Fasta,
          grch38Fasta=grch38Fasta,
          grch38InsideConfRegions=grch38InsideConfRegions,
          grch38OutsideConfRegions=grch38OutsideConfRegions,
          sampleID="Raw",
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
          sampleID="Polished",
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
            refAssembly=rawHap2Fasta,
            suffix="hap2PolToRaw",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    call project_blocks_t.project_blocks as projectFPKmersToRawHap1 {
        input:
            pafFile=alignHap1ToRaw.pafFile,
            bedFile=kmerPolishingEvalPolished.merquryAsmFPkmers,
            mode='asm2ref'
    }

    call project_blocks_t.project_blocks as projectFPKmersToRawHap2 {
        input:
            pafFile=alignHap2ToRaw.pafFile,
            bedFile=kmerPolishingEvalPolished.merquryAltHapFPkmers,
            mode='asm2ref'
    }

    call countEditsOverlappingFPKmers {
        input:
            polishingVcf=polishingVcf,
            hap1FPKmersProjectedBed=projectFPKmersToRawHap1.projectionBedFile,
            hap2FPKmersProjectedBed=projectFPKmersToRawHap2.projectionBedFile
    }

    call collateResults {
        input:
            editsIntersectingFPKmersTxt=countEditsOverlappingFPKmers.countsFile,
            totalEditsTxt=countEditsOverlappingFPKmers.totalEdits,

            wholeGenomeQVRawMerq=kmerPolishingEvalRaw.QV_whole_genome,
            insideConfQVRawMerq=kmerPolishingEvalRaw.QV_inside_conf,
            outsideConfQVRawMerq=kmerPolishingEvalRaw.QV_outside_conf,

            wholeGenomeQVPolMerq=kmerPolishingEvalPolished.QV_whole_genome,
            insideConfQVPolMerq=kmerPolishingEvalPolished.QV_inside_conf,
            outsideConfQVPolMerq=kmerPolishingEvalPolished.QV_outside_conf,

            yakTarBallWGRaw=kmerPolishingEvalRaw.yakTarBallWG,
            yakTarBallInsideConfRaw=kmerPolishingEvalRaw.yakTarBallInsideConf,
            yakTarBallOutsideConfRaw=kmerPolishingEvalRaw.yakTarBallOutsideConf,

            yakTarBallWGPol=kmerPolishingEvalPolished.yakTarBallWG,
            yakTarBallInsideConfPol=kmerPolishingEvalPolished.yakTarBallInsideConf,
            yakTarBallOutsideConfPol=kmerPolishingEvalPolished.yakTarBallOutsideConf,

            toilRunLog=toilRunLog,
            sampleID=sampleID
    }

    output {
      File editsIntersectingFPKmersTxt=countEditsOverlappingFPKmers.countsFile
      File totalEditsTxt=countEditsOverlappingFPKmers.totalEdits

      File hap1FPKmersProjectedBed=projectFPKmersToRawHap1.projectionBedFile
      File hap2FPKmersProjectedBed=projectFPKmersToRawHap2.projectionBedFile

      File wholeGenomeQVRawMerq=kmerPolishingEvalRaw.QV_whole_genome
      File insideConfQVRawMerq=kmerPolishingEvalRaw.QV_inside_conf
      File outsideConfQVRawMerq=kmerPolishingEvalRaw.QV_outside_conf

      File wholeGenomeQVPolMerq=kmerPolishingEvalPolished.QV_whole_genome
      File insideConfQVPolMerq=kmerPolishingEvalPolished.QV_inside_conf
      File outsideConfQVPolMerq=kmerPolishingEvalPolished.QV_outside_conf

      File yakTarBallWGRaw=kmerPolishingEvalRaw.yakTarBallWG
      File yakTarBallInsideConfRaw=kmerPolishingEvalRaw.yakTarBallInsideConf
      File yakTarBallOutsideConfRaw=kmerPolishingEvalRaw.yakTarBallOutsideConf

      File yakTarBallWGPol=kmerPolishingEvalPolished.yakTarBallWG
      File yakTarBallInsideConfPol=kmerPolishingEvalPolished.yakTarBallInsideConf
      File yakTarBallOutsideConfPol=kmerPolishingEvalPolished.yakTarBallOutsideConf

      File collatedQCResults=collateResults.QC_stats

      File hap1ToRawPaf=alignHap1ToRaw.pafFile
      File hap2ToRawPaf=alignHap2ToRaw.pafFile
      File FPkmersHap1=kmerPolishingEvalPolished.merquryAsmFPkmers
      File FPkmersHap2=kmerPolishingEvalPolished.merquryAltHapFPkmers
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

        gunzip -c ~{polishingVcf} > polisher_output.vcf

        cat ~{hap1FPKmersProjectedBed} ~{hap2FPKmersProjectedBed} | sed 's/[\t]*$//' > dip.FPkmers.projected
        bedtools intersect -a polisher_output.vcf -b dip.FPkmers.projected | sort | uniq | wc -l > edits_intersecting_FPkmers.txt
        grep -v "^#" polisher_output.vcf | sort | uniq | wc -l > total_edits.txt
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
      File outsideConfQVRawMerq

      File wholeGenomeQVPolMerq
      File insideConfQVPolMerq
      File outsideConfQVPolMerq

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
        echo "sampleID,Assembly,Runtime,total_edits,edits_overlapping_FPkmers,WholeGenomeQV_Merqury_Hap1,WholeGenomeQV_Merqury_Hap2,WholeGenomeQV_Merqury_Dip,WholeGenomeQV_Yak_Hap1,WholeGenomeQV_Yak_Hap2,WholeGenomeQV_Yak_Dip,WholeGenomeQV_Yak_Hap1_unNormalized,WholeGenomeQV_Yak_Hap2_unNormalized,WholeGenomeQV_Yak_Dip_unNormalized,InsideConfQV_Merqury_Hap1,InsideConfQV_Merqury_Hap2,InsideConfQV_Merqury_Dip,InsideConfQV_Yak_Hap1,InsideConfQV_Yak_Hap2,InsideConfQV_Yak_Dip,InsideConfQV_Yak_Hap1_unNormalized,InsideConfQV_Yak_Hap2_unNormalized,InsideConfQV_Yak_Dip_unNormalized,OutsideConfQV_Merqury_Hap1,OutsideConfQV_Merqury_Hap2,OutsideConfQV_Merqury_Dip,OutsideConfQV_Yak_Hap1,OutsideConfQV_Yak_Hap2,OutsideConfQV_Yak_Dip_unNormalized,OutsideConfQV_Yak_Hap1_unNormalized,OutsideConfQV_Yak_Hap2_unNormalized,OutsideConfQV_Yak_Dip_unNormalized,YakSwitchError,YakHammingError" > header.csv

        # add sample ID
        echo ~{sampleID},"polished" >> polished.sample.csv
        echo ~{sampleID},"raw","NA","NA","NA" >> raw.sample.csv

        ### Collate polished assembly results ###

        # get runtime
        grep "real" ~{toilRunLog} | cut -f2 > runtime.txt

        # get merqury WG results
        grep "asm" ~{wholeGenomeQVPolMerq} | cut -f4 > wholeGenomeQVPolMerq.Hap1.txt
        grep "altHap" ~{wholeGenomeQVPolMerq} | cut -f4 > wholeGenomeQVPolMerq.Hap2.txt
        grep "Both" ~{wholeGenomeQVPolMerq} | cut -f4 > wholeGenomeQVPolMerq.dip.txt

        # get Yak WG results
        mkdir yak_WG_polished
        tar -C yak_WG_polished -zxvf ~{yakTarBallWGPol}
        head -n 5 yak_WG_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_WG_pol.hap2.txt
        head -n 10 yak_WG_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_WG_pol.hap1.txt
        head -n 15 yak_WG_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_WG_pol.dip.txt

        # unnormalized QV
        head -n 5 yak_WG_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_WG_pol_unNorm.hap2.txt
        head -n 10 yak_WG_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_WG_pol_unNorm.hap1.txt
        head -n 15 yak_WG_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_WG_pol_unNorm.dip.txt

        # Get Merqury inside confidence regions results
        grep "asm" ~{insideConfQVPolMerq} | cut -f4 > insideConfQVPolMerq.Hap1.txt
        grep "altHap" ~{insideConfQVPolMerq} | cut -f4 > insideConfQVPolMerq.Hap2.txt
        grep "Both" ~{insideConfQVPolMerq} | cut -f4 > insideConfQVPolMerq.dip.txt

        # Get Yak inside confidence regions results
        mkdir yak_inside_polished
        tar -C yak_inside_polished -zxvf ~{yakTarBallInsideConfPol}
        head -n 5 yak_inside_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_inside_pol.hap2.txt
        head -n 10 yak_inside_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_inside_pol.hap1.txt
        head -n 15 yak_inside_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_inside_pol.dip.txt

        head -n 5 yak_inside_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_inside_pol_unNorm.hap2.txt
        head -n 10 yak_inside_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_inside_pol_unNorm.hap1.txt
        head -n 15 yak_inside_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_inside_pol_unNorm.dip.txt

        # Get Merqury outside confidence regions results
        grep "asm" ~{outsideConfQVPolMerq} | cut -f4 > outsideConfQVPolMerq.Hap1.txt
        grep "altHap" ~{outsideConfQVPolMerq} | cut -f4 > outsideConfQVPolMerq.Hap2.txt
        grep "Both" ~{outsideConfQVPolMerq} | cut -f4 > outsideConfQVPolMerq.dip.txt

        # Get Yak outside confidence regions results
        mkdir yak_outside_polished
        tar -C yak_outside_polished -zxvf ~{yakTarBallOutsideConfPol}
        head -n 5 yak_outside_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_outside_pol.hap2.txt
        head -n 10 yak_outside_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_outside_pol.hap1.txt
        head -n 15 yak_outside_polished/*.summary.txt | tail -n 1 | cut -f 3 > yak_outside_pol.dip.txt

        head -n 5 yak_outside_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_outside_pol_unNorm.hap2.txt
        head -n 10 yak_outside_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_outside_pol_unNorm.hap1.txt
        head -n 15 yak_outside_polished/*.summary.txt | tail -n 1 | cut -f 2 > yak_outside_pol_unNorm.dip.txt

        # Get Yak Switch / Hamming
        tail -n 3 yak_WG_polished/*.summary.txt | head -n 1 | cut -f4 > yak_switch.txt
        tail -n 2 yak_WG_polished/*.summary.txt | head -n 1 | cut -f4 > yak_hamming.txt

        # Paste polished results into one row
        paste -d "," polished.sample.csv \
        runtime.txt \
        ~{totalEditsTxt} \
        ~{editsIntersectingFPKmersTxt} \
        wholeGenomeQVPolMerq.Hap1.txt \
        wholeGenomeQVPolMerq.Hap2.txt \
        wholeGenomeQVPolMerq.dip.txt \
        yak_WG_pol.hap1.txt \
        yak_WG_pol.hap2.txt \
        yak_WG_pol.dip.txt \
        yak_WG_pol_unNorm.hap1.txt \
        yak_WG_pol_unNorm.hap2.txt \
        yak_WG_pol_unNorm.dip.txt \
        insideConfQVPolMerq.Hap1.txt \
        insideConfQVPolMerq.Hap2.txt \
        insideConfQVPolMerq.dip.txt \
        yak_inside_pol.hap1.txt \
        yak_inside_pol.hap2.txt \
        yak_inside_pol.dip.txt \
        yak_inside_pol_unNorm.hap1.txt \
        yak_inside_pol_unNorm.hap2.txt \
        yak_inside_pol_unNorm.dip.txt \
        outsideConfQVPolMerq.Hap1.txt \
        outsideConfQVPolMerq.Hap2.txt \
        outsideConfQVPolMerq.dip.txt \
        yak_outside_pol.hap1.txt \
        yak_outside_pol.hap2.txt \
        yak_outside_pol.dip.txt \
        yak_outside_pol_unNorm.hap1.txt \
        yak_outside_pol_unNorm.hap2.txt \
        yak_outside_pol_unNorm.dip.txt \
        yak_switch.txt \
        yak_hamming.txt \
        > tmp ; mv tmp polished.sample.csv

        ### Collate raw assembly results ###
        # get merqury WG results
        grep "asm" ~{wholeGenomeQVRawMerq} | cut -f4 > wholeGenomeQVRawMerq.Hap1.txt
        grep "altHap" ~{wholeGenomeQVRawMerq} | cut -f4 > wholeGenomeQVRawMerq.Hap2.txt
        grep "Both" ~{wholeGenomeQVRawMerq} | cut -f4 > wholeGenomeQVRawMerq.dip.txt

        # get Yak WG results
        mkdir yak_WG_Raw
        tar -C yak_WG_Raw -zxvf ~{yakTarBallWGRaw}
        head -n 5 yak_WG_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_WG_Raw.hap2.txt
        head -n 10 yak_WG_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_WG_Raw.hap1.txt
        head -n 15 yak_WG_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_WG_Raw.dip.txt

        head -n 5 yak_WG_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_WG_Raw_unNorm.hap2.txt
        head -n 10 yak_WG_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_WG_Raw_unNorm.hap1.txt
        head -n 15 yak_WG_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_WG_Raw_unNorm.dip.txt

        # Get Merqury inside confidence regions results
        grep "asm" ~{insideConfQVRawMerq} | cut -f4 > insideConfQVRawMerq.Hap1.txt
        grep "altHap" ~{insideConfQVRawMerq} | cut -f4 > insideConfQVRawMerq.Hap2.txt
        grep "Both" ~{insideConfQVRawMerq} | cut -f4 > insideConfQVRawMerq.dip.txt

        # Get Yak inside confidence regions results
        mkdir yak_inside_Raw
        tar -C yak_inside_Raw -zxvf ~{yakTarBallInsideConfRaw}
        head -n 5 yak_inside_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_inside_Raw.hap2.txt
        head -n 10 yak_inside_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_inside_Raw.hap1.txt
        head -n 15 yak_inside_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_inside_Raw.dip.txt

        head -n 5 yak_inside_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_inside_Raw_unNorm.hap2.txt
        head -n 10 yak_inside_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_inside_Raw_unNorm.hap1.txt
        head -n 15 yak_inside_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_inside_Raw_unNorm.dip.txt

        # Get Merqury outside confidence regions results
        grep "asm" ~{outsideConfQVRawMerq} | cut -f4 > outsideConfQVRawMerq.Hap1.txt
        grep "altHap" ~{outsideConfQVRawMerq} | cut -f4 > outsideConfQVRawMerq.Hap2.txt
        grep "Both" ~{outsideConfQVRawMerq} | cut -f4 > outsideConfQVRawMerq.dip.txt

        # Get Yak outside confidence regions results
        mkdir yak_outside_Raw
        tar -C yak_outside_Raw -zxvf ~{yakTarBallOutsideConfRaw}
        head -n 5 yak_outside_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_outside_Raw.hap2.txt
        head -n 10 yak_outside_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_outside_Raw.hap1.txt
        head -n 15 yak_outside_Raw/*.summary.txt | tail -n 1 | cut -f 3 > yak_outside_Raw.dip.txt

        head -n 5 yak_outside_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_outside_Raw_unNorm.hap2.txt
        head -n 10 yak_outside_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_outside_Raw_unNorm.hap1.txt
        head -n 15 yak_outside_Raw/*.summary.txt | tail -n 1 | cut -f 2 > yak_outside_Raw_unNorm.dip.txt

        # Get Yak Switch / Hamming
        tail -n 3 yak_WG_Raw/*.summary.txt | head -n 1 | cut -f4 > yak_switch.txt
        tail -n 2 yak_WG_Raw/*.summary.txt | head -n 1 | cut -f4 > yak_hamming.txt

        # Paste raw results into one row

        paste -d "," raw.sample.csv \
        wholeGenomeQVRawMerq.Hap1.txt \
        wholeGenomeQVRawMerq.Hap2.txt \
        wholeGenomeQVRawMerq.dip.txt \
        yak_WG_Raw.hap1.txt \
        yak_WG_Raw.hap2.txt \
        yak_WG_Raw.dip.txt \
        yak_WG_Raw_unNorm.hap1.txt \
        yak_WG_Raw_unNorm.hap2.txt \
        yak_WG_Raw_unNorm.dip.txt \
        insideConfQVRawMerq.Hap1.txt \
        insideConfQVRawMerq.Hap2.txt \
        insideConfQVRawMerq.dip.txt \
        yak_inside_Raw.hap1.txt \
        yak_inside_Raw.hap2.txt \
        yak_inside_Raw.dip.txt \
        yak_inside_Raw_unNorm.hap1.txt \
        yak_inside_Raw_unNorm.hap2.txt \
        yak_inside_Raw_unNorm.dip.txt \
        outsideConfQVRawMerq.Hap1.txt \
        outsideConfQVRawMerq.Hap2.txt \
        outsideConfQVRawMerq.dip.txt \
        yak_outside_Raw.hap1.txt \
        yak_outside_Raw.hap2.txt \
        yak_outside_Raw.dip.txt \
        yak_outside_Raw_unNorm.hap1.txt \
        yak_outside_Raw_unNorm.hap2.txt \
        yak_outside_Raw_unNorm.dip.txt \
        yak_switch.txt \
        yak_hamming.txt \
        > tmp ; mv tmp raw.sample.csv

        # combine results with header
        cat header.csv raw.sample.csv polished.sample.csv > ~{sampleID}.polishing.QC.csv
  >>>

  output {
      File QC_stats = glob("*polishing.QC.csv")[0]
  }

  runtime{
      memory: memSizeGB + " GB"
      cpu: threadCount
      disks: "local-disk " + diskSizeGB + " SSD"
      docker: dockerImage
  }
}
