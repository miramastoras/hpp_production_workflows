version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/project_blocks.wdl" as project_blocks_t
import "./merqury_stratifications.wdl" as merqury_stratifications_t
import "./align_asm_project_blocks.wdl" as align_asm_project_blocks_t

workflow merqury_GIAB_callable {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "For a polished diploid assembly, run Merqury whole genome and within GIAB confidence mappable regions"
    }

    input {
        File polishedHap1Fasta
        File polishedHap2Fasta

        File rawHap1Fasta
        File rawHap2Fasta

        File refFasta

        String sampleID
        String stratificationLabel="GIAB_conf_gt5x_MAPQ1"

        File dipFai
        File mosdepthBed
        File GIABConfBed
        File ilmMerylDBTarGz

      }

      # project GIAB confidence regions to polished assembly
      call align_asm_project_blocks_t.align_asm_project_blocks as projectGIABConfToAsm {
          input:
              asmHap1Fasta=polishedHap1Fasta,
              asmHap2Fasta=polishedHap2Fasta,
              refHap1Fasta=refFasta,
              refHap2Fasta=refFasta,
              projectionDirection="ref2asm",
              bedFile=GIABConfBed,
              sampleID="${sampleID}.GRCh38"
      }
      # project mosdepth bed file from raw to polished coordinates
      call align_asm_project_blocks_t.align_asm_project_blocks as projectMosdepthToPolishedAsm {
          input:
              asmHap1Fasta=polishedHap1Fasta,
              asmHap2Fasta=polishedHap2Fasta,
              refHap1Fasta=rawHap1Fasta,
              refHap2Fasta=rawHap2Fasta,
              projectionDirection="ref2asm",
              bedFile=mosdepthBed,
              sampleID="${sampleID}.polished2raw"
      }

      # subtract mosdepth bed from GIAB conf bed
      call bedtoolsSubtract as Hap1PolishedSubtractMosdepth {
          input:
              bedFileA=projectGIABConfToAsm.projectionBedFileHap1,
              bedFileB=projectMosdepthToPolishedAsm.projectionBedFileHap1,
              label="GIAB_conf_gt5x_MAPQ1_HAP1",
              sampleID=sampleID
      }
      call bedtoolsSubtract as Hap2PolishedSubtractMosdepth {
          input:
              bedFileA=projectGIABConfToAsm.projectionBedFileHap2,
              bedFileB=projectMosdepthToPolishedAsm.projectionBedFileHap2,
              label="GIAB_conf_gt5x_MAPQ1_HAP2",
              sampleID=sampleID
      }
      call concatBeds as concatHapBeds {
          input:
              bedFile1=Hap1PolishedSubtractMosdepth.subtractedBed,
              bedFile2=Hap2PolishedSubtractMosdepth.subtractedBed,
              Prefix=stratificationLabel,
              sampleID=sampleID
      }
      call merqury_stratifications_t.merqury_stratifications as runMerqury {
          input:
              Hap1Fasta=polishedHap1Fasta,
              Hap2Fasta=polishedHap2Fasta,
              dipFai=dipFai,
              bedFile=concatHapBeds.combinedBed,
              ilmMerylDBTarGz=ilmMerylDBTarGz,
              sampleID=sampleID
      }
      output {
          File insideQV = runMerqury.insideQV
          File insideMerquryTarball = runMerqury.insideMerquryTarball
          File outsideQV = runMerqury.outsideQV
          File outsideMerquryTarball = runMerqury.outsideMerquryTarball
          File Hap1PafToGRCh38=projectGIABConfToAsm.hap1ToRefPaf
          File Hap2PafToGRCh38=projectGIABConfToAsm.hap2ToRefPaf
          File GIABConfProjectionHap1=projectGIABConfToAsm.projectionBedFileHap1
          File GIABConfProjectionHap2=projectGIABConfToAsm.projectionBedFileHap2
          File Hap1PolToRawPaf=projectMosdepthToPolishedAsm.hap1ToRefPaf
          File Hap2PolToRawPaf=projectMosdepthToPolishedAsm.hap2ToRefPaf
          File MosdepthProjectionHap1=projectMosdepthToPolishedAsm.projectionBedFileHap1
          File MosdepthProjectionHap2=projectMosdepthToPolishedAsm.projectionBedFileHap2
      }

}

task bedtoolsSubtract {
  input {
      File bedFileA
      File bedFileB
      String label
      String sampleID

      Int memSizeGB = 12
      Int threadCount = 4
      Int diskSizeGB = 32
      String dockerImage = "pegi3s/bedtools"
  }

  command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        bedtools subtract -a ~{bedFileA} -b ~{bedFileB} > ~{sampleID}.~{label}.bed
  >>>

  output {
      File subtractedBed = glob("*.bed")[0]
  }

  runtime{
      memory: memSizeGB + " GB"
      cpu: threadCount
      disks: "local-disk " + diskSizeGB + " SSD"
      docker: dockerImage
  }
}


task concatBeds {
  input {
      File bedFile1
      File bedFile2
      String Prefix
      String sampleID

      Int memSizeGB = 12
      Int threadCount = 4
      Int diskSizeGB = 32
      String dockerImage = "pegi3s/bedtools"
  }

  command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cat ~{bedFile1} ~{bedFile2} > ~{sampleID}.~{Prefix}.bed
  >>>

  output {
      File combinedBed = glob("*.bed")[0]
  }

  runtime{
      memory: memSizeGB + " GB"
      cpu: threadCount
      disks: "local-disk " + diskSizeGB + " SSD"
      docker: dockerImage
  }
}
