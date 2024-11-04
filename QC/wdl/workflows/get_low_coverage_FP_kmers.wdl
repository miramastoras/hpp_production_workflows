version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/annotate_edit_with_fp_kmers.wdl" as annotate_edit_with_fp_kmers_t
import "./align_asm_project_blocks.wdl" as align_asm_project_blocks_t

workflow get_low_coverage_FP_kmers {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "get number of induced, fixed, unchanged error kmers intersecting low coverage regions on both haplotypes"
      }

    input {
      File Hap1RawFasta
      File Hap2RawFasta

      File Hap1PolFasta
      File Hap2PolFasta

      File rawMerquryTarGz
      File polishedMerquryTarGz
      File polishingVcf

      File lowCoverageBed

      String sampleID
      }

      # align hap1 pol to hap1 raw
      call long_read_aligner_t.alignmentPaf as alignPolToRawHap1{
          input:
              aligner="minimap2",
              preset="asm5",
              options="-L --eqx --cs -c",
              readFastq_or_queryAssembly=Hap1PolFasta,
              refAssembly=Hap1RawFasta,
              suffix="${sampleID}.PolToRawHap1",
              diskSize=300,
              threadCount=32,
              kmerSize=19,
              dockerImage="mobinasri/long_read_aligner:v0.3.3"
      }
      # align hap2 pol to hap2 raw
      call long_read_aligner_t.alignmentPaf as alignPolToRawHap2{
          input:
              aligner="minimap2",
              preset="asm5",
              options="-L --eqx --cs -c",
              readFastq_or_queryAssembly=Hap2PolFasta,
              refAssembly=Hap2RawFasta,
              suffix="${sampleID}.PolToRawHap2",
              diskSize=300,
              threadCount=32,
              kmerSize=19,
              dockerImage="mobinasri/long_read_aligner:v0.3.3"
      }

      # annotate edits with fp kmers
      call annotate_edit_with_fp_kmers_t.annotateVCFwithFPKmers as annotate_edits {
          input:
              hap1PolishedToRawPaf=alignPolToRawHap1.pafFile,
              hap2PolishedToRawPaf=alignPolToRawHap2.pafFile,
              polishedMerquryTarGz=polishedMerquryTarGz,
              rawMerquryTarGz=rawMerquryTarGz,
              deeppolisherVcfGz=polishingVcf,
              sample=sampleID
      }
      # align asm project blocks (low cov) hap1 to hap2 and vice versa
      call align_asm_project_blocks_t.align_asm_project_blocks as projectlowCovBlocks {
          input:
              asmHap1Fasta=Hap2RawFasta,
              asmHap2Fasta=Hap1RawFasta,
              refHap1Fasta=Hap1RawFasta,
              refHap2Fasta=Hap2RawFasta,
              bedFile=lowCoverageBed,
              projectionDirection="ref2asm",
              sampleID=sampleID

      }
      output {
          File inducedFPKmerBed = annotate_edits.inducedFPKmerBed
          File fixedFPKmerBed = annotate_edits.fixedFPKmerBed
          File unchangedFPKmerBed = annotate_edits.unchangedFPKmerBed
          File Hap1ToHap2ProjectedBed = projectlowCovBlocks.projectionBedFileHap1
          File Hap2ToHap1ProjectedBed = projectlowCovBlocks.projectionBedFileHap2
          File hap1RawToHap2RawPaf = projectlowCovBlocks.hap1ToRefPaf
          File hap1PolishedToRawPaf = alignPolToRawHap1.pafFile
          File hap2PolishedToRawPaf = alignPolToRawHap2.pafFile
      }
  }
