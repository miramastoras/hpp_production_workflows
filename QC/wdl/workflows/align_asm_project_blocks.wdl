version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/project_blocks.wdl" as project_blocks_t

workflow align_asm_project_blocks {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "align diploid assemblies and project blocks"
      }

    input {
      File asmHap1Fasta
      File asmHap2Fasta

      File refHap1Fasta
      File refHap2Fasta

      File bedFile

      String projectionDirection="ref2asm"

      String sampleID
      String pafAligner="minimap2"
      }

    # Align hap1 to hap1, in paf format
    call long_read_aligner_t.alignmentPaf as alignHap1ToRef{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=asmHap1Fasta,
            refAssembly=refHap1Fasta,
            suffix="${sampleID}.asmToRefHap1",
            diskSize=300,
            threadCount=32,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }
    call long_read_aligner_t.alignmentPaf as alignHap2ToRef{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=asmHap2Fasta,
            refAssembly=refHap2Fasta,
            suffix="${sampleID}.asmToRefHap2",
            diskSize=300,
            threadCount=32,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }
    # project blocks hap1
    call project_blocks_t.project_blocks as projectBedHap1 {
        input:
            pafFile=alignHap1ToRef.pafFile,
            bedFile=bedFile,
            mode=projectionDirection
    }
    call project_blocks_t.project_blocks as projectBedHap2 {
        input:
            pafFile=alignHap2ToRef.pafFile,
            bedFile=bedFile,
            mode=projectionDirection
    }
    output {
        File projectionBedFileHap1 = projectBedHap1.projectionBedFile
        File projectableBedFileHap1 = projectBedHap1.projectableBedFile
        File projectionBedFileHap2 = projectBedHap2.projectionBedFile
        File projectableBedFileHap2 = projectBedHap2.projectableBedFile
        File hap1ToRefPaf=alignHap1ToRef.pafFile
        File hap2ToRefPaf=alignHap2ToRef.pafFile

    }
}
