version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t

workflow projectBlocks {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Project bed coordinates from one genome assembly to another"
      }

    input {
      File assemblyFasta
      File refFasta
      File bedFile
      String pafAligner="minimap2"
      String mode

    }
    call long_read_aligner_t.alignmentPaf as alignAsm2Ref{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=assemblyFasta,
            refAssembly=refFasta,
            suffix="asmToRef",
            diskSize=512,
            threadCount=64,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    call project_blocks as projectBlocks {
        input:
            pafFile=alignAsm2Ref.pafFile,
            bedFile=bedFile,
            mode=mode
    }
    output {
        File projectionBedFile = projectBlocks.projectionBedFile
        File projectableBedFile = projectBlocks.projectableBedFile
    }
}

task project_blocks{
    input {
        File pafFile
        File bedFile
        String mode='ref2asm'
        Int memSizeGB = 12
        Int threadCount = 16
        Int diskSizeGB = 256
        String dockerImage = "mobinasri/flagger:latest"
      }

      command <<<
            set -o pipefail
            set -e
            set -u
            set -o xtrace

            pafBasename=$(basename ~{pafFile})
            bedBasename=$(basename ~{bedFile})

            python3 /home/programs/src/project_blocks_multi_thread.py \
            --threads ~{threadCount} \
            --mode ~{mode} \
            --paf {pafFile} \
            --blocks ~{bedFile} \
            --outputProjectable ${fileBasename}_${bedBasename}.projectable.bed \
            --outputProjection ${fileBasename}_${bedBasename}.projection.bed
    	>>>

      output {
          File projectionBedFile = glob("*.projection.bed")[0]
          File projectableBedFile = glob("*.projectable.bed")[0]
      }
}
