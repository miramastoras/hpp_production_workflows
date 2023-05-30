version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t

workflow asm2asmAlignerPaf{

  call long_read_aligner_t.alignmentPaf as alignmentPaf{
      input:
          aligner="winnowmap",
          preset="asm5",
          options="-L --eqx --cs -c",
          suffix="mat2pat",
          diskSize=512,
          threadCount=64,
          kmerSize=19,
          dockerImage="mobinasri/long_read_aligner:v0.2"
      }
  output {
      asm2asmPaf=alignmentPaf.pafFile
  }
}
