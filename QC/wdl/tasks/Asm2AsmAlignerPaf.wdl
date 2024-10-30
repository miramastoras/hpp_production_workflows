version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t

workflow asm2asmAlignerPaf{
  input {
      Int kmerSize
      String aligner
      String dockerImage
      File refAssembly
      File queryAssembly
      String preset

  }
  call long_read_aligner_t.alignmentPaf as alignmentPaf{
      input:
          preset=preset,
          aligner=aligner,
          readFastq_or_queryAssembly=queryAssembly,
          refAssembly=refAssembly,
          diskSize=512,
          threadCount=32,
          kmerSize=kmerSize,
          dockerImage=dockerImage
      }
  output {
      File asm2asmPaf=alignmentPaf.pafFile
  }
}
