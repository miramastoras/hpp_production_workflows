version 1.0

import ../tasks/asm2asm_aligner.wdl as asm2asm_aligner_t
workflow phasingHomozygous{

    input {
        File paternalFastaGz
        File maternalFastaGz
    }

    ## Align maternal to paternal assembly
    call asm2asm_aligner_t as asm2asm_aligner{
        input:
            preset="asm5",
            queryAssemblyFastaGz=maternalFastaGz,
            refAssemblyFastaGz=paternalFastaGz
    }

    ## Get Homozygous regions 

}
