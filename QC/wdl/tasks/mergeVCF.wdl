version 1.0

import "../tasks/hapDotPy.wdl" as runHappy
import "../workflows/snv_indel_assembly.wdl" as snvIndelAsm

workflow snv_indel_assembly {

    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "merge two VCFs using hap.py and Kishwar's python script"
    }

    input {
        File assembly
        File assemblyIndex
        String excludeExpr = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'"
        String excludeTypes = "indels"
        String sample
    }

    ## Compare filtered callsets (DeepVariant & PMDV) to see where they agree
    call runHappy.hapDotPy as happy_t {
        input:
            truthVCF      = filt_1.vcfOut,
            queryVCF      = filt_2.vcfOut,
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample
    }

    ## Output the union of the two filtered callsets (DeepVariant & PMDV)
    call mergeVCF {
        input:
            VCF1     = filt_1.vcfOut,
            VCF2     = filt_2.vcfOut,
            happyVCF = happy_t.vcfOut,
            sample   = sample
    }

    ## create text file w/ stats of final output from the last step
    call createVCFStats {
        input:
            inputVCF = mergeVCF.vcfOut
    }
