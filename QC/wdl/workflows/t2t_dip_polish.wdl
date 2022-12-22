version 1.0

import "../tasks/deepvariant.wdl" as runDeepVariant
import "../tasks/pepperMarginDeepVariant.wdl" as runPepperMarginDeepVariant
import "../tasks/hapDotPy.wdl" as runHappy

workflow snv_indel_assembly {

    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Application of small variant polishing workflow used on T2T-CHM13 to diploid assemblies [T2T Polishing Case Study](https://github.com/arangrhie/T2T-Polish/blob/master/doc/T2T_polishing_case_study.md)."
    }

    input {
        File dipAssembly
        File dipAssemblyIndex
        String excludeExpr = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'"
        String excludeTypes = "indels"
        String sample
    }
