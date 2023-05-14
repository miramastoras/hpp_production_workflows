version 1.0

import "pepperMarginDeepVariant.wdl" as PMDV


workflow filterVCF {

    call PMDV.bcftoolsFilter as filter_t

    output {
        File filtVcfOut=filter_t.vcfOut
        File filtVcfOutIdx=filter_t.vcfOutIdx
    }
}
