version 1.0

import "../tasks/meryl_hybrid.wdl" as meryl_t
import "../tasks/merqury.wdl" as merqury_t

workflow runMerquryHybrid {
    call meryl_t as runMerylHybrid
    call merqury_t as runMerqury

    output {
        File QV = runMerqury.QV
        File outputTarball = merqury.outputTarball
        File altHapFPkmers = merqury.altHapFPkmers
        File asmFPkmers = merqury.asmFPkmers
    }
}
