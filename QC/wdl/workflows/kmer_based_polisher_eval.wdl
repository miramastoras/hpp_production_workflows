version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/merqury.wdl" as merqury_t
import "../tasks/yak.wdl" as yak_t
import "../tasks/yak_non_trio.wdl" as yak_non_trio_t
import "../tasks/project_blocks.wdl" as project_blocks_t
import "../tasks/subFastaByBed.wdl" as subset_fasta_t

workflow kmerPolishingEval {
    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Evaluate diploid assembly with kmer-based metrics (QV,switch error,hamming error). Assumes meryl dbs and yak files already made"
      }

    input {
      File hap1Fasta
      File hap2Fasta

      File grch38Fasta
      File grch38InsideConfRegions
      File grch38OutsideConfRegions

      String sampleID
      String pafAligner="minimap2"

      File ilmMerylDBTarGz
      File sampleYak

      # if enableYakTrioEval = false, pass sampleYak as paternalYak and maternalYak
      Boolean enableYakTrioEval = false
      File paternalYak
      File maternalYak

      }

    # Align hap1 and hap2 to grch38, in paf format
    call long_read_aligner_t.alignmentPaf as alignHap1ToRef{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=hap1Fasta,
            refAssembly=grch38Fasta,
            suffix="asmToRef",
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
            readFastq_or_queryAssembly=hap2Fasta,
            refAssembly=grch38Fasta,
            suffix="asmToRef",
            diskSize=300,
            threadCount=32,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }
    # project inside and outside confidence regions to Hap1 and Hap2 coordinates
    call project_blocks_t.project_blocks as projectInsideConfHap1 {
        input:
            pafFile=alignHap1ToRef.pafFile,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t.project_blocks as projectInsideConfHap2 {
        input:
            pafFile=alignHap2ToRef.pafFile,
            bedFile=grch38InsideConfRegions
    }
    call project_blocks_t.project_blocks as projectOutsideConfHap1 {
        input:
            pafFile=alignHap1ToRef.pafFile,
            bedFile=grch38OutsideConfRegions
    }
    call project_blocks_t.project_blocks as projectOutsideConfHap2 {
        input:
            pafFile=alignHap2ToRef.pafFile,
            bedFile=grch38OutsideConfRegions
    }

    call subset_fasta_t.SubFastaByBed as subHap1InsideConf {
        input:
            Fasta=hap1Fasta,
            Bed=projectInsideConfHap1.projectionBedFile,
            outputLabel="hap1.insideConf",
            sampleID=sampleID
    }
    call subset_fasta_t.SubFastaByBed as subHap2InsideConf {
        input:
            Fasta=hap2Fasta,
            Bed=projectInsideConfHap2.projectionBedFile,
            outputLabel="hap2.insideConf",
            sampleID=sampleID
    }
    call subset_fasta_t.SubFastaByBed as subHap1OutsideConf {
        input:
            Fasta=hap1Fasta,
            Bed=projectOutsideConfHap1.projectionBedFile,
            outputLabel="hap1.outsideConf",
            sampleID=sampleID
    }
    call subset_fasta_t.SubFastaByBed as subHap2OutsideConf {
        input:
            Fasta=hap2Fasta,
            Bed=projectOutsideConfHap2.projectionBedFile,
            outputLabel="hap2.outsideConf",
            sampleID=sampleID
    }

    # Run merqury QV whole genome
    call merqury_t.merqury as merquryWholeGenome {
        input:
            assemblyFasta=hap1Fasta,
            altHapFasta=hap2Fasta,
            kmerTarball=ilmMerylDBTarGz
    }
    call merqury_t.merqury as merquryInsideConf {
        input:
            assemblyFasta=subHap1InsideConf.subFasta,
            altHapFasta=subHap2InsideConf.subFasta,
            kmerTarball=ilmMerylDBTarGz
    }
    call merqury_t.merqury as merquryOutsideConf {
        input:
            assemblyFasta=subHap1OutsideConf.subFasta,
            altHapFasta=subHap2OutsideConf.subFasta,
            kmerTarball=ilmMerylDBTarGz
    }

    # run yak trio eval on whole genome
    if (enableYakTrioEval) {
        call yak_t.yakAssemblyStats as yakQCWholeGenome {
            input:
                matYak=maternalYak,
                patYak=paternalYak,
                sampleYak=sampleYak,
                assemblyFastaPat=hap1Fasta,
                assemblyFastaMat=hap2Fasta,
                minSequenceLength="0",
                dockerImage="miramastoras/hpp_yak:latest"
        }
    }

    if (enableYakTrioEval == false) {
        call yak_non_trio_t.yakNonTrioAssemblyStats as yakQCWholeGenomeNonTrio {
            input:
                assemblyFastaHap2=hap2Fasta,
                assemblyFastaHap1=hap1Fasta,
                sampleYak=sampleYak,
                minSequenceLength="0",
                dockerImage="miramastoras/hpp_yak:latest"
        }
    }

    # Run Yak QV inside and outside conf

    call yak_non_trio_t.yakNonTrioAssemblyStats as yakQCInsideConf {
        input:
            sampleYak=sampleYak,
            assemblyFastaHap1=subHap1InsideConf.subFasta,
            assemblyFastaHap2=subHap2InsideConf.subFasta,
            minSequenceLength="0",
            dockerImage="miramastoras/hpp_yak:latest"
    }
    call yak_non_trio_t.yakNonTrioAssemblyStats as yakQCOutsideConf {
        input:
            sampleYak=sampleYak,
            assemblyFastaHap1=subHap1OutsideConf.subFasta,
            assemblyFastaHap2=subHap2OutsideConf.subFasta,
            minSequenceLength="0",
            dockerImage="miramastoras/hpp_yak:latest"
    }

    File yakWGTarBall = select_first([yakQCWholeGenome.outputTarball, yakQCWholeGenomeNonTrio.outputTarball])

    output {
        File QV_whole_genome = merquryWholeGenome.QV
        File QV_inside_conf = merquryInsideConf.QV
        File QV_outside_conf = merquryOutsideConf.QV
        File merquryAsmFPkmers = merquryWholeGenome.asmFPkmers
        File merquryAltHapFPkmers = merquryWholeGenome.altHapFPkmers
        File merquryWGTarBall=merquryWholeGenome.outputTarball
        File merquryInsideConfTarBall=merquryInsideConf.outputTarball
        File merquryOutsideConfTarBall=merquryOutsideConf.outputTarball
        File yakTarBallWG=yakWGTarBall
        File yakTarBallInsideConf=yakQCInsideConf.outputTarball
        File yakTarBallOutsideConf=yakQCOutsideConf.outputTarball
        File hap1InsideConfFasta=subHap1InsideConf.subFasta
        File hap2InsideConfFasta=subHap2InsideConf.subFasta
        File hap1OutsideConfFasta=subHap1OutsideConf.subFasta
        File hap2OutsideConfFasta=subHap2OutsideConf.subFasta
    }
}
