version 1.0

import "../tasks/sepReadsByHaplotype.wdl" as sepReadsByHap_t
import "../tasks/DeepPolisher.wdl" as deepPolisher_t

workflow diploid_DeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Workflow for running DeepPolisher on diploid assemblies"
    }

    input {
        File hifiPhasedDipBam
        File hifiPhasedDipBai

        File paternalRawFasta
        File paternalRawFastaIndex
        File maternalRawFasta
        File maternalRawFastaIdx

        File ModelFilesTarGZ
        File DeepPolisherDocker

        String sampleName

    }

    call sepReadsByHap_t.Separate as sepPhasedReadsByHap {
        input:
          dipBam=hifiPhasedDipBam,
          hap1Fai=paternalRawFastaIndex,
          hap2Fai=maternalRawFastaIndex
    }

    call deepPolisher_t.runDeepPolisher as runDeepPolisherHap1 {
        input:
          Bam=sepPhasedReadsByHap.hap1Bam,
          Bai=sepPhasedReadsByHap.hap1Bai,
          Fasta=paternalRawFasta,
          ModelFilesTarGZ=ModelFilesTarGZ,
          dockerImage=DeepPolisherDocker,
          sampleName=sampleName
    }

    call deepPolisher_t.runDeepPolisher as runDeepPolisherHap2 {
        input:
          Bam=sepPhasedReadsByHap.hap2Bam,
          Bai=sepPhasedReadsByHap.hap2Bai,
          Fasta=maternalRawFasta,
          ModelFilesTarGZ=ModelFilesTarGZ,
          dockerImage=DeepPolisherDocker,
          sampleName=sampleName
    }

    output {
        File Hap1_deepPolisherVcf=runDeepPolisherHap1.PolisherVcf
        File Hap1_deepPolisherVcfTbi=runDeepPolisherHap1.PolisherVcfTbi

        File Hap2_deepPolisherVcf=runDeepPolisherHap2.PolisherVcf
        File Hap2_deepPolisherVcfTbi=runDeepPolisherHap2.PolisherVcfTbi
    }
