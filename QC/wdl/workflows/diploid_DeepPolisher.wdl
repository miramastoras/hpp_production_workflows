version 1.0

import "../tasks/sepReadsByHaplotype.wdl" as sepReadsByHap_t
import "../tasks/DeepPolisher.wdl" as DeepPolisher_t

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
        File maternalRawFastaIndex

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

    call DeepPolisher_t.runDeepPolisher as runDeepPolisherHap1 {
        input:
          Bam=sepPhasedReadsByHap.hap1Bam,
          Bai=sepPhasedReadsByHap.hap1Bai,
          Fasta=paternalRawFasta,
          ModelFilesTarGZ=ModelFilesTarGZ,
          dockerImage=DeepPolisherDocker,
          sampleName=sampleName
    }

    call DeepPolisher_t.runDeepPolisher as runDeepPolisherHap2 {
        input:
          Bam=sepPhasedReadsByHap.hap2Bam,
          Bai=sepPhasedReadsByHap.hap2Bai,
          Fasta=maternalRawFasta,
          ModelFilesTarGZ=ModelFilesTarGZ,
          dockerImage=DeepPolisherDocker,
          sampleName=sampleName
    }

    output {
        File Hap1_DeepPolisherVcf=runDeepPolisherHap1.PolisherVcf
        File Hap1_DeepPolisherVcfTbi=runDeepPolisherHap1.PolisherVcfTbi

        File Hap2_DeepPolisherVcf=runDeepPolisherHap2.PolisherVcf
        File Hap2_DeepPolisherVcfTbi=runDeepPolisherHap2.PolisherVcfTbi
    }
}
