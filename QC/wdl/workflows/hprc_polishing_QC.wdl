version 1.0

import "./kmer_based_polisher_eval.wdl" as kmer_based_polisher_eval_wf
import "../tasks/project_blocks.wdl" as project_blocks_t
import "../"

workflow hprc_polishing_QC {

    meta {
      author: "Mira Mastoras"
      email: "mmastora@ucsc.edu"
      description: "Runs kmer based QC for hprc assemblies, before and after polishing"
      }

    input {

    }
    call meryl

    call kmer based qc {
        inputs:
    }
    call project fp kmers

    call intersect with polishing edits

    call collate all results to csv
}
