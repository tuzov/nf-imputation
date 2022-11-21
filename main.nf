// Imputation 1.0 : Pipeline for the imputation of a target dataset against a reference 
// main.nf -- the main imputation workflow that unites imputation itself, qc and etc.
nextflow.enable.dsl=2 

// Modules
include { RUN_IMPUTE; COLLECT_IMPUTE } from "./modules/impute.nf"
//include { RUN_QC; COLLECT_QC } from "./modules/qc.nf"

// Channels


workflow {
  // impute workflow

  // qc workflow


}

