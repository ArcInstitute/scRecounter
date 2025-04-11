// Subworkflows
include { DB_ACC_WF } from './workflows/db_acc.nf'
include { STAR_PARAMS_WF } from './workflows/star_params.nf'
include { STAR_FULL_WF } from './workflows/star_full.nf'
include { SRA_STAT } from './lib/utils.nf'
// util functions
include { readAccessions; addStats; } from './lib/utils.nf'

// Main workflow
workflow { 
    if (params.accessions == "" || params.accessions == true) {
        // Obtain accessions from the scRecoutner SQL database
        println "No accessions provided. Accessions will be obtained from the scRecounter database"
        ch_accessions = DB_ACC_WF()
    } else {
        // Use the provided accessions
        println "Using provided accessions."
        ch_accessions = Channel.fromPath(params.accessions, checkIfExists: true)
    }

    // read accessions file
    ch_accessions = readAccessions(ch_accessions)

    // run sra-stat on accessions
    ch_sra_stat = SRA_STAT(ch_accessions)
    ch_accessions = addStats(ch_accessions, ch_sra_stat)

    // filter out any accessions with max SRA file size greater than the user-specified size
    ch_accessions = ch_accessions.filter { it[3] <= params.max_sra_size }
    
    // determine best STAR parameters on a subset of reads
    ch_star_params = STAR_PARAMS_WF(ch_accessions, ch_sra_stat)

    // run STAR on all reads with selected parameters
    if (! params.define){
        STAR_FULL_WF(ch_accessions, ch_star_params)
    }
}

// On complete
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
