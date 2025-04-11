include { saveAsLog } from '../lib/utils.nf'
include { readAccessions } from '../lib/download.nf'

workflow DB_ACC_WF {
    main:
    // throw error if --use-database is false
    if (!params.use_database) {
        throw new Exception("--use-database is required to obtain accessions from the scRecounter SQL database")
    }

    // obtain accessions from the database
    ch_accessions = GET_DB_ACCESSIONS()
    ch_accessions.csv.ifEmpty { println 'No accessions found in the scRecounter database' }

    emit:
    ch_accessions.csv
}

// Save accessions csv
def saveAsFinalAcc(filename) {
    if (filename.endsWith(".csv")){
        filename = filename.tokenize("/").last()
        return "${filename}"
    }
    return null
}

process GET_DB_ACCESSIONS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, pattern: "*.csv"
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename) }
    label "download_env"
    disk 10.GB

    output:
    path "accessions.csv",      emit: "csv"
    path "${task.process}.log", emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    get-db-accessions.py \\
      --organisms "${params.organisms}" \\
      --max-srx ${params.max_samples} \\
      2>&1 | tee ${task.process}.log
    """
}

