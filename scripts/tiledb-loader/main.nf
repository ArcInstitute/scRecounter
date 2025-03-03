workflow { 
    // find target MTX files to add to the database
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out.csv.splitCsv( header: true )
        .map { row -> 
            // remove spaces from organism
            row["organism"] = row["organism"].replaceAll(" ", "_")
            // return tuple
            tuple( 
              row["matrix_type"], row["organism"], row["srx"], 
              file(row["matrix_path"]), file(row["features_path"]), file(row["barcodes_path"])
            )
        }

    //mtx_files.view()

    // aggregate mtx files as h5ad
    MTX_TO_H5AD( mtx_files )

    // group by h5ad files by matrix type and organism
    h5ad_files = MTX_TO_H5AD.out.h5ad.groupTuple(by:[0,1])

    // register
    H5AD_REGISTER( h5ad_files )

    // add the h5ad files to the database
    H5AD_TO_DB( MTX_TO_H5AD.out.h5ad, H5AD_REGISTER.out.pkl )
}

process H5AD_TO_DB {
    publishDir file(params.log_dir), mode: "copy", overwrite: true
    label "process_medium"
    maxForks 1

    input:
    tuple val(mtx_type), val(organism), val(srx), path(h5ad)
    each pkl

    output:
    path "h5ad-to-db_${srx}.log", emit: log

    script:
    """
    h5ad-to-db.py \\
      --feature-type ${params.feature_type} \\
      --organism ${organism} \\
      --db-uri ${params.db_uri} \\
      --registration-plan ${pkl} \\
      $h5ad 2>&1 | tee h5ad-to-db_${srx}.log
    """
}

process H5AD_REGISTER {
    publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
    label "process_medium"
    maxForks 1

    input:
    tuple val(organism), val(srx), path(h5ad)

    output:
    tuple "registration-plan.pkl", emit: pkl
    path "h5ad-register.log",     emit: log

    script:
    """
    h5ad-register.py \\
      --db-uri ${params.db_uri} \\
      --feature-type ${params.feature_type} \\
      --matrix-type ${mtx_type} \\
      --organism ${organism} \\
      $h5ad 2>&1 | tee h5ad-register.log
    """
}

process MTX_TO_H5AD {
    publishDir file(params.log_dir) , mode: "copy", overwrite: true, pattern: "*.log"
    label "process_high"
    maxForks 200

    input:
    tuple val(mtx_type), val(organism), val(srx), path("matrix.mtx.gz"), path("features.tsv.gz"), path("barcodes.tsv.gz")

    output:
    tuple val(mtx_type), val(organism), val(srx), path("${srx}.h5ad"), emit: h5ad
    path "mtx-to-h5ad_${srx}.log", emit: log

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    mtx-to-h5ad.py \\
      --feature-type "${params.feature_type}" \\
      --missing-metadata "${params.missing_metadata}" \\
      --srx $srx \\
      --mtx-path matrix.mtx.gz \\
      2>&1 | tee mtx-to-h5ad_${srx}.log
    """
}

process FIND_MTX {
    publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    output:
    path "mtx_files.csv", emit: csv
    path "find-mtx.log",  emit: log

    script:
    def organisms = params.organisms != "" ? "--organisms \"${params.organisms}\"" : ""
    def redo_processed = params.redo_processed.toString() == "true" ? "--redo-processed" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    find-mtx.py ${organisms} ${redo_processed} \\
      --feature-type ${params.feature_type} \\
      --max-datasets ${params.max_datasets} \\
      --batch-size ${params.mtx_batch_size} \\
      --db-uri ${params.db_uri} \\
      ${params.input_dir} \\
      2>&1 | tee find-mtx.log
    """
}