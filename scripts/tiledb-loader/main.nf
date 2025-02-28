workflow { 
    // find target MTX files to add to the database
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out.csv
        .splitCsv( header: true )
        .map { row -> 
            tuple( 
              row["batch"], row["matrix_type"], row["organism"],
              row["srx"], file(row["matrix_path"]), file(row["features_path"]), file(row["barcodes_path"])
            )
        }.groupTuple(by:[0,1,2])
        
    //mtx_files.view()

    // group Velocyto MTX files by SRX
    /*
    if( params.feature_type == "Velocyto"){
      mtx_files = mtx_files.groupTuple(by: [0,1]).map{ group -> 
        tuple(group[0], group[1], group[2], group[3][0], group[4][0])
      }
    } else {
      mtx_files = mtx_files.groupTuple()
    }
    */
    
    // aggregate mtx files as h5ad
    MTX_TO_H5AD( mtx_files )

    // add the h5ad files to the database
    H5AD_TO_DB( MTX_TO_H5AD.out.h5ad )
}

process H5AD_TO_DB {
    publishDir file(params.log_dir), mode: "copy", overwrite: true
    label "process_medium"
    maxForks 1

    input:
    path h5ad

    output:
    path "h5ad-to-db.log", emit: log

    script:
    """
    h5ad-to-db.py \\
      --feature-type ${params.feature_type} \\
      --db-uri ${params.db_uri} \\
      $h5ad 2>&1 | tee h5ad-to-db.log
    """
}

process MTX_TO_H5AD {
    publishDir file(params.log_dir) , mode: "copy", overwrite: true, pattern: "*.log"
    label "process_high"
    maxForks 4

    input:
    tuple val(batch), val(mtx_type), val(organism), val(srx), path("*_matrix.mtx.gz"), path("*_features.tsv.gz"), path("*_barcodes.tsv.gz")

    output:
    path "data.h5ad",                      emit: h5ad
    path "mtx-to-h5ad_batch-${batch}.log", emit: log

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    mtx-to-h5ad.py \\
      --threads ${task.cpus} \\
      --feature-type "${params.feature_type}" \\
      --missing-metadata "${params.missing_metadata}" \\
      --srx "$srx" \\
      --mtx-path *_matrix.mtx \\
      2>&1 | tee mtx-to-h5ad_batch-${batch}.log
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