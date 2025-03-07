workflow { 
    // find target MTX files to add to the database
    FIND_MTX()

    // list target MTX files
    mtx_files = FIND_MTX.out.csv
      .splitCsv( header: true )
      .map{ row -> 
        tuple(row.srx, file(row.matrix_path), file(row.features_path), file(row.barcodes_path))
      }

    // group Velocyto MTX files by SRX
    if( params.feature_type == "Velocyto"){
      mtx_files = mtx_files.groupTuple().map{ group -> 
        tuple(group[0], group[1], group[2][0], group[3][0])
      }
    }

    // convert to h5ad and publish
    MTX_TO_H5AD( mtx_files, Channel.fromPath(params.tissue_categories) )

    // write parquet after all MTX_TO_H5AD jobs complete
    if( params.update_db ){
      DB_TO_PARQUET( MTX_TO_H5AD.out.h5ad.collect() )
    }
    
    // aggregate obs metadata
    AGG_OBS_METADATA( MTX_TO_H5AD.out.csv.collate(100) )
}

process AGG_OBS_METADATA {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, pattern: "metadata_TMP/${params.feature_type}/*.csv.gz"
    publishDir file(params.log_dir) / params.feature_type, mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    input:
    path csv_files

    output:
    path "metadata_TMP/${params.feature_type}/*.csv.gz", emit: obs_meta
    path "agg-obs-metadata.log",                         emit: log

    script:
    """
    agg-obs-metadata.py \\
      --feature-type ${params.feature_type} \\
      ${csv_files} 2>&1 | tee agg-obs-metadata.log
    """
}

process DB_TO_PARQUET {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, pattern: "metadata/${params.feature_type}/*/sample_metadata.parquet.gz"
    publishDir file(params.log_dir) / params.feature_type, mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    input:
    path csv_files

    output:
    path "metadata/${params.feature_type}/*/sample_metadata.parquet.gz", emit: samp_meta
    path "db-to-parquet.log",                                            emit: log

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    db-to-parquet.py \\
      --feature-type ${params.feature_type} \\
      2>&1 | tee db-to-parquet.log
    """
}

process MTX_TO_H5AD {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, pattern: "h5ad/${params.feature_type}/*/*.h5ad.gz"
    publishDir file(params.log_dir) / params.feature_type, mode: "copy", overwrite: true, pattern: "*.log"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    label "process_low"
    maxForks 200

    input:
    tuple val(srx), path(mtx_path), path(features_path), path(barcodes_path)
    each path(tissue_categories)

    output:
    path "h5ad/${params.feature_type}/*/${srx}.h5ad.gz",  emit: h5ad
    path "metadata/${srx}.csv.gz", emit: csv
    path "mtx-to-h5ad_${srx}.log", emit: log

    script:
    def update_db = params.update_db ? "--update-database" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    mtx-to-h5ad.py ${update_db} \\
      --feature-type ${params.feature_type} \\
      --missing-metadata "${params.missing_metadata}" \\
      --tissue-categories "${tissue_categories}" \\
      --srx ${srx} \\
      --matrix ${mtx_path} \\
      --publish-path "${params.output_dir}" \\
      2>&1 | tee mtx-to-h5ad_${srx}.log
    """
}

process FIND_MTX {
    publishDir file(params.log_dir) / params.feature_type, mode: "copy", overwrite: true, pattern: "*.log"
    label "process_low"

    output:
    path "mtx_files.csv", emit: csv
    path "find-mtx.log",  emit: log

    script:
    def multi_mapper = params.multi_mapper.toString().tokenize(",").collect{ it.trim() }.join(" ")
    def organisms = params.organisms != "" ? "--organisms \"${params.organisms}\"" : ""
    def redo_processed = params.redo_processed.toString() == "true" ? "--redo-processed" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    find-mtx.py ${organisms} ${redo_processed} \\
      --feature-type ${params.feature_type} \\
      --multi-mapper ${multi_mapper} \\
      --max-datasets ${params.max_datasets} \\
      ${params.input_dir} \\
      2>&1 | tee find-mtx.log
    """
}