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
              row["organism"], row["matrix_type"], row["srx"], 
              file(row["matrix_path"]), 
              file(row["features_path"]), 
              file(row["barcodes_path"])
            )
        }

    // filter to just one organism
    //mtx_files = mtx_files.filter { organism,matrix_type,srx,m_path,f_path,b_path -> organism == "Mus_musculus" }
    //mtx_files = mtx_files.filter { organism,matrix_type,srx,m_path,f_path,b_path -> organism == "Homo_sapiens" }

    // aggregate mtx files as h5ad
    MTX_TO_H5AD( mtx_files )

    // register
    H5AD_REGISTER( MTX_TO_H5AD.out.h5ad.groupTuple(by:[0,1]) )

    // join MTX_TO_H5AD.out.h5ad and H5AD_REGISTER.out.pkl on `organism`
    h5ad_files = MTX_TO_H5AD.out.h5ad.combine( H5AD_REGISTER.out.pkl, by: 0 )

    // add the h5ad files to the database
    H5AD_TO_DB( h5ad_files )
}

process H5AD_TO_DB {
    publishDir file(params.log_dir), mode: "copy", overwrite: true
    label "process_medium"
    maxForks 1

    input:
    tuple val(organism), val(mtx_type), val(srx), path(h5ad), path(pkl)

    output:
    path "h5ad-to-db_${srx}.log", emit: log

    script:
    """
    set -o pipefail
    h5ad-to-db.py \\
      --db-uri ${params.db_uri} \\
      --organism ${organism} \\
      --matrix-type ${mtx_type} \\
      --registration-plan ${pkl} \\
      $h5ad 2>&1 | tee h5ad-to-db_${srx}.log
    """
}

process H5AD_REGISTER {
  publishDir file(params.log_dir), mode: "copy", overwrite: true, pattern: "*.log"
  label "process_medium"
  maxForks 1

  input:
  tuple val(organism), val(mtx_type), val(srx), path(h5ad)

  output:
  tuple val(organism), path("registration-plan.pkl"), emit: pkl
  path "h5ad-register.log", emit: log

  script:
  """
  set -o pipefail
  h5ad-register.py \\
    --db-uri ${params.db_uri} \\
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
    tuple val(organism), val(mtx_type), val(srx), path("matrix.mtx.gz"), path("features.tsv.gz"), path("barcodes.tsv.gz")

    output:
    tuple val(organism), val(mtx_type), val(srx), path("${srx}.h5ad"), emit: h5ad
    path "mtx-to-h5ad_${srx}.log", emit: log

    script:
    """
    set -o pipefail
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    mtx-to-h5ad.py \\
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
    set -o pipefail
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