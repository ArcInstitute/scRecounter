include { joinReads; saveAsLog; } from '../lib/utils.nf'
include { addSaSize; starFullMem; } from '../lib/star_full.nf'
// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_FULL_WF{
    take:
    ch_accessions
    ch_star_params
    
    main:
    //-- Download all reads --//
    // filter out samples that lack a set of selected parameters
    ch_accessions_filt = ch_accessions.combine(
        ch_star_params.map{ it[0] }.unique(), by: 0
    )
    
    // group by sample
    ch_accessions_filt_group = ch_accessions_filt.map{
        sample, accession, metadata, star_params -> [sample, accession]
    }.groupTuple()

    // Use xsra to download all reads
    XSRA(ch_accessions_filt_group)
    ch_fastq = XSRA.out.R1.join(XSRA.out.R2)
    ch_fastq.count().view{ count -> "XSRA accession count: $count" }

    //-- For accessions lacking paired reads from fasterq-dump, fallback to fastq-dump --//
    ch_accessions_fallback = XSRA.out.R1
        .join(
            XSRA.out.R2.map{ sample,fastq -> [sample,fastq,true] },
            remainder: true
        )
        .filter{ it -> it[3] != true }
        .map{ it -> it[0] } 
        .combine(ch_accessions_filt_group, by: 0)

    // run fastq-dump on the fallback accessions
    FASTQ_DUMP(ch_accessions_fallback)
    ch_fastq_fallback = FASTQ_DUMP.out.R1.join(FASTQ_DUMP.out.R2)
    ch_fastq_fallback.count().view{ count -> "fastq-dump (fallback) accession count: $count" }
    
    // combine the fasterq-dump and fastq-dump results
    ch_fastq = ch_fastq.mix(ch_fastq_fallback)
    ch_fastq.count().view{ count -> "Total (XSRA + fastq-dump) accession count: $count" }

    //-- group by reads by sample and join with star params --//
    ch_fastq = ch_fastq.groupTuple().join(ch_star_params)

    //-- Run STAR with the selected parameters on all reads --//
    // add the SA file size to the parameters channel
    ch_fastq = addSaSize(ch_fastq)

    // run STAR
    STAR_FULL(ch_fastq)

    // summarize the STAR results
    STAR_FULL_SUMMARY(
        STAR_FULL.out.summary
    )
}


process STAR_FULL_SUMMARY {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "star_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), path(summary_csv)

    output:
    tuple val(sample), path("combined.csv"), emit: "csv"
    path "${task.process}.log",              emit: "log"

    script:
    def use_database = params.use_database ? "--use-database" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    star-summary.py ${use_database} \\
      --sample ${sample} \\
      ${summary_csv} \\
      2>&1 | tee ${task.process}.log
    """
}

process STAR_FULL {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "star_env"
    maxRetries 3
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 12
    memory { starFullMem(sa_size, task.attempt) }
    time { 10.h * task.attempt }
    disk { [request: (375 * task.attempt).GB, type: 'local-ssd'] }
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*']
        return options[new Random().nextInt(options.size())]
    }

    input:
    tuple val(sample), path("input*_R1.fq.zst"), path("input*_R2.fq.zst"), 
          path(barcodes_file), path(star_index),
          val(cell_barcode_length), val(umi_length), val(strand), val(sa_size)

    output: 
    tuple val(sample), path("summary/*.csv"),                emit: summary
    tuple val(sample), path("h5ad/filtered/*.h5ad"),         emit: h5ad_filtered
    tuple val(sample), path("h5ad/raw/*.h5ad"),              emit: h5ad_raw, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.stats.gz"), emit: stats, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.txt.gz"),   emit: txt, optional: true
    path "${task.process}.log",                              emit: "log"

    script:
    def use_database = params.use_database ? "--use-database" : ""
    def keep_raw_h5ad = params.keep_raw_h5ad ? "--keep-raw-h5ad" : ""
    """
    echo "TEST"
    echo "# Running STAR for ${sample}" | tee ${task.process}.log

    # Format R1 and R2 file paths for STAR
    R1=\$(printf "%s," input*_R1.fq.zst)
    R1=\${R1%,} 
    R2=\$(printf "%s," input*_R2.fq.zst)
    R2=\${R2%,}
    
    # Define feature types as an array
    FEATURE_TYPES=("Gene" "GeneFull" "GeneFull_ExonOverIntron" "GeneFull_Ex50pAS" "Velocyto")
    FEATURE_STR=\$(printf "%s " "\${FEATURE_TYPES[@]}")
    
    # Run STAR
    STAR \\
      --readFilesIn \$R2 \$R1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${star_index} \\
      --soloCBwhitelist ${barcodes_file} \\
      --soloUMIlen ${umi_length} \\
      --soloStrand ${strand} \\
      --soloCBlen ${cell_barcode_length} \\
      --soloType CB_UMI_Simple \\
      --clipAdapterType CellRanger4 \\
      --outFilterScoreMin 30 \\
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
      --soloCellFilter EmptyDrops_CR \\
      --soloUMIfiltering MultiGeneUMI_CR \\
      --soloUMIdedup 1MM_CR \\
      --soloFeatures \${FEATURE_STR} \\
      --soloMultiMappers EM Uniform \\
      --outSAMtype None \\
      --soloBarcodeReadLength 0 \\
      --outFileNamePrefix results \\
      --readFilesCommand zstd -dcf \\
      2>&1 | tee -a ${task.process}.log

    echo "# Renaming the summary files" | tee -a ${task.process}.log
    mkdir -p summary/
    for feature in "\${FEATURE_TYPES[@]}"; do
        mv resultsSolo.out/\${feature}/Summary.csv summary/\${feature}.csv
    done

    echo "# Compressing the output for the sake of scanpy" | tee -a ${task.process}.log
    find resultsSolo.out -type f -name "*.mtx" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.tsv" | xargs -P ${task.cpus} gzip

    echo "# Converting mtx to h5ad" | tee -a ${task.process}.log
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"
    mtx-to-h5ad.py ${use_database} ${keep_raw_h5ad} \\
      --sample ${sample} \\
      resultsSolo.out 2>&1 | \\
      tee -a ${task.process}.log 
    """
}

def saveAsSTAR(sample, filename) {
    def extensions = [".h5ad", ".txt.gz", ".stats.gz", ".csv"]
    if (extensions.any { filename.endsWith(it) }) {
        def parts = filename.tokenize("/")
        if (parts.size() > 1) {
            parts = parts[1..-1]
        } 
        def org_part = null
        if (filename.endsWith(".csv")) {
            org_part = "summary"
        } 
        if (org_part != null) {
            if (parts.size() > 1) {
                parts = parts[0..-2] + [org_part] + [parts[-1]]
            } else {
                parts = [org_part] + parts
            }
        }   
        return "STAR/${sample}/" + parts.join('/')
    } 
    return null
}

process XSRA {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "download_env"
    maxRetries 3
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 6
    memory { 8.GB * (task.attempt > 2 ? task.attempt - 1 : 1) }
    disk { [request: (375 * (task.attempt > 2 ? task.attempt - 1 : 1)).GB, type: 'local-ssd'] }
    time { 2.h * task.attempt }
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accessions)

    output:
    tuple val(sample), path("reads/read_1.fq.zst"), emit: "R1"
    tuple val(sample), path("reads/read_2.fq.zst"), emit: "R2", optional: true
    path "${task.process}.log",                     emit: "log"

    script:
    def accessions_str = accessions.join(" ")
    def use_database = params.use_database ? "--use-database" : ""
    def gcp_project_id = params.gcp_project_id ? "--gcp-project-id ${params.gcp_project_id}" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    xsra-full.py ${use_database} ${gcp_project_id} \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --provider ${params.sra_provider} \\
      --output-dir reads \\
      ${accessions_str} \\
      2>&1 | tee ${task.process}.log
    """
}

process FASTQ_DUMP {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 1
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 4
    memory { 4.GB * task.attempt }
    time { 6.h * task.attempt }
    disk {[request: 375.GB, type: 'local-ssd']}
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accessions)

    output:
    tuple val(sample), path("reads/read_1.fq.zst"), emit: "R1"
    tuple val(sample), path("reads/read_2.fq.zst"), emit: "R2", optional: true
    path "${task.process}.log",                     emit: "log"

    script:
    def accessions_str = accessions.join(" ")
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    fastq-dump.py \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --outdir reads \\
      --maxSpotId ${params.fallback_max_spots} \\
      ${accessions_str} \\
      2>&1 | tee -a ${task.process}.log
    """
}
