include { joinReads; saveAsLog; subsampleByGroup; } from '../lib/utils.nf'
include { makeParamSets; validateRequiredColumns; loadBarcodes; loadStarIndices; expandStarParams; addSaSizeToParams } from '../lib/star_params.nf'

// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_PARAMS_WF{
    take:
    ch_accessions
    ch_sra_stat

    main:
    //-- Download subset of reads reads --//    
    // Run `xsra dump --limit`; only testing on subset of SRR accessions
    XSRA(
        subsampleByGroup(ch_accessions, params.max_accessions, 23482)
    )
    
    // Join R1 and R2 channels, which will filter out empty R2 records
    ch_fasta = joinReads(XSRA.out.R1, XSRA.out.R2)

    // Get read lengths
    SEQKIT_STATS(ch_fasta)

    //-- STAR param search on subsampled reads --//

    // Load barcodes file
    ch_barcodes = loadBarcodes(params)

    // Load star indices
    ch_star_indices = loadStarIndices(params)

    // Pairwise combine samples with barcodes, strand, and star index
    ch_params = makeParamSets(ch_fasta, ch_barcodes, ch_star_indices)

    // Add SA file size to the parameters channel
    ch_params = addSaSizeToParams(ch_params)

    // Run STAR on subsampled reads, for all pairwise parameter combinations
    STAR_PARAM_SEARCH(ch_params)

    // remove the last value of ch_params (sa_size)
    ch_params = ch_params.map { tuple -> tuple[0..-2] }

    // Format the STAR parameters into a CSV file
    STAR_FORMAT_PARAMS(STAR_PARAM_SEARCH.out.csv)
    
    // Get best parameters
    ch_params_all = STAR_FORMAT_PARAMS.out.csv
        .groupTuple(by: [0,1])
        .join(SEQKIT_STATS.out, by: [0,1])
        .join(ch_sra_stat, by: [0,1])
    STAR_SELECT_PARAMS(ch_params_all)

    // Filter empty params
    ch_star_params_json = STAR_SELECT_PARAMS.out.json
        .filter { sample, accession, json_file -> 
            if(json_file.size() < 5) {
                log.warn "WARNING: No valid STAR parameters found for ${sample}; skipping"
            }
            return json_file.size() > 5
        }

    // Extract selected parameters from the JSON files
    ch_star_params = expandStarParams(ch_fasta, ch_star_params_json)

    // Merge STAR parameters to a single set for each sample
    def majorityRule = { list ->
        list.groupBy { it }
            .collectEntries { k, v -> [(k): v.size()] }
            .max { it.value }
            .key
    }

    ch_star_params = ch_star_params.groupTuple()
        .map { sample, accessions, barcodes, star_indices, cell_barcode_lengths, umi_lengths, strands ->
            def combined_params = [barcodes, star_indices, cell_barcode_lengths, umi_lengths, strands].transpose()
            [sample] + majorityRule(combined_params)
        }

    ch_star_params.view{ sample,barcodes_file,star_index,cb_len,umi_len,strand -> 
        def barcodes_name = barcodes_file.tokenize("/").last().tokenize(".")[0]
        def star_index_name = star_index.tokenize("/").last().tokenize(".")[0]
        "Selected parameter set for ${sample}:\n - barcodes: ${barcodes_name}\n - STAR index: ${star_index_name}\n - Barcode length: ${cb_len}\n - UMI length: ${umi_len}\n - Strand: ${strand}\n"
    }

    // Save the final STAR parameters to the database
    STAR_SAVE_FINAL_PARAMS(ch_star_params)

    emit:
    star_params = ch_star_params
}

// Set STAR parameters based on valid barcodes
def saveAsFinalParams(sample, filename) {
    if (filename.endsWith(".csv") || filename.endsWith(".json")){
        filename = filename.tokenize("/").last()
        return "STAR/${sample}/${filename}"
    }
    return null
}

process STAR_SAVE_FINAL_PARAMS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsFinalParams(sample, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "star_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), val(barcodes), val(star_index), val(cell_barcode_length), val(umi_length), val(strand)

    output:
    path "star_params.csv",     emit: "csv"
    path "${task.process}.log", emit: "log"

    script:
    def use_database = params.use_database ? "--use-database" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    upload-final-star-params.py ${use_database} \\
      --sample ${sample} \\
      --barcodes ${barcodes} \\
      --star-index ${star_index} \\
      --cell-barcode-length ${cell_barcode_length} \\
      --umi-length ${umi_length} \\
      --strand ${strand} \\
      2>&1 | tee ${task.process}.log
    """

    stub:
    """
    touch star_params.csv ${task.process}.log
    """
}

// Set STAR parameters based on valid barcodes
def saveAsParams(sample, accession, filename) {
    if (filename.endsWith(".csv") || filename.endsWith(".json")){
        filename = filename.tokenize("/").last()
        return "STAR/${sample}/${accession}/${filename}"
    }
    return null
}

process STAR_SELECT_PARAMS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsParams(sample, accession, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "star_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), val(accession), path("star_params*.csv"), path(read_stats), path(sra_stats)

    output:
    tuple val(sample), val(accession), path("results/merged_star_params.csv"),    emit: "csv"
    tuple val(sample), val(accession), path("results/selected_star_params.json"), emit: "json"
    path "${task.process}.log",  emit: "log"
    
    script:
    def use_database = params.use_database ? "--use-database" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"
    
    select-star-params.py ${use_database} \\
      --sample ${sample} \\
      --accession ${accession} \\
      $read_stats $sra_stats star_params*.csv \\
      2>&1 | tee ${task.process}.log
    """

    stub:
    """
    touch star_params.json merged_star_params.csv ${task.process}.log
    """
}

process STAR_FORMAT_PARAMS {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "star_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), val(accession), val(metadata), val(params), path(star_summary)

    output:
    tuple val(sample), val(accession), path("star_params.csv"), emit: "csv"
    path "${task.process}*.log",                                emit: "log"

    script:
    """
    BARCODES_FILE=\$(basename ${params.barcodes_file})
    STAR_INDEX=\$(basename ${params.star_index})

    format-star-params.py \\
      --sample ${params.sample} \\
      --accession ${params.accession} \\
      --strand ${params.strand} \\
      --barcodes-name ${params.barcodes_name} \\
      --barcodes-file ${params.barcodes_file} \\
      --cell-barcode-length ${params.cell_barcode_length} \\
      --umi-length ${params.umi_length} \\
      --organism ${params.organism} \\
      --star-index ${params.star_index} \\
      --outfile star_params.csv \\
      $star_summary \\
      2>&1 | tee "${task.process}:\${STAR_INDEX}:\${BARCODES_FILE}:${params.strand}.log"
    """
}

// Set STAR parameters based on valid barcodes
def saveAsValid(sample, filename) {
    return "param_search/${filename}"
}

// Run STAR alignment on subsampled reads with various parameters to determine which parameters produce the most valid barcodes
process STAR_PARAM_SEARCH {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "star_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 8
    memory { (sa_size > 0 ? Math.round(nextflow.util.MemoryUnit.of(sa_size).toGiga()) : 30.GB) * task.attempt }
    time { 4.h * task.attempt }
    disk 10.GB

    input:
    tuple val(sample), val(accession), val(metadata), path(fastq_1), path(fastq_2), 
          path(barcodes_file), path(star_index), val(params), val(sa_size)

    output:
    tuple val(sample), val(accession), val(metadata), val(params), path("star_summary.csv"), emit: "csv"
    path "${task.process}*.log", emit: "log"

    script:
    """
    BARCODES_FILE=\$(basename "${params.barcodes_file}")
    STAR_INDEX=\$(basename "${params.star_index}")

    # run STAR
    STAR \\
      --readFilesIn $fastq_2 $fastq_1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${star_index} \\
      --soloCBwhitelist ${barcodes_file} \\
      --soloCBlen ${params.cell_barcode_length} \\
      --soloUMIlen ${params.umi_length} \\
      --soloStrand ${params.strand} \\
      --soloType CB_UMI_Simple \\
      --outFilterScoreMin 30 \\
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
      --soloCellFilter EmptyDrops_CR \\
      --soloUMIfiltering MultiGeneUMI_CR \\
      --soloUMIdedup 1MM_CR \\
      --soloFeatures GeneFull \\
      --soloMultiMappers EM \\
      --outSAMtype None \\
      --soloBarcodeReadLength 0 \\
      --outFileNamePrefix results \\
      --readFilesCommand zstd -dcf \\
      2>&1 | tee "${task.process}:\${STAR_INDEX}:\${BARCODES_FILE}:${params.strand}.log"
    
    # rename output
    mv -f resultsSolo.out/GeneFull/Summary.csv star_summary.csv
    """

    stub:
    """
    touch star_summary.csv ${task.process}.log
    """
}

// Get read lengths via `seqkit stats`
process SEQKIT_STATS {
    label "download_env"
    disk 10.GB

    input:
    tuple val(sample), val(accession), val(metadata), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), val(accession), path("${sample}_${accession}_stats.tsv")

    script:
    """
    seqkit stats -T $fastq_1 $fastq_2 \\
      > ${sample}_${accession}_stats.tsv
    """

    stub:
    """
    touch ${sample}_${accession}_stats.tsv
    """
}

// Dump a subset of reads from the SRA
process XSRA {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 3
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 4
    memory { 4.GB * (task.attempt > 2 ? task.attempt - 1 : 1) }
    disk { 10.GB * (task.attempt > 2 ? task.attempt - 1 : 1) }
 
    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fa.zst"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fa.zst"), emit: "R2", optional: true
    path "${task.process}.log",                                                    emit: "log"

    script:
    def use_database = params.use_database ? "--use-database" : ""
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    xsra-limit.py ${use_database} \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --limit ${params.max_spots} \\
      --output-dir reads \\
      ${accession} \\
      2>&1 | tee ${task.process}.log
    """
}
