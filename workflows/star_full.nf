include { joinReads; saveAsLog; } from '../lib/utils.groovy'

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

    // combine reads and star params
    ch_fastq = XSRA.out.R1
        .map{ sample,fastq -> [sample,fastq] }
        .join(
            XSRA.out.R2.map{ sample,fastq -> [sample,fastq] }
        )
        .groupTuple()
        .join(ch_star_params)

    //-- Run STAR with the selected parameters on all reads --//
    // run STAR
    STAR_FULL(ch_fastq)

    // summarize the STAR results
    STAR_FULL_SUMMARY(
        STAR_FULL.out.gene_summary,
        STAR_FULL.out.gene_full_summary,
        STAR_FULL.out.gene_ex50_summary,
        STAR_FULL.out.gene_ex_int_summary,
        STAR_FULL.out.velocyto_summary
    )
}


process STAR_FULL_SUMMARY {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "star_env"
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk 10.GB

    input:
    tuple val(sample), path("gene_summary.csv")
    tuple val(sample), path("gene_full_summary.csv")
    tuple val(sample), path("gene_ex50_summary.csv")
    tuple val(sample), path("gene_ex_int_summary.csv")
    tuple val(sample), path("velocyto_summary.csv")

    output:
    tuple val(sample), path("Summary.csv"), emit: "csv"
    path "${task.process}.log",             emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    star-summary.py \\
      --sample ${sample} \\
      gene_summary.csv \\
      gene_full_summary.csv \\
      gene_ex50_summary.csv \\
      gene_ex_int_summary.csv \\
      velocyto_summary.csv \\
      2>&1 | tee ${task.process}.log
    """
}

process STAR_FULL {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, filename) }
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "star_env"
    label "process_high"
    //errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk { [request: (375 * (task.attempt > 1 ? 2 : 1)).GB, type: 'local-ssd'] }
    machineType { 
        def options = ['n2-*', 'n2d-*']
        return options[new Random().nextInt(options.size())]
    }

    input:
    tuple val(sample), path("input*_R1.fq.zst"), path("input*_R2.fq.zst"), 
          path(barcodes_file), path(star_index),
          val(cell_barcode_length), val(umi_length), val(strand)

    output: 
    tuple val(sample), path("resultsSolo.out/Gene/Summary.csv"),                    emit: gene_summary
    tuple val(sample), path("resultsSolo.out/GeneFull/Summary.csv"),                emit: gene_full_summary
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/Summary.csv"),        emit: gene_ex50_summary
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/Summary.csv"), emit: gene_ex_int_summary
    tuple val(sample), path("resultsSolo.out/Velocyto/Summary.csv"),                emit: velocyto_summary
    tuple val(sample), path("resultsSolo.out/*/raw/*"),                             emit: raw
    tuple val(sample), path("resultsSolo.out/*/filtered/*"),                        emit: filt, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.stats.gz"),                        emit: stats, optional: true
    tuple val(sample), path("resultsSolo.out/*/*.txt.gz"),                          emit: txt, optional: true
    path "${task.process}.log",                                                     emit: "log"

    script:
    """
    echo "Running STAR for ${sample}" > ${task.process}.log

    R1=\$(printf "%s," input*_R1.fq.zst)
    R1=\${R1%,} 
    R2=\$(printf "%s," input*_R2.fq.zst)
    R2=\${R2%,}
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
      --soloFeatures Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS Velocyto \\
      --soloMultiMappers EM Uniform \\
      --outSAMtype None \\
      --soloBarcodeReadLength 0 \\
      --outFileNamePrefix results \\
      --readFilesCommand zstd -dcf \\
      2>&1 | tee -a ${task.process}.log

    # gzip the results
    mkdir -p resultsSolo.out
    find resultsSolo.out -type f -name "*.stats" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.txt" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.tsv" | xargs -P ${task.cpus} gzip
    find resultsSolo.out -type f -name "*.mtx" | xargs -P ${task.cpus} gzip
    """
}

def saveAsSTAR(sample, filename) {
    def extensions = [".mtx.gz", ".tsv.gz", ".txt.gz", ".stats.gz", ".csv"]
    if (extensions.any { filename.endsWith(it) }) {
        def parts = filename.tokenize("/")
        if (parts.size() > 1) {
            return "STAR/${sample}/" + parts[1..-1].join('/')
        } else {
            return "STAR/${sample}/" + parts[0]
        }
    } 
    return null
}

process XSRA {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample) }
    label "download_env"
    maxRetries 1
    //errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }   // TODO: uncomment this
    cpus 6
    memory { 12.GB * task.attempt }
    /*
    time { (4.h + (sra_file_size_gb * 0.8).h) * task.attempt }
    disk { 
        def disk_size = 
            sra_file_size_gb > 360 ? 7 * 375.GB :
            sra_file_size_gb > 300 ? 6 * 375.GB :
            sra_file_size_gb > 240 ? 5 * 375.GB :
            sra_file_size_gb > 180 ? 4 * 375.GB :
            sra_file_size_gb > 120 ? 3 * 375.GB :
            sra_file_size_gb > 50 ? 2 * 375.GB :
            375.GB
        disk_size = disk_size + (375 * (task.attempt - 1)).GB
        [re
        quest: disk_size, type: 'local-ssd'] 
    }
    */
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
    accessions = accessions.join(" ")
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"
    export GCP_PROJECT_ID="${params.gcp_project_id}"

    xsra-full.py \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --output-dir reads \\
      ${accessions} \\
      2>&1 | tee ${task.process}.log
    """
}

