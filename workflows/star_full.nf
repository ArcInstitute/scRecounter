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

    // Use xsra to download all reads
    ch_fastq = XSRA(ch_accessions_filt)
    ch_fastq = joinReads(ch_fastq.R1, ch_fastq.R2)

    // For accessions lacking paired reads from fasterq-dump, fallback to fastq-dump
    ch_accessions_fallback = ch_accessions_filt
        .join(
            ch_fastq.map{ it -> [it[0], it[1], true] }, 
            by: [0,1],
            remainder: true
        )
        .filter{ it -> it[4] != true }
        .map{ it -> it[0..3] }

    // run fastq-dump on the fallback accessions
    ch_fastq_fallback = FASTQ_DUMP(ch_accessions_fallback)
    ch_fastq_fallback = joinReads(ch_fastq_fallback.R1, ch_fastq_fallback.R2)
    ch_fastq_fallback.count().view{ count -> "No. of fastq-dump fallback accessions: $count" }

    // combine the fasterq-dump and fastq-dump results
    ch_fastq = ch_fastq.mix(ch_fastq_fallback)
    ch_fastq.count().view{ count -> "No. of fast(er)q-dump accessions: $count" }

    // combine reads and star params
    ch_fastq = ch_fastq
        .map{ sample, accession, metadata, R1, R2 -> [sample, R1, R2] }
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
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    disk { [request: (375 * (task.attempt > 1 ? 2 : 1)).GB, type: 'local-ssd'] }
    machineType { 
        def options = ['n2-*', 'n2d-*']
        return options[new Random().nextInt(options.size())]
    }

    input:
    tuple val(sample), path("input*_R1.fastq"), path("input*_R2.fastq"), 
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

    R1=\$(printf "%s," input*_R1.fastq)
    R1=\${R1%,} 
    R2=\$(printf "%s," input*_R2.fastq)
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

process FASTQ_DUMP {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 1
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 4
    memory { 4.GB * task.attempt }
    time { (6.h + (sra_file_size_gb * 0.8).h) * task.attempt }
    disk {[request: 375.GB, type: 'local-ssd']}
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "${task.process}.log",                                                   emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    echo "Downloading ${accession} for ${sample}" > ${task.process}.log
    echo "sra-stat file size: ${sra_file_size_gb} GB" >> ${task.process}.log

    echo "Running fastq-dump as backup for fasterq-dump" >> ${task.process}.log
    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --outdir reads \\
      --maxSpotId ${params.fallback_max_spots} \\
      ${accession} \\
      2>&1 | tee -a ${task.process}.log
    """

    stub:
    """
    mkdir -p reads
    touch reads/read_1.fastq reads/read_2.fastq ${task.process}.log
    """
}

process XSRA {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 1
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 6
    memory { 16.GB * task.attempt }
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
        [request: disk_size, type: 'local-ssd'] 
    }
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1"
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "${task.process}.log",                                                   emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    echo "Downloading ${accession} for ${sample}" > ${task.process}.log
    echo "sra-stat file size: ${sra_file_size_gb} GB" >> ${task.process}.log

    xsra.py \\
      --sample ${sample} \\
      --threads ${task.cpus} \\
      --min-read-length ${params.min_read_len} \\
      --outdir reads \\
      ${accession}
    """
}

process FASTERQ_DUMP {
    publishDir file(params.output_dir), mode: "copy", overwrite: true, saveAs: { filename -> saveAsLog(filename, sample, accession) }
    label "download_env"
    maxRetries 1
    errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
    cpus 6
    memory { 16.GB * task.attempt }
    time { (10.h + (sra_file_size_gb * 0.8).h) * task.attempt }
    disk { 
        def disk_size = 
            sra_file_size_gb > 260 ? 2625.GB :
            sra_file_size_gb > 220 ? 2250.GB :
            sra_file_size_gb > 170 ? 1875.GB :
            sra_file_size_gb > 120 ? 1500.GB :
            sra_file_size_gb > 60 ? 1125.GB :
            sra_file_size_gb > 30 ? 750.GB :
            375.GB
        disk_size = disk_size + (375 * (task.attempt - 1)).GB
        [request: disk_size, type: 'local-ssd'] 
    }
    machineType { 
        def options = ['n2-*', 'c2-*', 'n2d-*', 'c2d-*'] //, 'n1-*']
        return options[new Random().nextInt(options.size())]
    }
    
    input:
    tuple val(sample), val(accession), val(metadata), val(sra_file_size_gb)

    output:
    tuple val(sample), val(accession), val(metadata), path("reads/read_1.fastq"), emit: "R1", optional: true
    tuple val(sample), val(accession), val(metadata), path("reads/read_2.fastq"), emit: "R2", optional: true
    path "${task.process}.log",                                                   emit: "log"

    script:
    """
    export GCP_SQL_DB_HOST="${params.db_host}"
    export GCP_SQL_DB_NAME="${params.db_name}"
    export GCP_SQL_DB_USERNAME="${params.db_username}"

    echo "Downloading ${accession} for ${sample}" > ${task.process}.log
    echo "sra-stat file size: ${sra_file_size_gb} GB" >> ${task.process}.log

    # run prefetch and fasterq-dump
    fq-dump.py \\
      --sample ${sample} \\
      --accession ${accession} \\
      --threads ${task.cpus} \\
      --bufsize 200MB \\
      --curcache 1GB \\
      --mem 12GB \\
      --max-size-gb ${params.max_sra_size} \\
      --min-read-length ${params.min_read_len} \\
      --temp ${params.fasterq_tmp} \\
      --outdir reads \\
      ${accession} \\
      2>&1 | tee -a ${task.process}.log
    """

    stub:
    """
    mkdir -p reads
    touch reads/read_1.fastq reads/read_2.fastq ${task.process}.log
    """
}
