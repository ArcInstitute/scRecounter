// Workflow to run STAR alignment on scRNA-seq data
workflow STAR_WF{
    take:
    ch_fastq

    main:
    // for each sample (accession), run STAR on subset of reads with various parameters to determine which params produce the most valid barcodes
    
    // load barcodes file
    ch_barcodes = Channel
        .fromPath(params.barcodes, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def req_columns = ["name", "cell_barcode_length", "umi_length", "file_path"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            // remove special characters
            row.name = row.name.replaceAll("\\s", "_")
            return [
                row.name, 
                row.cell_barcode_length.toInteger(), 
                row.umi_length.toInteger(),
                file(row.file_path)  
            ]
        }

    // load star indices
    ch_star_indices = Channel
        .fromPath(params.star_indices, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def req_columns = ["organism", "star_index"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            // remove special characters
            row.organism = row.organism.replaceAll("\\s", "_")
            return [row.organism, file(row.star_index)]
        }

    // Subsample reads
    SUBSAMPLE_READS(ch_fastq)

    // Get read lengths
    SEQKIT_STATS(SUBSAMPLE_READS.out)

    // Pairwise combine samples with barcodes, strand, and star index
    ch_fastq_barcodes = SUBSAMPLE_READS.out
        .combine(Channel.of("Forward", "Reverse"))
        .combine(ch_barcodes)
        .combine(ch_star_indices)

    // Run STAR on subsampled reads, for all pairwise parameter combinations
    STAR_GET_VALID_BARCODES(ch_fastq_barcodes)

    // Set STAR parameters
    ch_valid_barcodes = STAR_GET_VALID_BARCODES.out
        .groupTuple()
        .join(SEQKIT_STATS.out, by: 0)
    STAR_SET_PARAMS(
        ch_valid_barcodes, 
        Channel.fromPath(params.barcodes, checkIfExists: true),
        Channel.fromPath(params.star_indices, checkIfExists: true)
    )

    // Run STAR with best parameters
    if (! params.define){
        STAR_FULL(ch_fastq.join(STAR_SET_PARAMS.out, by: 0))
    }
}

// STAR alignment with all reads and selected parameters
def saveAsSTAR(sample, filename) {
    if (filename.endsWith(".mtx") || filename.endsWith(".tsv")){
        def parts = filename.tokenize("/")
        return "${sample}/" + parts[1..-1].join('/')
    } 
    return null
}

process STAR_FULL {
    publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true, saveAs: { filename -> saveAsSTAR(sample, filename) }
    conda "envs/star.yml"
    label "process_high"

    input:
    tuple val(sample), path("input*_R1.fastq"), path("input*_R2.fastq"), path(star_params)

    output:
    tuple val(sample), path("resultsSolo.out/Gene/raw/*"),                         emit: gene_raw
    tuple val(sample), path("resultsSolo.out/Gene/filtered/*"),                    emit: gene_filt, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull/raw/*"),                     emit: gene_full_raw
    tuple val(sample), path("resultsSolo.out/GeneFull/filtered/*"),                emit: gene_full_filt, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/raw/*"),             emit: gene_ex50_raw
    tuple val(sample), path("resultsSolo.out/GeneFull_Ex50pAS/filtered/*"),        emit: gene_ex50_filt, optional: true
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/raw/*"),      emit: gene_ex_int_raw
    tuple val(sample), path("resultsSolo.out/GeneFull_ExonOverIntron/filtered/*"), emit: gene_ex_int_filt, optional: true
    tuple val(sample), path("resultsSolo.out/Velocyto/raw/*"),                     emit: velocyto_raw
    tuple val(sample), path("resultsSolo.out/Velocyto/filtered/*"),                emit: velocyto_filt, optional: true

    script:
    """
    # load parameters
    json2env.py \\
      --params BARCODE_FILE CELL_BARCODE_LENGTH UMI_LENGTH STRAND STAR_INDEX \\
      -- $star_params > params.env
    source params.env

    # run STAR
    R1=\$(printf "%s," input*_R1.fastq)
    R1=\${R1%,} 
    R2=\$(printf "%s," input*_R2.fastq)
    R2=\${R2%,}
    STAR \\
      --readFilesIn \$R2 \$R1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir \$STAR_INDEX \\
      --soloCBwhitelist \$BARCODE_FILE \\
      --soloUMIlen \$UMI_LENGTH \\
      --soloStrand \$STRAND \\
      --soloCBlen \$CELL_BARCODE_LENGTH \\
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
      --outFileNamePrefix results
    """
}
//TODO: remove `--soloBarcodeReadLength 0`?

// Set STAR parameters based on valid barcodes
def saveAsParams(sample, filename) {
    return "${sample}/${filename}"
}

process STAR_SET_PARAMS {
    publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true, saveAs: { filename -> saveAsParams(sample, filename) }
    conda "envs/star.yml"

    input:
    tuple val(sample), path(star_summary), path(stats_tsv)
    each path(barcodes)
    each path(star_index)

    output:
    tuple val(sample), path("star_params.json")

    script:
    """
    set_star_params.py \\
      --sample $sample \\
      --stats $stats_tsv \\
      --barcodes $barcodes \\
      --star-index $star_index \\
      $star_summary
    """

    stub:
    """
    touch star_params.json
    """
}

// Set STAR parameters based on valid barcodes
def saveAsValid(sample, filename) {
    return "${sample}/param_search/${filename}"
}

// Run STAR alignment on subsampled reads with various parameters to determine which parameters produce the most valid barcodes
process STAR_GET_VALID_BARCODES {
    publishDir file(params.outdir) / "STAR", mode: "copy", overwrite: true, saveAs: { filename -> saveAsValid(sample, filename) }
    conda "envs/star.yml"
    label "process_medium"

    input:
    tuple val(sample), path(fastq_1), path(fastq_2), 
          val(strand), val(barcode_name), val(cell_barcode_length), val(umi_length), path(barcodes), 
          val(organism), path(star_index)

    output:
    tuple val(sample), path("${sample}_${organism}_${strand}_${barcode_name}_star-summary.csv")

    script:
    """
    # run STAR
    STAR \\
      --readFilesIn $fastq_2 $fastq_1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${star_index} \\
      --soloCBwhitelist $barcodes \\
      --soloCBlen ${cell_barcode_length} \\
      --soloUMIlen ${umi_length} \\
      --soloStrand ${strand} \\
      --soloType CB_UMI_Simple \\
      --clipAdapterType CellRanger4 \\
      --outFilterScoreMin 30 \\
      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \\
      --soloCellFilter EmptyDrops_CR \\
      --soloUMIfiltering MultiGeneUMI_CR \\
      --soloUMIdedup 1MM_CR \\
      --soloFeatures GeneFull \\
      --soloMultiMappers EM \\
      --outSAMtype None \\
      --soloBarcodeReadLength 0 \\
      --outFileNamePrefix results 

    # format output
    OUTNAME="${sample}_${organism}_${strand}_${barcode_name}_star-summary.csv"
    mv -f resultsSolo.out/GeneFull/Summary.csv \$OUTNAME
    echo "STRAND,${strand}" >> \$OUTNAME
    echo "BARCODE_NAME,${barcode_name}" >> \$OUTNAME
    echo "CELL_BARCODE_LENGTH,${cell_barcode_length}" >> \$OUTNAME
    echo "UMI_LENGTH,${umi_length}" >> \$OUTNAME
    echo "ORGANISM,${organism}" >> \$OUTNAME
    """

    stub:
    """
    touch ${sample}_${organism}_${strand}_${barcode_name}_star-summary.csv
    """
}

// Get read lengths via `seqkit stats`
process SEQKIT_STATS {
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)

    output:
    tuple val(sample), path("${sample}_stats.tsv")

    script:
    """
    seqkit -j $task.cpus stats -T $fastq_1 $fastq_2 > ${sample}_stats.tsv
    """

    stub:
    """
    touch ${sample}_stats.tsv
    """
}

// Subsample reads
process SUBSAMPLE_READS {
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), path("input*_R1.fq"), path("input*_R2.fq")

    output:
    tuple val(sample), path("${sample}_R1.fq"), path("${sample}_R2.fq")

    script: 
    """
    subsample.py --num-seqs ${params.subsample} --out-file ${sample}_R1.fq input*_R1.fq
    subsample.py --num-seqs ${params.subsample} --out-file ${sample}_R2.fq input*_R2.fq
    """
    
    stub:
    """
    touch ${sample}_R1.fq ${sample}_R2.fq
    """
}
