workflow CHARACTERIZE_WF{
    take:
    ch_fastq

    main:
    // for each barcode, run STAR & determine the number of reads aligned
    //  align two reads against three possible chemistry configurations, v3, multiome GEX, and v2/5'.  
    
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
            return [
                row.name, 
                row.cell_barcode_length.toInteger(), 
                row.umi_length.toInteger(),
                file(row.file_path)  
            ]
        }

    // subsample reads
    SUBSAMPLE_READS(ch_fastq)

    // get read lengths
    SEQKIT_STATS(SUBSAMPLE_READS.out)

    // pairwise combine samples with barcodes and strand
    ch_fastq_barcodes = SUBSAMPLE_READS.out
        .combine(Channel.of("Forward", "Reverse"))
        .combine(ch_barcodes)
        
    // run STAR on subsampled reads, for all pairwise parameter combinations
    STAR_GET_VALID_BARCODES(ch_fastq_barcodes)

    // Set STAR parameters
    ch_valid_barcodes = STAR_GET_VALID_BARCODES.out
        .groupTuple()
        .join(SEQKIT_STATS.out, by: 0)

    SET_STAR_PARAMS(ch_valid_barcodes)
}

process SET_STAR_PARAMS {
    conda "envs/star.yml"

    input:
    tuple val(sample), path(summary_csv), path(stats_tsv)

    output:
    tuple val(sample), path("star_params.json")

    script:
    """
    set_star_params.py \\
      --sample $sample \\
      --stats $stats_tsv \\
      $summary_csv > star_params.json
    """

    stub:
    """
    touch star_params.json
    """
}

process STAR_GET_VALID_BARCODES {
    conda "envs/star.yml"
    label "process_medium"
    //scratch true

    input:
    tuple val(sample), path(fastq_1), path(fastq_2), val(strand), val(barcode_name), val(cell_barcode_length), val(umi_length), path(barcodes)

    output:
    tuple val(sample), path("${sample}_${strand}_${barcode_name}.csv")

    script:
    """
    # run STAR
    STAR \\
      --readFilesIn $fastq_2 $fastq_1 \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${params.star_index} \\
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
    OUTNAME="${sample}_${strand}_${barcode_name}.csv"
    mv -f resultsSolo.out/GeneFull/Summary.csv \$OUTNAME
    echo "STRAND,${strand}" >> \$OUTNAME
    echo "BARCODE_NAME,${barcode_name}" >> \$OUTNAME
    echo "CELL_BARCODE_LENGTH,${cell_barcode_length}" >> \$OUTNAME
    echo "UMI_LENGTH,${umi_length}" >> \$OUTNAME
    """

    stub:
    """
    touch ${sample}_${strand}_${barcode_name}.csv
    """
}

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

process SUBSAMPLE_READS {
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), path("R1.fq"), path("R2.fq")

    output:
    tuple val(sample), path("${sample}_R1.fq"), path("${sample}_R2.fq")

    script:
    """
    # subsample to 0.1 mil reads
    head -n 400000 R1.fq > ${sample}_R1.fq
    head -n 400000 R2.fq > ${sample}_R2.fq
    """
    
    stub:
    """
    touch ${sample}_R1.fq ${sample}_R2.fq
    """
}
