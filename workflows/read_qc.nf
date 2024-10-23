// Read QC workflow for fastq files
workflow READ_QC_WF{
    take:
    ch_fastq

    main:
    // Flatten the channel to process each read file separately
    ch_fastq_flat = ch_fastq.flatMap { sample, fastq_1, fastq_2 -> 
        [ [sample, "R1", fastq_1], [sample, "R2", fastq_2] ] 
    }

    // Run seqkit stats
    SEQKIT_STATS(ch_fastq_flat)
        .collectFile(
          name: "seqkit-stats.tsv", 
          storeDir:  file(params.outdir) / "read_qc", 
          newLine: false, keepHeader: true
        )
}

// Run `seqkit stats` on fastq files
process SEQKIT_STATS {
    conda "envs/read_qc.yml"
    label "process_low"

    input:
    tuple val(sample), val(read), path("${sample}_${read}_*.fastq")

    output:
    path "${sample}_${read}.tsv"

    script:
    """
    seqkit -j $task.cpus stats -a -T ${sample}_${read}_*.fastq > ${sample}_${read}.tsv
    """

    stub:
    """
    touch ${sample}_${read}.tsv
    """
}
