def addSaSize(ch_fastq) {
    // Add SA file size to the parameters channel
    def ch_fastq_with_size = ch_fastq.map { sample, fastq_1, fastq_2, barcodes_file, star_index, cell_barcode_length, umi_length, strand ->
        def sa_path = file("${star_index}/SA")
        def sa_size = 30L * 1024 * 1024 * 1024 // Default size if file not found
        try {
            if (sa_path.exists()) {
                sa_size = sa_path.size()
            } else {
                log.warn "${sample} SA file not found at ${sa_path}. Using default memory overhead (30 GB)."
            }
        } catch (Exception e) {
            log.warn "${sample} Error accessing SA file ${sa_path}: ${e.getMessage()}. Using default memory overhead (30 GB)."
        }
        sa_size = Math.round(nextflow.util.MemoryUnit.of(sa_size).toGiga()).GB
        tuple(sample, fastq_1, fastq_2, barcodes_file, star_index, cell_barcode_length, umi_length, strand, sa_size)
    }
    return ch_fastq_with_size
}
