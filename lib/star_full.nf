def addSaSize(ch_fastq) {
    // Add SA file size to the parameters channel
    def ch_fastq_with_size = ch_fastq.map { sample, fastq_1, fastq_2, barcodes_file, star_index, cell_barcode_length, umi_length, strand ->
        def sa_path = file("${star_index}/SA")
        def sa_size = 0L // Default size if file not found
        try {
            if (sa_path.exists()) {
                sa_size = sa_path.size()
            } else {
                log.warn "[${sample}/${accession}] SA file not found at ${sa_path}. Using default memory overhead (20 GB)."
            }
        } catch (Exception e) {
            log.warn "[${sample}/${accession}] Error accessing SA file ${sa_path}: ${e.getMessage()}. Using default memory overhead (20 GB)."
        }
        def size_gb = nextflow.util.MemoryUnit.of(sa_size)
        tuple(sample, fastq_1, fastq_2, barcodes_file, star_index, cell_barcode_length, umi_length, strand, sa_size)
    }
    return ch_fastq_with_size
}

// Calculate memory for STAR_FULL based on SA file size
def starFullMem(sa_size, task_attempt) {
    if (sa_size < 0){
        return 32.GB * task_attempt
    }
    def mem_gb = Math.round(nextflow.util.MemoryUnit.of(sa_size).toGiga())
    return mem_gb.GB * task_attempt
}