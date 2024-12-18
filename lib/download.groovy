def readAccessions(accessions_file){
    return Channel
        .fromPath(accessions_file, checkIfExists: true)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def req_columns = ["sample", "accession"]
            def opt_columns = ["organism"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            // remove special characters from the sample name
            row.sample = row.sample.replaceAll("\\s", "_")
            def result = [row.sample, row.accession]
            if (row.containsKey("organism")) {
                result << row.organism
            }
            return result
        }
}

def joinReads(ch_read1, ch_read2){
    return ch_read1.join(ch_read2, by: [0, 1], remainder: true)
        .filter { sample, accession, r1, r2 -> 
            if(r2 == null) {
                println "Warning: Read 2 is empty for ${sample}; skipping"
            }
            return r2 != null
        }
        .groupTuple()
        .map { sample, accession, fastq_1, fastq_2 ->
            return [sample, fastq_1.flatten(), fastq_2.flatten()]
        }
}