def readAccessions(accessions_input){
    // Read the input CSV file with the sample names and SRA accessions
    ch_acc = accessions_input
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def req_columns = ["sample", "accession"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
                error "Missing columns in the input CSV file: ${miss_columns}"
            }
            // remove special characters from the sample name
            row.sample = row.sample.replaceAll("\\s", "_")
            def result = [row.sample, row.accession]
            // add optional, metadata columns
            def metadata = [:]
            ["organism", "tech_10x"].each { col ->
                metadata[col] = row.containsKey(col) ? row[col] : ""
            }
            result << metadata
            return result
        }

    // print srx values
    ch_acc
        .map{ sample, accession, metadata -> sample }
        .distinct()
        .collect() 
        .map{ it.join(',') }
        .view{ "SRX accessions: ${it}" } 

    return ch_acc
}

