def readStarParams(star_params_file){
    // read the input CSV file and check if all required columns are present
    return Channel
        .fromPath(star_params_file, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def req_columns = ["sample", "fastq_1", "fastq_2", "barcodes_file", "star_index", 
                                "cell_barcode_length", "umi_length", "strand"]
            def miss_columns = req_columns.findAll { !row.containsKey(it) }
            if (miss_columns) {
            error "Missing columns in the input CSV file: ${miss_columns}"
        }
        // remove special characters from the sample name
        row.sample = row.sample.replaceAll("\\s", "_")
        return [row.sample, row.fastq_1, row.fastq_2, row.barcodes_file, row.star_index, 
                row.cell_barcode_length, row.umi_length, row.strand]
    }
}

def readAccessions(accessions_input){
    // read the input accessions CSV file and check if all required columns are present
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
                metadata[col] = row.containsKey(col) ? row[col].replaceAll("\\s", "_") : ""
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

def addStats(ch_accessions, ch_sra_stat){
    // add file size information to the accessions
    ch_stats = ch_sra_stat
        .map{ sample,acc,csv -> csv }
        .splitCsv(header: true, sep: ",")
        .map{ row -> [row.accession, row.file_size_gb.toDouble()] }
        .join(ch_sra_stat.map{ sample,acc,csv -> [acc, sample] }, by: [0])
        .map{ acc, size, sample -> [sample, acc, size] }
    return ch_accessions.join(ch_stats, by: [0,1]) // sample, acc, metadata, size
}

def joinReads(ch_read1, ch_read2){
    // extract metadata to prevent incorrect joining
    ch_metadata = ch_read1.map{ sample,accession,metadata,fastq -> [sample,accession,metadata] }

    // join the read1 and read2 channels
    return ch_read1
        .map{ sample,accession,metadata,fastq -> [sample,accession,fastq] }
        .join(
            ch_read2.map{ sample,accession,metadata,fastq -> [sample,accession,fastq] }, 
            by: [0,1]
        )
        .join(
            ch_metadata, by: [0,1]
        )
        .map{ 
            sample,accession,fastq1,fastq2,metadata -> [sample,accession,metadata,fastq1,fastq2] 
        }
}

def saveAsLog(filename, sample=null, accession=null) {
    if (filename.endsWith(".log")) {
        def basename = filename.tokenize("/")[-1]
        def path = "logs"
        if (sample){
            path = "${path}/${sample}"
        }
        if (accession) {
            path = "${path}/${accession}"
        } 
        path = "${path}/${basename}"
        return path
    }
    return null
}

def subsampleByGroup(ch_accessions, max_per_group, seed) {
    ch_accessions
        .groupTuple()
        .map { samples, accessions, meta, sra_stat ->
            accessions = accessions.toList()
            meta = meta.toList()
            sra_stat = sra_stat.toList()
        
            if (accessions) { // Ensure lists are not empty
                def indices = (0..<accessions.size()).toList() // Create indices
                indices.shuffle(new Random(seed)) // Shuffle indices with seed
            
                def shuffledAcc = indices.collect { accessions[it] } // Shuffle accessions based on indices
                def shuffledMeta = indices.collect { meta[it] } // Shuffle meta based on the same indices
                def shuffledStat = indices.collect { sra_stat[it] } // Shuffle meta based on the same indices
            
                def maxSize = Math.min(shuffledAcc.size(), max_per_group)
                shuffledAcc = maxSize > 0 ? shuffledAcc[0..<maxSize] : []
                shuffledMeta = maxSize > 0 ? shuffledMeta[0..<maxSize] : []
                shuffledStat = maxSize > 0 ? shuffledStat[0..<maxSize] : []
            
                [samples, shuffledAcc, shuffledMeta, shuffledStat]
            } else {
                [samples, [], [], []] // Handle empty groups gracefully
            }
        }
        .flatMap { samples, accessions, meta, sra_stat ->
            def flattened = []
            for (int i = 0; i < accessions.size(); i++) {
                flattened << [samples, accessions[i], meta[i], sra_stat[i]]
            }
            flattened
        }
}