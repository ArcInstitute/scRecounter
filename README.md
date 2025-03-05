scRecouter
==========

A Nextflow pipeline to re-process single-cell RNA-seq data from the Sequence Read Archive.

# Workflow

* **User provides:**
  * A table of samples & associated accessions
    * Alternatively, the pipeline can pull accessions from the scRecounter SQL database
  * Associated files required:
    * A table of barcodes to use for cell barcode and UMI identification
    * A table of STAR index directories to use for mapping
* **Pipeline:**
  * Load accessions from provided table or SQL database
  * For each accession:
    * Use `fastq-dump` to download a subset of reads as fastq files from the SRA
    * Determine the "best" STAR parameters by mapping the reads using various parameter combinations
      * Parameters: version of cell barcodes, cell barcode length, UMI length, strand, STAR reference index
      * The STAR parameters are selected based on the fraction of valid barcodes
    * Download all reads with `fasterq-dump`
      * If download fails, try again with `fastq-dump` using a max of `fallback_max_spots` reads (see `nextflow.config`).
    * Map the reads with STARsolo using the "best" STAR parameters

# Manuscript

**scBaseCamp: An AI agent-curated, uniformly processed, and continually expanding single cell data repository**.
Nicholas D Youngblut, Christopher Carpenter, Jaanak Prashar, Chiara Ricci-Tam, Rajesh Ilango, Noam Teyssier,
Silvana Konermann, Patrick Hsu, Alexander Dobin, David P Burke, Hani Goodarzi, Yusuf H Roohani.
bioRxiv 2025.02.27.640494; doi: [https://doi.org/10.1101/2025.02.27.640494](https://doi.org/10.1101/2025.02.27.640494)

# Installation

## Conda & mamba install

`mamba` is needed to run the pipeline. 
It is a faster version of `conda`. 
`mamba` can be installed via `conda`. 
You can use `conda` instead of `mamba` if you prefer.

## Nextflow install

It is easiest to install Nextflow using `mamba` (or `conda`).

```bash
mamba create -n nextflow_env -c bioconda nextflow
```

Make sure to activate the environment before running the pipeline:

```bash
mamba activate nextflow_env
```

All other dependencies will be installed by Nextflow.


## Pipeline install

### Clone the repo

```bash
git clone https://github.com/ArcInstitute/scRecounter.git \
  && cd scRecounter
```

### Pipeline conda environments (if running locally)

The pipeline uses conda environments to manage dependencies. 
Nextflow will automatically create the environments as long as `mamba` is installed.

**Note:** it can take a while to create the environments, even with `mamba`.

### Pipeline Docker containers (if running on GCP) 

The pipeline defaults to using custom Docker containers hosted on Google Artifact Registry.

You can build the Docker containers yourself. See [./docker/README.md](./docker/README.md) for details.
Be sure to update the [profiles.config](./config/profiles.config) file to point to the new containers.

# Usage

## Input

### Accessions table

Lists the samples and their associated SRA experiment accessions.

> This table is not required if the pipeline is pulling accessions from the scRecounter SQL database.
  To pull accessions from the database, do not provide `--accessions` via the command line.

Example:

| sample      | accession   | organism |
|-------------|-------------|----------|
| SRX22716300 | SRR27024456 | human    |
| SRX25994842 | SRR30571763 | mouse    |

> `organism` is optional. It will determine the STAR index to use for mapping. Otherwise all indexes will be used for parameter selection.

### Barcode table

Lists all of the possible barcodes that will be used to determine the cell barcode and UMI for the samples.

Example:

| name             | cell_barcode_length | umi_length | file_path                                                                |
|------------------|---------------------|------------|--------------------------------------------------------------------------|
| 737K-arc-v1      | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-arc-v1.txt      |
| 737K-august-2016 | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-august-2016.txt |
| 3M-february-2018 | 16                  | 10         | /large_storage/goodarzilab/public/scRecount/genomes/3M-february-2018.txt |


### STAR index table

Lists the STAR index files that will be used to map the reads.

Example:

| Organism | Star Index Path                                                                   |
|----------|-----------------------------------------------------------------------------------|
| human    | /large_storage/goodarzilab/public/scRecount/genomes/star_refData_2020_hg38        |
| mouse    | /large_storage/goodarzilab/public/scRecount/genomes/star2.7.11_refData_2020_mm10  |


> If `organism` is provided in the `Accessions` table, the STAR index will be selected based on the `organism` column.
  Thus, it reduces the number of parameter combinations that need to be tested.

## Nextflow run 

### Test runs

Local run with provided accessions:

```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev
```

Local run with provided accessions (problematic datasets)

```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev_problems
```

With conda, accessions pulled from scRecounter database:

```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,no_acc_dev
```

GCP run with provided accessions:

```bash
nextflow run main.nf \
  -profile docker,trace,report,gcp,gcp_dev,dev,acc_dev
```

GCP run with accessions pulled from scRecounter SQL database:

```bash
nextflow run main.nf \
  -profile docker,trace,report,gcp,gcp_dev,dev,no_acc_dev
```

### Characterize datasets

Use just a small subset of reads in the dataset to identify library prep method, species, etc.

```bash
nextflow run /home/nickyoungblut/dev/nextflow/scRecounter/main.nf \
  -work-dir gs://arc-ctc-nextflow/scRecounter/work \
  -profile docker,gcp \
  -ansi-log false \
  --max_spots 100000 \
  --output_dir gs://arc-ctc-nextflow/scRecounter/results/ \
  --accessions TMP/SRX22716300.csv
```

### Deploy to GCP Cloud Run

See [./docker/sc-recounter-run/README.md](./docker/sc-recounter-run/README.md) for details.


# Contributing

Feel free to fork the repository and submit a pull request.
However, the top priority is to keep SRAgent functioning 
for the ongoing scBaseCamp project.