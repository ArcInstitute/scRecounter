scRecouter
==========

Nextflow pipeline to re-process all public single-cell RNA-seq data.


# Installation

## Conda & mamba install

`mamba` is needed to run the pipeline. 
It is a faster version of `conda`. `mamba` can be installed via `conda`. 

To install both `conda` and `mamba`, 
see the [conda/mamba Notion docs](https://www.notion.so/arcinstitute/Conda-Mamba-8106bed9553d46cca1af4e10f486bec2).

## Nextflow install

It is easiest to install Nextflow using `mamba`:

```bash
mamba install -n nextflow_env -c bioconda nextflow
```

Make sure to activate the environment before running the pipeline:

```bash
mamba activate nextflow_env
```

## Pipeline install

### Add ssh key to GitHub

> This is only needed if you have not already added your ssh key to GitHub.

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```

* change `your_email@example.com` to your Arc email

```bash
cat ~/.ssh/id_ed25519.pub
```

* copy the output
* GoTo: `GitHub => Settings > SSH and GPG keys > New SSH key`
* Paste the output into the key field
* Add a title (e.g., `Chimera`)
* Click `Add SSH key`

### Clone the repo

```bash
git clone git@github.com:ArcInstitute/scRecounter.git \
  && cd scRecounter
```

### Pipeline conda environments 

The pipeline uses conda environments to manage dependencies. 
Nextflow will automatically create the environments as long as `mamba` is installed.

**Note:** it can take a while to create the environments, even with `mamba`.


# Usage

## Input

### Samples table

Lists the samples and their associated read (fastq) files.

Example:


| sample                     | fastq_1                                                                                                       | fastq_2                                                                                                       |
|-----------------------------|----------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| SRX10188997_SRR13806043     | path/to/reads/SRX10188997_SRR13806043_1.fastq.gz               | path/to/reads/SRX10188997_SRR13806043_2.fastq.gz               |
| SRX10188963_SRR13806077     | path/to/reads/SRX10188963_SRR13806077_1.fastq.gz               | path/to/reads/SRX10188963_SRR13806077_2.fastq.gz               |


### Barcode table

Lists all of the possible barcodes that will be used to determine the cell barcode and UMI for the samples.

Example:

| name               | cell_barcode_length | umi_length | file_path                                                                      |
|--------------------|---------------------|------------|--------------------------------------------------------------------------------|
| 737K-arc-v1        | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-arc-v1.txt             |
| 737K-august-2016   | 16                  | 12         | /large_storage/goodarzilab/public/scRecount/genomes/737K-august-2016.txt        |
| 3M-february-2018   | 16                  | 10         | /large_storage/goodarzilab/public/scRecount/genomes/3M-february-2018.txt        |


## Nextflow run

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --samples
```


# Dev

## Run

```bash
nextflow run main.nf -profile conda,vm
```

## Convert Docker container to Apptainer

Pull the docker image (e.g., `ubuntu:latest`) and convert it to an Apptainer container:

```bash
apptainer pull ubuntu_latest.sif docker://ubuntu:latest
```

## Resources

### Software

* [sra-tools](https://github.com/ncbi/sra-tools)
  * download data from the SRA
* [gget](https://github.com/pachterlab/gget)
  * efficient querying of genomic databases
* [ffq](https://github.com/pachterlab/ffq)
  * Fetch metadata information from the SRA and other databases
  * Can be used to get SRA study accessions from paper DOIs
* [gencube](https://github.com/snu-cdrc/gencube)
  * Efficient retrieval, download, and unification of genomic data from leading biodiversity databases
* [nf-core/fetchngs](https://nf-co.re/fetchngs/1.12.0/)
  * Nextflow pipeline for downloading NGS data

## Workflow

* Input
  * csv of SRA experiment accessions
* Download
  * adapted from [nf-core/fetchngs](https://nf-co.re/fetchngs/1.12.0/)
    * https://github.com/nf-core/fetchngs/blob/master/workflows/sra/main.nf
* Read QC 
  * seqkit stats
* Characterize datasets
  * See https://github.com/ArcInstitute/scRecount/blob/chris/chris_scripts/dsub_solution/process_srr.sh
* Map to reference
  * STARsolo
  * See https://github.com/ArcInstitute/scRecount/blob/chris/chris_scripts/dsub_solution/process_srr.sh