scRecouter
==========

Nextflow pipeline to re-process public single-cell RNA-seq data.


# Workflow

* User provides:
  * A table of samples & associated accessions
    * Alternatively, the pipeline can pull accessions from the scRecounter database
  * Associated files required:
    * A table of barcodes to use for cell barcode and UMI identification
    * A table of STAR index directories to use for mapping
* Pipeline:
  * Read file formatting and QC
    * includes assessing read length
  * Cell barcode, UMI, and strand identification
    * via STARsolo with a subset of reads
      * assess fraction of valid barcodes
  * The STAR parameter setting with the highest fraction of valid barcodes is used for the full STAR run
  * STAR run with selected parameters
    * All final count tables are "published" to the output directory

# Installation

## Conda & mamba install

`mamba` is needed to run the pipeline. 
It is a faster version of `conda`. `mamba` can be installed via `conda`. 

To install both `conda` and `mamba`, 
see the [conda/mamba Notion docs](https://www.notion.so/arcinstitute/Conda-Mamba-8106bed9553d46cca1af4e10f486bec2).

## Nextflow install

It is easiest to install Nextflow using `mamba`:

```bash
mamba create -n nextflow_env -c bioconda nextflow
```

Make sure to activate the environment before running the pipeline:

```bash
mamba activate nextflow_env
```

## Pipeline install

### Clone the repo

You will need a GitHub PAT to clone the repo.

To obtain a GitHub PAT:
 - Go to your GitHub account settings
 - Click on Developer settings
 - Click on Personal access tokens
 - Click on Generate new token

```bash
git clone https://github.com/ArcInstitute/scRecounter.git \
  && cd scRecounter
```

> Provide the GitHub PAT when prompted for a password

To cache your GitHub credentials:

```bash
git config --global credential.helper cache
```

### Pipeline conda environments (if running locally)

The pipeline uses conda environments to manage dependencies. 
Nextflow will automatically create the environments as long as `mamba` is installed.

**Note:** it can take a while to create the environments, even with `mamba`.

### Pipeline Docker containers (if running on GCP) 

The containers are hosted on the Google Artifact Registry.

Be sure to set the `GOOGLE_APPLICATION_CREDENTIALS` environment variable to a service
account with the required permissions.


# Usage

## Input

### Accessions table

Lists the samples and their associated SRA experiment accessions.

> This table is not required if the pipeline is pulling accessions from the scRecounter database.
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

## Nextflow run 

### Test runs

Local run with provided accessions:

```bash
nextflow run main.nf \
  -work-dir tmp/work \
  -profile conda,trace,report,vm,vm_dev,dev,acc_dev
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

GCP run with accessions pulled from scRecounter database:

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

***

# Dev

## deploy new GCP Cloud Run batch revision

See [./docker/sc-recounter/README.md](./docker/sc-recounter/README.md) for details.

## set env variables

```bash
export GOOGLE_APPLICATION_CREDENTIALS="${HOME}/.gcp/c-tc-429521-6f6f5b8ccd93.json"
```

## Run locally

Basic run:

```bash
nextflow run main.nf \
  -profile conda,vm,dev,vm_dev,acc_dev \
  -resume 
```

Run with problematic accessions

```bash
nextflow run main.nf \
  -profile conda,vm  \
  --define \
  --max_spots 50000 \
  --accessions data/accessions_problems.csv \
  --output_dir tmp/accessions_problems \
  -resume
```


## SLURM, with defined accessions

A couple of small-data accessions

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --accessions data/accessions_small_n2.csv \
  --output_dir tmp/results_small_n2
```

Many small-data accessions, subsampled

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --max_spots 500000 \
  --accessions data/accessions_small_n10.csv \
  --output_dir tmp/results_small_n10
```

## Convert accessions

Use the `acc2srr.py` script in this repo. An example:

```bash
./scripts/acc2srr.py --email YOUR_EMAIL@arcinstitute.org accessions.txt 
```

The `accessions.txt` file should contain a list of GEO or SRA accessions (one per line). 

## Convert Docker container to Apptainer

Pull the docker image (e.g., `ubuntu:latest`) and convert it to an Apptainer container:

```bash
apptainer pull ubuntu_latest.sif docker://ubuntu:latest
```

## GCP VM setup

Create VM

```bash
gcloud compute instances create sc-recounter-vm \
    --project=c-tc-429521 \
    --zone=us-east1-c \
    --machine-type=e2-standard-8 \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --boot-disk-size=50GB \
    --boot-disk-type=pd-balanced \
    --tags=allow-http,allow-https \
    --scopes=storage-full
```

ssh onto the VM

```bash
gcloud compute ssh sc-recounter-vm \
  --zone=us-east1-c \
  --project=c-tc-429521 \
  --impersonate-service-account=${HOME}/.gcp/c-tc-429521-6f6f5b8ccd93.json
```

Install micromamba

```bash
curl -L \
  -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
  && bash Miniforge3-$(uname)-$(uname -m).sh
```

Add the required channels

```bash
conda config --add channels nodefaults \
  && conda config --add channels pytorch \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge
```

Install nextflow

```bash
conda install -n nextflow-env nextflow \
  && conda activate nextflow-env
```

> Be sure that the `GOOGLE_APPLICATION_CREDENTIALS` env variable is set to the path of the service account key

The service account needs the following roles:

```bash
# batch
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/batch.serviceAgent"

# compute
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/compute.admin"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/compute.instanceAdmin.v1"

# storage
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/storage.objectViewer"

# Artifact Registry
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/artifactregistry.reader"

# network and Monitoring Access
gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/iam.serviceAccountUser"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/compute.networkUser"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/logging.viewer"

gcloud projects add-iam-policy-binding $GCP_PROJECT \
    --member="serviceAccount:$SERVICE_ACCOUNT" \
    --role="roles/monitoring.viewer"
```

Run the pipeline

```bash
nextflow run main.nf -profile docker,gcp,acc_dev
```

To stop the VM:

```bash
gcloud compute instances stop sc-recounter-vm --zone=us-east1-b
```

### Deploy to GCP Cloud Run

See [./docker/sc-recounter-run/README.md](./docker/sc-recounter-run/README.md) for details.

## Resources

### Software

* [sra-tools](https://github.com/ncbi/sra-tools)
  * download data from the SRA
* [gget](https://github.com/pachterlab/gget)
  * efficient querying of genomic databases
* [ffq](https://github.com/pachterlab/ffq)
  * Fetch metadata information from the SRA and other databases
  * Can be used to get SRA study accessions from paper DOIs
* [pysradb](https://github.com/saketkc/pysradb)
  * A Python package for retrieving metadata from SRA/ENA/GEO
* [gencube](https://github.com/snu-cdrc/gencube)
  * Efficient retrieval, download, and unification of genomic data from leading biodiversity databases
* [geofetch](https://pep.databio.org/geofetch/)
  * Downloads and organizes data and metadata from GEO and SRA
* [nf-core/fetchngs](https://nf-co.re/fetchngs/1.12.0/)
  * Nextflow pipeline for downloading NGS data

### Databases

* [Single cell studies database](https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/edit?gid=0#gid=0)

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


