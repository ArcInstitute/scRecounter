gcp-loader
==========

A simple Nextflow pipeline for efficiently loading single-cell data as h5ad files onto GCP


# Dev

Local run

```bash
nextflow run main.nf -profile conda,vm,dev --feature_type GeneFull -resume
```

Slurm run

```bash
nextflow run main.nf -profile conda,slurm,dev -resume 
```

## prod

### GeneFull_Ex50pAS

#### CZI

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --feature_type GeneFull_Ex50pAS \
  --organisms "Mus musculus,Homo sapiens,Macaca mulatta" \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

#### SRA

```bash
nextflow run main.nf -profile conda,slurm --feature_type GeneFull_Ex50pAS
```

## Velocyto

### CZI

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --feature_type Velocyto \
  --organisms "Mus musculus,Homo sapiens,Macaca mulatta" \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

### SRA

```bash
nextflow run main.nf -profile conda,slurm --feature_type Velocyto
```


### Gene 

### CZI

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --feature_type Gene \
  --organisms "Mus musculus,Homo sapiens,Macaca mulatta" \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

### SRA

```bash
nextflow run main.nf -profile conda,slurm --feature_type Gene
```

### GeneFull

### CZI

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --feature_type GeneFull \
  --organisms "Mus musculus,Homo sapiens,Macaca mulatta" \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

### SRA

```bash
nextflow run main.nf -profile conda,slurm --feature_type GeneFull
```


### GeneFull_ExonOverIntron

### CZI

```bash
nextflow run main.nf \
  -profile conda,slurm \
  --feature_type GeneFull_ExonOverIntron \
  --organisms "Mus musculus,Homo sapiens,Macaca mulatta" \
  --input_dir /processed_datasets/scRecount/cellxgene/counted_SRXs
```

### SRA

```bash
nextflow run main.nf -profile conda,slurm --feature_type GeneFull_ExonOverIntron
```

