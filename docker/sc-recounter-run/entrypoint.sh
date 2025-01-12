#!/bin/bash

# create run name
RUN_NAME="SCRECOUNTER_$(date +"%Y-%m-%d_%H-%M-%S")"

# Activate the micromamba environment and run the pipeline
micromamba run -n sc-recounter-run \
  nextflow run main.nf \
    -profile docker,gcp,report,trace \
    -name $RUN_NAME \
    -work-dir "gs://arc-ctc-nextflow/scRecounter/prod/work/${RUN_NAME}" \
    --output_dir "gs://arc-ctc-screcounter/prod2/${RUN_NAME}" \
    -ansi-log false "$@"

# Delete output directory if only nf-report and nf-trace
export GCP_SQL_DB_HOST="35.243.133.29"
export GCP_SQL_DB_NAME="sragent-prod"
export GCP_SQL_DB_USERNAME="postgres"
micromamba run -n sc-recounter-run \
  python cleanup.py \
    "gs://arc-ctc-nextflow/scRecounter/prod/work/${RUN_NAME}" \
    "gs://arc-ctc-screcounter/prod2/${RUN_NAME}"
