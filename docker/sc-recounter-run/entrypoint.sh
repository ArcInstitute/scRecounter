#!/bin/bash

# create run name
RUN_NAME="SCRECOUNTER_$(date +"%Y-%m-%d_%H-%M-%S")"

# Set the profile list from the command line arguments
PROFILE_LIST=$(IFS=,; echo "$*")

# Run the pipeline
if [ "$DEPLOYMENT" == "test" ]; then
  WORK_DIR="gs://arc-ctc-nextflow/scRecounter/test/work/${RUN_NAME}"
  OUTPUT_DIR="gs://arc-ctc-screcounter/test/${RUN_NAME}"
  export GCP_SQL_DB_NAME="sragent-test"
else  
  WORK_DIR="gs://arc-ctc-nextflow/scRecounter/prod/work/${RUN_NAME}"
  OUTPUT_DIR="gs://arc-ctc-screcounter/prod4/${RUN_NAME}"
  export GCP_SQL_DB_NAME="sragent-prod"
fi

micromamba run -n sc-recounter-run \
  nextflow run main.nf \
    -profile $PROFILE_LIST \
    -name $RUN_NAME \
    -work-dir $WORK_DIR \
    --output_dir $OUTPUT_DIR \
    -ansi-log false "$@"

# Delete output directory if only nf-report and nf-trace
export GCP_SQL_DB_HOST="35.243.133.29"
export GCP_SQL_DB_USERNAME="postgres"
micromamba run -n sc-recounter-run \
  python cleanup.py $WORK_DIR $OUTPUT_DIR
