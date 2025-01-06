## sc-recounter-run

Env vars

```bash
IMG_NAME=sc-recounter-run
IMG_VERSION=0.1.2
REGION="us-east1"
GCP_PROJECT_ID="c-tc-429521"
SERVICE_ACCOUNT_EMAIL="nick-nextflow@c-tc-429521.iam.gserviceaccount.com"
SERVICE_ACCOUNT_JSON="c-tc-429521-6f6f5b8ccd93.json"
```

### Docker

Build

> from the base directory of the repository

```bash
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --build-arg CONDA_ENV_YAML=docker/${IMG_NAME}/environment.yml \
  --platform linux/amd64 \
  --tag ${IMG_NAME}:${IMG_VERSION} \
  .
```

Run the image (`-help`)

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  -v ${HOME}/.gcp/:/.gcp \
  --env GOOGLE_APPLICATION_CREDENTIALS=/.gcp/${SERVICE_ACCOUNT_JSON} \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION} \
  -help
```

Run the image (`-profile`)

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  -v ${HOME}/.gcp/:/.gcp \
  --env GOOGLE_APPLICATION_CREDENTIALS=/.gcp/${SERVICE_ACCOUNT_JSON} \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION} \
  -profile docker,gcp,gcp_dev,dev,no_acc_dev
```

Run with bash entrypoint

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  -v ${HOME}/.gcp/:/.gcp \
  --env GOOGLE_APPLICATION_CREDENTIALS=/.gcp/${SERVICE_ACCOUNT_JSON} \
  --entrypoint /bin/bash \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION}
```

### GCP Artifact Registry

Create (if needed)

```bash
DESCRIPTION="Run the scRecounter nextflow pipeline"
gcloud artifacts repositories create ${IMG_NAME} \
  --repository-format=docker \
  --project=${GCP_PROJECT_ID} \
  --location=${REGION} \
  --description="${DESCRIPTION}" \
  --async
```

Push

```bash
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${GCP_PROJECT_ID}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${GCP_PROJECT_ID}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```

### GCP Cloud Run Jobs

Create/update the job

```bash
JOB_NAME="${IMG_NAME}"
gcloud beta run jobs update ${JOB_NAME} \
  --service-account=${SERVICE_ACCOUNT_EMAIL} \
  --project=${GCP_PROJECT_ID} \
  --region=${REGION} \
  --image=${REGION}-docker.pkg.dev/${GCP_PROJECT_ID}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  --set-env-vars=TZ=America/Los_Angeles \
  --cpu=2 \
  --memory=6Gi \
  --task-timeout=4320m \
  --args=""
```

