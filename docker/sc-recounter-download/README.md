sc-recounter-download container
===============================

# Build and push to GCP Container Registry

Env vars

```bash
IMG_NAME=sc-recounter-download
IMG_VERSION=0.1.0
REGION="us-east1"
PROJECT="c-tc-429521"
```

Build

> from the base directory of the repository

```bash
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --build-arg CONDA_ENV_YAML=envs/download.yml \
  --platform linux/amd64 \
  --tag ${IMG_NAME}:${IMG_VERSION} \
  .
```

Run

```bash
docker run -it --rm \
  -u $(id -u):$(id -g) \
  -v ${PWD}:/data \
  -v ${HOME}/.gcp/:/.gcp \
  --platform linux/amd64 \
  ${IMG_NAME}:${IMG_VERSION}
```

Push

```bash
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```
