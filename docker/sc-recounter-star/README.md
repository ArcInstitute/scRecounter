sc-recounter-star container
===========================

# Build and push to GCP Container Registry

Env vars

```bash
IMG_NAME=sc-recounter-star
IMG_VERSION=0.1.0
REGION="us-east1"
PROJECT="c-tc-429521"
```

Build

> from the base directory of the repository

```bash
docker build \
  --file docker/${IMG_NAME}/Dockerfile \
  --build-arg CONDA_ENV_YAML=envs/star.yml \
  --platform linux/amd64 \
  --tag ${IMG_NAME}:${IMG_VERSION} \
  .
```

Push

```bash
docker tag ${IMG_NAME}:${IMG_VERSION} \
  ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION} \
  && docker push ${REGION}-docker.pkg.dev/${PROJECT}/${IMG_NAME}/${IMG_NAME}:${IMG_VERSION}
```
