# docker run -it --rm gcr.io/ucsd-medicine-cast/nonlin_prs_prs_aml:latest /bin/bash

# Definitions
project_id = ucsd-medicine-cast
identifier = nonlin_prs_prs_aml
version = 0.0.1

# GCR prefix
gcr_prefix = gcr.io/$(project_id)

# Image name and tag
name = $(gcr_prefix)/$(identifier)
tag = $(version)

# Build
build:
	# cp ../../resources/plink2 .
	docker build \
		--progress=plain \
		--platform linux/amd64 \
		--build-arg CACHEBUST=$(date +%s) \
		-t $(name):$(tag) .
	docker tag $(name):$(tag) $(name):latest

# Push
push:
	docker push $(name):$(tag)
	docker push $(name):latest