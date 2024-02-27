ARG UBUNTU_VERSION=22.04
ARG BASE_IMAGE=ubuntu:${UBUNTU_VERSION}

# Base image
FROM ${BASE_IMAGE}

# Install system dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    ca-certificates \
    ccache \
    cmake \
    curl \
    git \
    libjpeg-dev \
    libpng-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy in plink2 binary from local directory
COPY plink2 /usr/local/bin/plink2

# Copy in PRSice2 linux binary from local directory
COPY PRSice_linux /usr/local/bin/PRSice_linux
RUN chmod +x /usr/local/bin/PRSice_linux