FROM ubuntu:20.04

LABEL maintainer="Carlo de Franchis <carlodef@gmail.com>"

# https://goo.gl/aypXVx
ARG DEBIAN_FRONTEND=noninteractive
# https://github.com/mapbox/rasterio#ssl-certs
ENV CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    cmake \
    libfftw3-dev \
    libgdal-dev \
    libgeotiff-dev \
    libtiff5-dev \
    python3-dev \
    python3-pip && \
    rm -fr /var/lib/apt/lists/*

# Copy files needed to install s2p
WORKDIR /root
COPY LICENSE.txt MANIFEST.in README.md setup.py makefile s2p/
WORKDIR /root/s2p
COPY 3rdparty/ 3rdparty/
COPY c/ c/
COPY s2p/ s2p/
COPY bin/ bin/
COPY lib/ lib/

# Install s2p
RUN pip install --no-cache-dir -e .

