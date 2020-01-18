FROM ubuntu:latest
MAINTAINER Carlo de Franchis <carlodef@gmail.com>
# https://goo.gl/aypXVx
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    cmake \
    gdal-bin \
    geographiclib-tools \
    libfftw3-dev \
    libgdal-dev \
    libgeographic-dev \
    libgeotiff-dev \
    libtiff5-dev \
    python3-dev

# Install pip
RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && python3 get-pip.py

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
RUN pip install -e .

# https://github.com/mapbox/rasterio#ssl-certs
ENV CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
