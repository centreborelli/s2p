FROM ubuntu:latest
MAINTAINER Carlo de Franchis <carlodef@gmail.com>
# https://goo.gl/aypXVx
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gdal-bin \
    geographiclib-tools \
    libfftw3-dev \
    libgdal-dev \
    libgeographic-dev \
    libgeotiff-dev \
    libgsl-dev \
    libtiff5-dev \
    python3 \
    python3-numpy \
    python3-pip

# Copy files needed to install s2p
WORKDIR /root
COPY LICENSE.txt MANIFEST.in README.md setup.py makefile s2p/
WORKDIR /root/s2p
COPY 3rdparty/ 3rdparty/
COPY c/ c/
COPY s2p/ s2p/

# Install s2p
RUN pip3 install -e .
