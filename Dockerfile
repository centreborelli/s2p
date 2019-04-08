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
    libtiff5-dev \
    python3 \
    python3-numpy \
    python3-pip

# Install s2p
RUN pip3 install s2p
