# docker for s2p
# Carlo de Franchis <carlodef@gmail.com>

FROM ubuntu:16.04
MAINTAINER Carlo de Franchis <carlodef@gmail.com>
RUN apt-get update && apt-get install -y \
    build-essential \
    gdal-bin \
    geographiclib-tools \
    git \
    libfftw3-dev \
    libgdal-dev \
    libgeographic-dev \
    libgeotiff-dev \
    libtiff5-dev \
    python \
    python-gdal \
    python-numpy \
    python-pip \
    cmake \
    vim
RUN pip install -U pip
RUN pip install utm bs4 lxml requests

# install s2p
RUN git clone https://github.com/MISS3D/s2p.git --branch master --single-branch --depth 1 --recursive
RUN cd s2p && make

WORKDIR /s2p
