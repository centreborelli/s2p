# docker for s2p
# Carlo de Franchis <carlodef@gmail.com>

FROM ubuntu:16.04
MAINTAINER Carlo de Franchis <carlodef@gmail.com>
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gdal-bin \
    git \
    libfftw3-dev \
    libgdal-dev \
    libgeotiff-dev \
    libtiff5-dev \
    libtiff-tools \
    libxslt1-dev \
    python \
    python-gdal \
    python-numpy \
    python-pip
RUN pip install -U pip
RUN pip install utm

# install s2p
RUN git clone https://github.com/carlodef/s2p.git --branch master --single-branch --depth 1 s2p
RUN cd s2p && make

WORKDIR /s2p
