# docker for s2p
# Carlo de Franchis <carlodef@gmail.com>

FROM ubuntu:16.04
MAINTAINER Carlo de Franchis <carlodef@gmail.com>
# https://goo.gl/aypXVx
ARG DEBIAN_FRONTEND=noninteractive 
RUN apt-get update && apt-get install -y \
    build-essential \
    geographiclib-tools \
    git \
    libfftw3-dev \
    libgdal-dev \
    libgeographic-dev \
    libgeotiff-dev \
    libtiff5-dev \
    python \
    python-numpy \
    python-pip \
    cmake \
    software-properties-common \
    python-software-properties \
    unzip
RUN pip install -U pip
RUN pip install utm bs4 lxml requests

# Install GDAL 2.x from ubuntugis-stable
RUN apt-get install dialog apt-utils -y
RUN add-apt-repository ppa:ubuntugis/ppa -y 
RUN apt-get update && apt-get install -y \
    gdal-bin \
    python-gdal

# Install s2p from MISS3D/s2p
RUN git clone https://github.com/MISS3D/s2p.git --recursive
RUN cd s2p && make all

WORKDIR /s2p
