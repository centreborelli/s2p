# S2P - Satellite Stereo Pipeline

[![Build Status](https://travis-ci.com/cmla/s2p.svg?branch=master)](https://travis-ci.com/cmla/s2p)
[![PyPI version](https://img.shields.io/pypi/v/s2p)](https://pypi.org/project/s2p)

S2P is a Python library and command line tool that implements a stereo
pipeline which produces elevation models from images taken by high resolution
optical satellites such as Pléiades, WorldView, QuickBird, Spot or Ikonos. It
generates 3D point clouds and digital surface models from stereo pairs (two
images) or tri-stereo sets (three images) in a completely automatic fashion.

S2P was used to win the 2016 [IARPA Multi-View Stereo 3D Mapping Challenge](https://www.iarpa.gov/challenges/3dchallenge.html).

A wide variety of stereo correlation algorithms are supported, including several
flavors of semi-global matching (SGM), TV-L1 optical flow, etc.

The main language is Python, although several operations are handled by
binaries written in C.

The pipeline is implemented in the Python package `s2p`. It can be used
to produce surface models and 3D point clouds from arbitrarily large regions
of interest or from complete images. If needed, it cuts the region of interest
in several small tiles and process them in parallel.

Its main source code repository is https://github.com/cmla/s2p.


# Dependencies

## GDAL
The main dependency is GDAL. Version 2.1.0 or newer is required.

### On Ubuntu
`gdal` can be installed with `apt-get`. In order to get a recent version we
recommend adding the PPA `ubuntugis-unstable` (first command below):

    sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable
    sudo apt-get update
    sudo apt-get install libgdal-dev gdal-bin

### On macOS
[Download](http://www.kyngchaos.com/files/software/frameworks/GDAL_Complete-2.4.dmg)
and install the `.dmg` file. Update your `PATH` after the installation by
running this command:

    export PATH="$PATH:/Library/Frameworks/GDAL.framework/Programs"

Copy it in your `~/.profile`.

## Other dependencies (geographiclib, fftw)

On Ubuntu:

    apt-get install build-essential geographiclib-tools libgeographic-dev libfftw3-dev libgeotiff-dev libtiff5-dev

On macOS:

    brew install geographiclib fftw


# Installation

    pip install s2p

Alternatively, if you want to get the latest commit or want to edit the
sources, install it in editable mode from a git clone:

    git clone https://github.com/cmla/s2p.git --recursive
    cd s2p
    pip install -e .

The `--recursive` option for `git clone` allows to clone all git submodules, such
as the [iio](https://github.com/mnhrdt/iio) library.

If the `--recursive` option wasn't used when cloning, the submodules can now be
retrieved with

    git submodule update --init

All `s2p` python submodules are located in the `s2p` package. Some python
functions of these modules rely on external binaries. Most of these binaries
were written on purpose for the needs of the pipeline, and their source code is
provided here in the `c` folder. For the other binaries, the source code is
provided in the `3rdparty` folder.

All the sources (ours and 3rdparties) are compiled from the same makefile. Just
run `make all` from the `s2p` folder to compile them. This will create a `bin`
directory containing all the needed binaries. This makefile is used when
running `pip install .`

You can test if S2P is correctly working using:

    make test

## Docker image
[![Docker Status](http://dockeri.co/image/cmla/s2p)](https://hub.docker.com/r/cmla/s2p/)

A precompiled docker image is available and ready to use:

    docker pull cmla/s2p


# Usage

`s2p` is a Python library that can be imported into other applications. It also
comes with a Command Line Interface (CLI).

## From the command line

The `s2p` CLI usage instructions can be printed with the `-h` and `--help` switches.

    $ s2p -h
    usage: s2p.py [-h] config.json

    S2P: Satellite Stereo Pipeline

    positional arguments:
      config.json           path to a json file containing the paths to input and
                            output files and the algorithm parameters

    optional arguments:
      -h, --help            show this help message and exit

To run the whole pipeline, call `s2p` with a json configuration file as unique argument:

    s2p tests/data/input_pair/config.json

All the parameters of the algorithm, paths to input and output data are stored
in the json file. See the provided `test.json` file for an example, and the
comments in the file `s2p/config.py` for some explanations about the roles
of these parameters.

Notice that each input image must have RPC coefficients, either in its GeoTIFF
tags or in a companion `.xml` or `.txt` file.

#### ROI definition

The processed Region of interest (ROI) is defined by the image coordinates (x,
y) of its top-left corner, and its dimensions (w, h) in pixels. These four
numbers must be given in the `json` configuration file, as in the `test.json`
example file. They are ignored if the parameter `'full_img'` is set to `true`.
In that case the full image will be processed.

#### File paths in json configuration files

In the json configuration files, input and output paths are relative to the json
file location, not to the current working directory.


### MicMac (optional)

If you want to use MicMac for the stereo matching step, you must install it
first and create a symlink to the micmac directory (the one containing a 'bin'
folder with a bunch of executables in it, among with 'MICMAC' and 'mm3d') in
the 'bin' folder:

    ln -s PATH_TO_YOUR_MICMAC_DIR bin/micmac


## References

If you use this software please cite the following papers:

[*An automatic and modular stereo pipeline for pushbroom
images*](http://dx.doi.org/10.5194/isprsannals-II-3-49-2014), Carlo de
Franchis, Enric Meinhardt-Llopis, Julien Michel, Jean-Michel Morel, Gabriele
Facciolo. ISPRS Annals 2014.

[*On Stereo-Rectification of Pushbroom
Images*](http://dx.doi.org/10.1109/ICIP.2014.7026102), Carlo de Franchis, Enric
Meinhardt-Llopis, Julien Michel, Jean-Michel Morel, Gabriele Facciolo.  ICIP
2014.

[*Automatic sensor orientation refinement of Pléiades stereo
images*](http://dx.doi.org/10.1109/IGARSS.2014.6946762), Carlo de Franchis,
Enric Meinhardt-Llopis, Julien Michel, Jean-Michel Morel, Gabriele Facciolo.
IGARSS 2014.
