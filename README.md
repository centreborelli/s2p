# S2P - Satellite Stereo Pipeline

[![Build Status](https://travis-ci.org/MISS3D/s2p.svg?branch=master)](https://travis-ci.org/MISS3D/s2p)

This software implements a stereo pipeline which produces elevation models from
images taken by high resolution optical satellites such as Pléiades, WorldView,
QuickBird, Spot or Ikonos. It generates 3D point clouds and digital surface
models from stereo pairs (two images) or tri-stereo sets (three images) in a
completely automatic fashion.

S2P was used to win the 2016 [IARPA Multi-View Stereo 3D Mapping Challenge](https://www.iarpa.gov/challenges/3dchallenge.html).

A wide variety of stereo correlation algorithms are supported, including several
flavors of semi-global matching (SGM), TV-L1 optical flow, etc.

The main language is Python, although several operations are handled by
binaries written in C.

The pipeline is implemented in the file `s2p.py`. The `s2p` module can be used
to produce surface models and 3D point clouds from arbitrarily large regions
of interest or from complete images. If needed, it cuts the region of interest
in several small tiles and process them in parallel.


## Docker image
[![Docker Status](http://dockeri.co/image/carlodef/s2p)](https://hub.docker.com/r/carlodef/s2p/)

A precompiled docker image is available and ready to use:

    docker pull carlodef/s2p


## Usage

Run the python script `s2p.py` with a json configuration file as unique
argument:

    python s2p.py testdata/input_pair/config.json

All the parameters of the algorithm, paths to input and output data are stored
in the json file. See the provided `test.json` file for an example, and the
comments in the file `s2plib/config.py` for some explanations about the roles
of these parameters.

Notice that each input image must be accompanied by an `xml` file containing rpc
coefficients.

#### ROI definition

The processed Region of interest (ROI) is defined by the image coordinates (x,
y) of its top-left corner, and its dimensions (w, h) in pixels. These four
numbers must be given in the `json` configuration file, as in the `test.json`
example file. They are ignored if the parameter `'full_img'` is set to `true`.
In that case the full image will be processed.

If neither the ROI definition or the `'full_img'` flag are present in the
configuration file, then a preview of the reference image must be provided. The
ROI will be selected interactively on that preview. The path of the preview
file must be given by the key `'prv'` of the `'images'[0]` dictionary (as in
the example).

#### File paths in json configuration files

In the json configuration files, input and output paths are relative to the json
file location, not to the current working directory.

## Dependencies

Required dependencies (Ubuntu 16.04):

    add-apt-repository -y ppa:ubuntugis/ppa  # The repository is added so that the version >= 2.1 of gdal is installed (requirement)
    apt-get update
    apt-get install build-essential cmake gdal-bin geographiclib-tools libgeographic-dev libfftw3-dev libgdal-dev libgeotiff-dev libtiff5-dev python python-gdal python-numpy python-pip

and

    pip install utm bs4 lxml requests

`gdal` version must be 2.1.0 or newer.


## Installation

    git clone https://github.com/MISS3D/s2p.git --recursive
    cd s2p
    make all

The `--recursive` option for `git clone` allows to clone all submodules, such
as the [iio](https://github.com/mnhrdt/iio) library.

If the `--recursive` option wasn't used when cloning, the submodules can now be
retrieved with

    git submodule update --init

All the python modules are located in the `s2plib` folder. Some python
functions of these modules rely on external binaries. Most of these binaries
were written on purpose for the needs of the pipeline, and their source code is
provided here in the `c` folder. For the other binaries, the source code is
provided in the `3rdparty` folder.

All the sources (ours and 3rdparties) are compiled from the same makefile. Just
run `make all` from the `s2p` folder to compile them.  This will create a `bin`
directory containing all the needed binaries.

You can test if S2P is correctly working using:

    make test

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
