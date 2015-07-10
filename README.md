# S2P - Satellite Stereo Pipeline

This code implements a stereo pipeline which produces elevation models from
images taken by high resolution optical satellites such as Pléiades, WorldView,
QuickBird, Spot or Ikonos. It generates 3D point clouds and digital surface
models from stereo pairs (two images) or tri-stereo sets (three images) in a
complete automatic fashion.

The main language is Python, although several operations are handled by
binaries written in C.

The pipeline is implemented in the file `s2p.py`. The `s2p` module can be used
to produce surface models and 3D point clouds from arbitrarily large regions
of interest or from complete images. If needed, it cuts the region of interest
in several small tiles and process them in parallel.

## Usage

### From a shell
The easiest way to use `s2p` is to run the binary `s2p.py` from a shell, with a
json configuration file as unique argument:

    $ ./s2p.py test.json

All the parameters of the algorithm, paths to input and output data are stored
in the json file. See the provided `test.json` file for an example, and the
comments in the file `python/config.py` for some explanations about the roles
of these parameters.

#### ROI definition

The processed Region of interest (ROI) is defined by the image coordinates (x,
y) of its top-left corner, and its dimensions (w, h) in pixels. These four
numbers must be given in the `json` configuration file (as in the example).
They are ignored if the parameter `'full_img'` is set to `true`. In that case
the full image will be processed. If neither the ROI definition or the
`'full_img'` flag are present in the configuration file, then a preview of the
reference image must be provided. The ROI will be selected interactively on
that preview. The path of the preview file must be given by the key `'prv'` of
the `'images'[0]` dictionary (as in the example).

### From a python interpreter
Another way is to import the `s2p` module in a python session, and run the
functions `process_pair` or `process_triplet`, depending on the kind of dataset
you have (stereo pair or triplet), to generate an altitude map.  A 3D point
cloud can be generated from this altitude map using the function
`generate_cloud`, and then a raster digital surface model map can be generated
from the 3D point cloud.

    python
    >>> import s2p
    >>> alt_map = s2p.process_triplet('test', 'pan1.tif', 'rpc1.xml', 'pan2.tif', 'rpc2.xml', 'pan3.tif', 'rpc3.xml', 25150, 24250, 300, 300, 3)
    >>> s2p.generate_cloud('test', 'pan1.tif', 'rpc1.xml', 'xs1.tif', 25150, 24250, 300, 300, alt_map)
    >>> s2p.generate_dsm('test/dsm.tif', ['test/cloud.ply'], 4)

See the docstrings of the functions `process_pair`, `process_triplet`,
`generate_cloud` and `generate_dsm` for a complete description of their
arguments.

## Installation

All the python modules are located in the `python` folder. Some python
functions of these modules rely on external binaries. Most of these binaries
were written on purpose for the needs of the pipeline, and their source code is
provided here in the `c` folder. For the other binaries, the source code is
provided in the `3rdparty` folder.

All the sources (ours and 3rdparties) are compiled from the same makefile. Just
run `make` from the `s2p` folder to compile them.  This will create a `bin`
directory containing all the needed binaries.

## Dependencies

Required dependencies: `cmake, libtiff, libpng, libjpeg, libfftw3, libgeotiff,
gdal`. These can be installed on Debian 7 (wheezy) with

    apt-get install cmake libfftw3-dev libtiff5-dev libgeotiff-dev gdal-bin

`gdal` version must be 1.10 or newer. Comments about dependencies that were
required in previous versions of `s2p` are still here just in case.

### GDAL >= 1.10

    cd 3rdparty
    wget http://download.osgeo.org/gdal/1.10.1/gdal-1.10.1.tar.gz
    tar xzf gdal-1.10.1.tar.gz
    cd gdal-1.10.1
    ./configure --prefix=$HOME/local
    make install

### Old dependencies. Not required anymore

Here are the instuctions to install them in a local folder, without root
privileges. First create a local folder:

    mkdir ~/local

If you want to use a different path, do the appropriate changes in the
following instructions.

#### GeographicLib

    cd 3rdparty/GeographicLib-1.32
    mkdir build
    cd build
    cmake -D CMAKE_INSTALL_PREFIX=~/local -D GEOGRAPHICLIB_DATA=~/local/share/GeographicLib ..
    make
    make install

Since we use GeographicLib to evaluate geoid heights we must also install the
geoids data files by running the following script, which has been configured by
cmake:

    ~/local/sbin/geographiclib-get-geoids


#### SGBM (Semi-Global Block-Matching)

It is a wrapper around the OpenCV implementation of semi-global block-matching,
which is a variant of Hirschmuller famous semi-global matching (SGM). You must
first compile and install a few OpenCV modules:

    git clone https://github.com/Itseez/opencv.git
    cd opencv
    git checkout 2.4
    mkdir build
    cd build
    cmake -D CMAKE_INSTALL_PREFIX=~/local -D CMAKE_BUILD_TYPE=RELEASE -D BUILD_opencv_ml=OFF -D BUILD_opencv_objdetect=OFF -D BUILD_opencv_video=OFF -D BUILD_opencv_photo=OFF -D BUILD_opencv_nonfree=OFF -D BUILD_opencv_java=OFF -D BUILD_opencv_ts=OFF ..
    make
    make install

Please note the importance of compiling the code of the 2.4 branch, and not the
master branch, since the OpenCV API may change and break the compatibility with
our wrapper.

Now you can compile the SGBM wrapper:

    cd 3rdparty/stereo_hirschmuller_2008
    mkdir build
    cd build
    cmake -D CMAKE_PREFIX_PATH=~/local ..
    make

#### numpy >= 1.6.1

    cd 3rdparty
    wget http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1.tar.gz
    tar xzf numpy-1.6.1.tar.gz
    cd numpy-1.6.1
    python setup.py build --fcompiler=gnu95 <!---check the README.txt of numpy for the compiler options-->
    python setup.py install --prefix=~/local/

#### piio

    cd 3rdparty
    git clone https://github.com/carlodef/iio.git
    cd piio_packaged
    python setup.py install --prefix=~/local
    export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.6/site-packages

on python 2.7, linux x86_64 (Fedora)

    export PYTHONPATH=$PYTHONPATH:$HOME/local/lib64/python2.7/site-packages

#### Copy 3rd party binaries into bin

To make the 3rd party binaries available to the s2p system run the following script

    cd bin
    . copy_from_3rdparty.sh

If GDAL or Geographic LIB was installed locally then also run:

    cp ~/local/bin/{CartConvert,GeoidEval,gdal_translate} .


## Satellite images datasets

Each image must be accompanied by an `xml` file containing rpc coefficients.

Due to storage limitations, no images are available on this repository.  If you
want to run the s2p code, **you have to get a copy of Pleiades, WorldView of
QuickBird datasets by other means.** The typical size of a Pléiades image is
around 40000 x 40000 pixels, covering an area of 20km x 20km. The files weigh
approximately 2GB each.

Run the script `data/pleiades/create_links.sh` to generate symbolic links to
the images files and copy the `xml` files. The links will be located in the
`data/pleiades/*` subfolders and are easier to use than the actual real paths
to the image files.
