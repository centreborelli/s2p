# S2P - Satellite Stereo Pipeline

This code implements a stereo pipeline for producing elevation models from
Pleiades satellite images. It aims at automatically generating digital
elevation models from tristereo sets (three images).

The main language is Python, although several operations are handled by
binaries written in C.

The pipeline is implemented in the file `s2p.py`. The `s2p` module can be
used to produce elevation models and 3D point clouds of arbitrarily large
regions of interest. If needed, it cuts the region of interest in several
small tiles and process them in parallel.

## Usage

Import `s2p` the module in a python session, and run the functions
`process_pair` or `process_triplet`, depending on the kind of dataset you have
(stereo pair or triplet)

    python
    >>> import s2p
    >>> s2p.process_pair('s2p_test', 'toulouse', 2, 1, 16000, 12000, 300, 300)

See the docstrings of the functions `process_pair` and `process_triplet` for a
complete description of their arguments.

## Installation

All the python modules are located in the `python` folder.  Some python
functions of these modules rely on external binaries. Most of these binaries
were written on purpose for the needs of the pipeline, and their source code is
provided here in the `c` folder.

### S2P binaries

The source code is in the `c` folder. To compile it:

    cd c
    make

This will create a `bin` directory containing all the s2p binaries.

### 3rd party binaries

A few other binaries are 3rd party. Their source code is located in the
`3rdparty` directory. You must compile them before launching the script, using
the provided makefiles. For example, for GeographicLib, do:

    cd 3rdparty/GeographicLib-1.32
    make
    sudo make install

Since we use GeographicLib to evaluate geoid heights we must also install the
geoids data files by running the script:

    3rdparty/GeographicLib-1.32/tools/geographiclib-get-geoids.sh



For SGBM (Semi-Global Block-Matching), do:

    cd 3rdparty/stereo_hirschmuller_2008
    mkdir build
    cd build
    cmake ..
    make

This binary uses OpenCV implementation of Hirschmuller Semi-Global Matching.
You must have OpenCV 2.4.x installed on your system to compile it.

    git clone https://github.com/Itseez/opencv.git
    cd opencv
    git checkout 2.4
    mkdir build_2.4
    cd build_2.4
    cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=~/local ..
    make

For sift, do:

    cd 3rdparty/sift_20130403
    make

In addition, the `gdal_translate` binary is needed, and can be installed
through your favourite package manager.


## Pleiades data

Several Pleiades stereoscopic datasets are available. We have pairs and
triplets. Each image is accompanied by an `xml` file containing rpc
coefficients.

Due to storage limitations, the images are not available on this repository.
If you want to run the s2p code, **you have to get a copy of our Pleiades
dataset by other means.** The size of these images is around 40000 x 40000
pixels, covering an area of 20km x 20km. The files weigh approximately 2GB
each.

Only the `xml` files containing the calibration data (encoded by rpc
coefficients) are provided. They are located in the folder `pleiades_data/rpc`.
The folder `pleiades_data/images` contains the list of the relative paths to
the full images of our dataset. These paths are relative to the location of
your copy of the dataset. Run the script `create_links.sh` to generate symbolic
links to the images files. These links will be located in the
`pleiades_data/images/*` subfolders and are easier to use than the actual real
paths to the image files.
