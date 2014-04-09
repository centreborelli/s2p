# S2P - Satellite Stereo Pipeline

This code implements a stereo pipeline which produces elevation models from
images taken by high resolution satellites such as Pléiades, WorldView and
QuickBird. It generates automatically digital elevation models from stereo
pairs (two images) or tri-stereo sets (three images).

The main language is Python, although several operations are handled by
binaries written in C.

The pipeline is implemented in the file `s2p.py`. The `s2p` module can be used
to produce elevation models and 3D point clouds from arbitrarily large regions
of interest or from complete images. If needed, it cuts the region of interest
in several small tiles and process them in parallel.

## Usage

The easiest way to use `s2p` is to run the binary `s2p.py` from a shell, with a
json configuration file as unique argument:

    $ ./s2p.py config.json

All the parameters of the algorithm, paths to input and output data are stored
in the json file. See the provided `config.json.example` file for an example.

An other way is to import the `s2p` module in a python session, and run the
functions `process_pair` or `process_triplet`, depending on the kind of dataset
you have (stereo pair or triplet), to generate a digital elevation model (DEM),
and then generate a 3D point cloud from this DEM using the function
`generate_cloud`.

    python
    >>> import s2p
    >>> dem = s2p.process_triplet('test', 'pan1.tif', 'rpc1.xml', 'pan2.tif', 'rpc2.xml', 'pan3.tif', 'rpc3.xml', 25150, 24250, 300, 300, 3)
    >>> s2p.generate_cloud('test', 'pan1.tif', 'rpc1.xml', 'xs1.tif', 25150, 24250, 300, 300, dem)

See the docstrings of the functions `process_pair`, `process_triplet` and
`generate_cloud` for a complete description of their arguments.

## Installation

All the python modules are located in the `python` folder. Some python
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
`3rdparty` directory. You must compile and install them before launching the
script, using the provided makefiles. Here are the instuctions to install them
in a local folder, without root privileges. First create a local folder:

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

#### GDAL >= 1.10

    cd 3rdparty
    wget http://download.osgeo.org/gdal/1.10.1/gdal-1.10.1.tar.gz
    tar xzf gdal-1.10.1.tar.gz
    cd gdal-1.10.1
    ./configure --prefix=$HOME/local
    make install

#### Sift

    cd 3rdparty/sift_20130403
    make

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

#### piio

    cd 3rdparty
    git clone https://github.com/carlodef/iio.git
    cd piio_packaged
    python setup.py install --prefix=~/local
    export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python2.6/site-packages

#### numpy >= 1.6.1

    cd 3rdparty
    wget http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1.tar.gz
    tar xzf numpy-1.6.1.tar.gz
    cd numpy-1.6.1
    python setup.py build --fcompiler=gnu95 <!---check the README.txt of numpy for the compiler options-->
    python setup.py install --prefix=~/local/


### Copy 3rd party binaries into bin

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
