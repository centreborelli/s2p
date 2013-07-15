# S2P - Satellite Stereo Pipeline

This code implements a stereo pipeline for producing elevation models from
Pleiades satellite images. It aims at automatically generating digital
elevation models from tristereo sets (three images).

The main language is Python, although several operations are handled by
binaries written in C.

## Usage

The pipeline is implemented in the file `main_script_pairs.py`. Its parameters
(paths to Pleiades data, definition of a region of interest, paths for output
data), are defined at the beginning of that file. To launch it, simply run:

    ./main_script_pairs.py



## Dependencies

Some python functions of the S2P modules rely on external binaries. Most of
these binaries were written on purpose for the needs of the pipeline, and their
source code is provided here in the `c` folder.

### S2P binaries

The source code is in the `c` folder. To compile it:

    cd c
    make

This will create a `bin` directory containing all the S2P binaries.

### 3rd party binaries

A few other binaries are 3rd party. Their source code is located in the
`3rdparty` directory. You must compile them before launching the script, using
the provided makefiles. For example, for GeographicLib, do:

    cd 3rdparty/GeographicLib-1.32
    make
    sudo make install

In addition, the following binaries must be available on your system:

    gdal_translate
    tiffcp
    java

You can install them through a package manager.


## Pleiades data

Several Pleiades stereoscopic datasets are available. We have pairs and
triplets. Each image is accompanied by an `xml` file containing rpc
coefficients.

All the data is stored in the `pleiades_data` folder, which contains two
subfolders, `rpc` and `images`.
