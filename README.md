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


## Pleiades data

Several Pleiades stereoscopic datasets are available. We have pairs and
triplets. Each image is accompanied by an `xml` file containing rpc
coefficients.

All the data is stored in the `pleiades_data` folder, which contains two
subfolders, `rpc` and `images`.

## Dependencies

All the binaries called by the python script are either located in the `bin` or
the `3rdparty` directories. You must compile them before launching the script,
using the provided makefiles.

In addition, the following binaries must be available on your system:

    gdal_translate
    tiffcp
    java
