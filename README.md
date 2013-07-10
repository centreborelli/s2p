# S2P - Satellite Stereo Pipeline

This code implements a stereo pipeline for producing elevation models from
Pleiades satellite images. It aims at automatically generating digital
elevation models from tristereo sets (three images).


## Pleiades data

Several Pleiades stereoscopic datasets are available. We have pairs and
triplets. Each image is accompanied by an xml file containing rpc coefficients.

All the data is stored in the `pleiades_data` folder, which contains two
subfolders, `rpc` and `images`.

## Needed binaries

The following binaries are used by the python routines and should be in your
PATH:
    srtm4
    rectify_mindistortion
    gdal_translate
    sift_keypoints
    siftu
    homography
    zoom_zeropadding
    fftconvolve
