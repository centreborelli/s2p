# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from __future__ import print_function
import os
import sys
import shutil
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

from s2plib import common
from s2plib import piio
from s2plib.config import cfg

def register_heights(im1, im2):
    """
    Affine registration of heights.

    Args:
        im1: first height map
        im2: second height map, to be registered on the first one

    Returns
        path to the registered second height map
    """
    # remove high frequencies with a morphological zoom out
    im1_low_freq = common.image_zoom_out_morpho(im1, 4)
    im2_low_freq = common.image_zoom_out_morpho(im2, 4)

    # first read the images and store them as numpy 1D arrays, removing all the
    # nans and inf
    i1 = piio.read(im1_low_freq).ravel() #np.ravel() gives a 1D view
    i2 = piio.read(im2_low_freq).ravel()
    ind = np.logical_and(np.isfinite(i1), np.isfinite(i2))
    h1 = i1[ind]
    h2 = i2[ind]

    # for debug
    print(np.shape(i1))
    print(np.shape(h1))

#    # 1st option: affine
#    # we search the (u, v) vector that minimizes the following sum (over
#    # all the pixels):
#    #\sum (im1[i] - (u*im2[i]+v))^2
#    # it is a least squares minimization problem
#    A = np.vstack((h2, h2*0+1)).T
#    b = h1
#    z = np.linalg.lstsq(A, b)[0]
#    u = z[0]
#    v = z[1]
#
#    # apply the affine transform and return the modified im2
#    out = common.tmpfile('.tif')
#    common.run('plambda %s "x %f * %f +" > %s' % (im2, u, v, out))

    # 2nd option: translation only
    v = np.mean(h1 - h2)
    out = common.tmpfile('.tif')
    common.run('plambda %s "x %f +" -o %s' % (im2, v, out))

    return out


def merge(im1, im2, thresh, out, conservative=False):
    """
    Args:
        im1, im2: paths to the two input images
        thresh: distance threshold on the intensity values
        out: path to the output image
        conservative (optional, default is False): if True, keep only the
            pixels where the two height map agree

    This function merges two images. They are supposed to be two height maps,
    sampled on the same grid. If a pixel has a valid height (ie not inf) in
    only one of the two maps, then we keep this height (if the 'conservative'
    option is set to False). When two heights are available, if they differ
    less than the threshold we take the mean, if not we discard the pixel (ie
    assign NAN to the output pixel).
    """
    # first register the second image on the first
    im2 = register_heights(im1, im2)

    if conservative:
        # then merge
        # the following plambda expression implements:
        # if isfinite x
        #   if isfinite y
        #     if fabs(x - y) < t
        #       return (x+y)/2
        #     return nan
        #   return nan
        # return nan
        common.run("""
            plambda %s %s "x isfinite y isfinite x y - fabs %f < x y + 2 / nan if nan
            if nan if" -o %s
            """ % ( im1, im2, thresh, out))
    else:
        # then merge
        # the following plambda expression implements:
        # if isfinite x
        #   if isfinite y
        #     if fabs(x - y) < t
        #       return (x+y)/2
        #     return nan
        #   return x
        # return y
        common.run("""
            plambda %s %s "x isfinite y isfinite x y - fabs %f < x y + 2 / nan if x
            if y if" -o %s
            """ % ( im1, im2, thresh, out))


def average_if_close(x, threshold):
    """
    """
    if np.nanmax(x) - np.nanmin(x) > threshold:
        return np.nan
    else:
        return np.nanmedian(x)


def merge_n(output, inputs, offsets, averaging='average_if_close', threshold=1):
    """
    Merge n images of equal sizes by taking the median/mean/min/max pixelwise.

    Args:
        inputs: list of paths to the input images
        output: path to the output image
        averaging: string containing the name of a function that accepts
            1D arrays. It is applied to 1D slices of the stack of images along
            the last axis. Possible values are, for instance np.min, np.max,
            np.mean, np.median and their nanproof counterparts, ie np.nanmin,
            np.nanmax, np.nanmean, np.nanmedian
    """
    assert(len(inputs) == len(offsets))

    # get input images size
    if inputs:
        f = gdal.Open(inputs[0])
        w, h = f.RasterXSize, f.RasterYSize
        f = None  # this is the gdal way of closing files

    # read input images and apply offsets
    x = np.empty((h, w, len(inputs)))
    for i, img in enumerate(inputs):
        f = gdal.Open(img)
        x[:, :, i] = f.GetRasterBand(1).ReadAsArray() - offsets[i]
        f = None
        if cfg['debug']:
            piio.write('{}_registered.tif'.format(os.path.splitext(img)[0]),
                       x[:, :, i] + np.mean(offsets))

    # apply the averaging operator
    if averaging.startswith(('np.', 'numpy.')):
        avg = np.apply_along_axis(getattr(sys.modules['numpy'], averaging.split('.')[1]),
                                  axis=2, arr=x)
    elif averaging == 'average_if_close':
        avg = np.apply_along_axis(average_if_close, 2, x, threshold)

    # add the mean offset
    avg += np.mean(offsets)

    # write the average to output
    if inputs:
        shutil.copy(inputs[0], output)  # copy an input file to get the metadata
        f = gdal.Open(output, gdal.GA_Update)
        f.GetRasterBand(1).WriteArray(avg)  # update the output file content
        f = None
