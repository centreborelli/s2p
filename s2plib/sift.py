# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>


from __future__ import print_function
import os
import numpy as np

import rasterio as rio
import ctypes
from ctypes import c_uint, c_float, byref, POINTER
from numpy.ctypeslib import ndpointer

from s2plib import common
from s2plib import rpc_utils
from s2plib import estimation
from s2plib.config import cfg


def image_keypoints(im, x, y, w, h, max_nb=None, extra_params=''):
    """
    Runs SIFT (the keypoints detection and description only, no matching).

    It uses Ives Rey Otero's implementation published in IPOL:
    http://www.ipol.im/pub/pre/82/

    Args:
        im: path to the input image
        max_nb (optional): maximal number of keypoints. If more keypoints are
            detected, those at smallest scales are discarded
        extra_params (optional): extra parameters to be passed to the sift
            binary

    Returns:
        path to the file containing the list of descriptors
    """
    # Read file with rasterio
    with rio.open(im) as ds:
       in_buffer = ds.read(window=((x,x+w),(y,y+h)))

    # load shared library 
    lib = ctypes.CDLL('sift4ctypes.so')

    # retrieve numpy buffer dimensions  
    (band,h,w) = in_buffer.shape

    lib.sift.argtypes = (ndpointer(dtype=c_float, shape=(band,w,h)),c_uint, c_uint, c_float,c_uint, c_uint,ctypes.POINTER(c_uint),ctypes.POINTER(c_uint))
    lib.sift.restype = POINTER(c_float)
    nb_points = c_uint()
    desc_size = c_uint()
    keypoints_ptr = lib.sift(in_buffer.astype(np.float32),w,h,0.0133,8,3,byref(desc_size),byref(nb_points))

    

    print("Number of points:{}".format(nb_points))
    print("Descriptor size:{}".format(desc_size))
    print(repr(keypoints_ptr))

    keypoints = numpy.array(keypoints_ptr.contents,(nb_points.value,desc_size.value))
    #keypoints = keypoints.reshape((nb_points.value,desc_size.value))

   # print(keypoints.shape())

    keyfile = common.tmpfile('.txt')
    
    return keyfile


def keypoints_match(k1, k2, method='relative', sift_thresh=0.6, F=None,
                    model=None, epipolar_threshold=10):
    """
    Find matches among two lists of sift keypoints.

    Args:
        k1, k2: paths to text files containing the lists of sift descriptors
        method (optional, default is 'relative'): flag ('relative' or
            'absolute') indicating wether to use absolute distance or relative
            distance
        sift_thresh (optional, default is 0.6): threshold for distance between SIFT
            descriptors. These descriptors are 128-vectors, whose coefficients
            range from 0 to 255, thus with absolute distance a reasonable value
            for this threshold is between 200 and 300. With relative distance
            (ie ratio between distance to nearest and distance to second
            nearest), the commonly used value for the threshold is 0.6.
        F (optional): affine fundamental matrix
        model (optional, default is None): model imposed by RANSAC when
            searching the set of inliers. If None all matches are considered as
            inliers.
        epipolar_threshold (optional, default is 10): maximum distance allowed for
            a point to the epipolar line of its match.

    Returns:
        if any, a numpy 2D array containing the list of inliers matches.
    """
    # compute matches
    mfile = common.tmpfile('.txt')
    cmd = "matching %s %s -o %s --sift-threshold %f" % (k1, k2, mfile, sift_thresh)
    if method == 'absolute':
        cmd += " --absolute"
    if F is not None:
        fij = ' '.join(str(x) for x in [F[0, 2], F[1, 2], F[2, 0],
                                        F[2, 1], F[2, 2]])
        cmd = "%s -f \"%s\"" % (cmd, fij)
        cmd += " --epipolar-threshold {}".format(epipolar_threshold)
    common.run(cmd)

    matches = np.loadtxt(mfile)
    if matches.ndim == 2:  # filter outliers with ransac
        if model == 'fundamental' and len(matches) >= 7:
            common.run("ransac fmn 1000 .3 7 %s < %s" % (mfile, mfile))
        elif model == 'homography' and len(matches) >= 4:
            common.run("ransac hom 1000 1 4 /dev/null /dev/null %s < %s" % (mfile,
                                                                            mfile))
        elif model == 'hom_fund' and len(matches) >= 7:
            common.run("ransac hom 1000 2 4 /dev/null /dev/null %s < %s" % (mfile,
                                                                            mfile))
            common.run("ransac fmn 1000 .2 7 %s < %s" % (mfile, mfile))

    if os.stat(mfile).st_size > 0:  # return numpy array of matches
        return np.loadtxt(mfile)


def matches_on_rpc_roi(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Compute a list of SIFT matches between two images on a given roi.

    The corresponding roi in the second image is determined using the rpc
    functions.

    Args:
        im1, im2: paths to two large tif images
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers defining the rectangular ROI in the first
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the full images.
    """
    x2, y2, w2, h2 = rpc_utils.corresponding_roi(rpc1, rpc2, x, y, w, h)

    # estimate an approximate affine fundamental matrix from the rpcs
    rpc_matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, 5)
    F = estimation.affine_fundamental_matrix(rpc_matches)

    # sift matching method:
    method = 'relative' if cfg['relative_sift_match_thresh'] is True else 'absolute'

    # if less than 10 matches, lower thresh_dog. An alternative would be ASIFT
    thresh_dog = 0.0133
    for i in range(2):
        p1 = image_keypoints(im1, x, y, w, h, extra_params='--thresh-dog %f' % thresh_dog)
        p2 = image_keypoints(im2, x2, y2, w2, h2, extra_params='--thresh-dog %f' % thresh_dog)
        matches = keypoints_match(p1, p2, method, cfg['sift_match_thresh'],
                                  F, model='fundamental',
                                  epipolar_threshold=cfg['max_pointing_error'])
        if matches is not None and matches.ndim == 2 and matches.shape[0] > 10:
            break
        thresh_dog /= 2.0
    else:
        print("WARNING: sift.matches_on_rpc_roi: found no matches.")
        return None
    return matches
