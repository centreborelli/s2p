# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os
import ctypes
import warnings

import numpy as np
import rasterio as rio
from numpy.ctypeslib import ndpointer
import ransac

from s2p import rpc_utils
from s2p import estimation

# Locate sift4ctypes library and raise an ImportError if it can not be
# found This call will raise an exception if library can not be found,
# at import time

# TODO: This is kind of ugly. Cleaner way to do this is to update
# LD_LIBRARY_PATH, which we should do once we have a proper config file
here = os.path.dirname(os.path.abspath(__file__))
sift4ctypes = os.path.join(os.path.dirname(here), 'lib', 'libsift4ctypes.so')
lib = ctypes.CDLL(sift4ctypes)


# Filter warnings from rasterio reading files wihtout georeferencing
warnings.filterwarnings("ignore", category=rio.errors.NotGeoreferencedWarning)


def keypoints_from_nparray(arr, thresh_dog=0.0133, nb_octaves=8, nb_scales=3, offset=None):
    """
    Runs SIFT (the keypoints detection and description only, no matching) on an image stored in a 2D numpy array

    It uses Ives Rey Otero's implementation published in IPOL:
    http://www.ipol.im/pub/pre/82/

    Args:
        arr: A 2D numpy array respresenting the input image
        thresh_dog (optional): Threshold on gaussian derivative
        nb_octaves (optional): Number of octaves
        nb_scales (optional): Number of scales
        offset (optional): offset to apply to sift position in case arr is an extract of a bigger image

    Returns:
        A numpy array of shape (nb_points,132) containing for each row (y,x,scale,orientation, sift_descriptor)
    """
    # retrieve numpy buffer dimensions
    h, w = arr.shape

    # Set expected args and return types
    lib.sift.argtypes = (ndpointer(dtype=ctypes.c_float, shape=(h, w)), ctypes.c_uint, ctypes.c_uint, ctypes.c_float,
                         ctypes.c_uint, ctypes.c_uint, ctypes.POINTER(ctypes.c_uint), ctypes.POINTER(ctypes.c_uint))
    lib.sift.restype = ctypes.POINTER(ctypes.c_float)

    # Create variables to be updated by function call
    nb_points = ctypes.c_uint()
    desc_size = ctypes.c_uint()

    # Call sift fonction from sift4ctypes.so
    keypoints_ptr = lib.sift(arr.astype(np.float32), w, h, thresh_dog,
                             nb_octaves, nb_scales, ctypes.byref(desc_size), ctypes.byref(nb_points))

    # Transform result into a numpy array
    keypoints = np.asarray([keypoints_ptr[i]
                            for i in range(nb_points.value*desc_size.value)])

    # Delete results to release memory
    lib.delete_buffer.argtypes = (ctypes.POINTER(ctypes.c_float)),
    lib.delete_buffer(keypoints_ptr)

    # Reshape keypoints array
    keypoints = keypoints.reshape((nb_points.value, desc_size.value))

    if offset is not None:
        x, y = offset
        keypoints[:, 0] += x
        keypoints[:, 1] += y

    return keypoints


def image_keypoints(im, x, y, w, h, max_nb=None, thresh_dog=0.0133, nb_octaves=8, nb_scales=3):
    """
    Runs SIFT (the keypoints detection and description only, no matching).

    It uses Ives Rey Otero's implementation published in IPOL:
    http://www.ipol.im/pub/pre/82/

    Args:
        im (str): path to the input image
        max_nb (optional): maximal number of keypoints. If more keypoints are
            detected, those at smallest scales are discarded

    Returns:
        numpy array of shape (n, 132) containing, on each row: (y, x, s, o, 128-descriptor)
    """
    # Read file with rasterio
    with rio.open(im) as ds:
        # clip roi to stay inside the image boundaries
        if x < 0:  # if x is negative then replace it with 0 and reduce w
            w += x
            x = 0
        if y < 0:
            h += y
            y = 0
        # if extract not completely inside the full image then resize (w, h)
        w = min(w, ds.width - x)
        h = min(h, ds.height - y)
        in_buffer = ds.read(window=rio.windows.Window(x, y, w, h))

    # Detect keypoints on first band
    keypoints = keypoints_from_nparray(in_buffer[0], thresh_dog=thresh_dog,
                                       nb_octaves=nb_octaves,
                                       nb_scales=nb_scales, offset=(x, y))

    # Limit number of keypoints if needed
    if max_nb is not None:
        keypoints = keypoints[:max_nb]

    return keypoints


def string_dump_of_keypoint_and_descriptor(k):
    """
    Return a string representing a keypoint and its descriptor.

    Args:
        k (array_like): list of 132 floats, the first four elements are the
            keypoint (x, y, scale, orientation), the 128 following elements are
            the coefficients of the SIFT descriptor and take integer values
            between 0 and 255.

    Return:
        string dump of the descriptor, such as for example
        "342.254 003.570 0.91346 2.36788 000 001 005 000 000 000 028 029 179..."
    """
    s = "{:8.3f} {:8.3f} {:7.3f} {: 5.3f} ".format(*k[:4])
    s += " ".join("{:3d}".format(int(x)) for x in k[4:])
    return s


def keypoints_match(k1, k2, method='relative', sift_thresh=0.6, F=None,
                    epipolar_threshold=10, model=None, ransac_max_err=0.3):
    """
    Find matches among two lists of sift keypoints.

    Args:
        k1 (array): numpy array of shape (n, 132), where each row represents a
            sift keypoint with (y, x, scale, orientation, 128-descriptor)
        k2 (array): numpy array of shape (m, 132), where each row represents a
            sift keypoint
        method (optional, default is 'relative'): flag ('relative' or
            'absolute') indicating whether to use absolute distance or relative
            distance
        sift_thresh (optional, default is 0.6): threshold for distance between SIFT
            descriptors. These descriptors are 128-vectors, whose coefficients
            range from 0 to 255, thus with absolute distance a reasonable value
            for this threshold is between 200 and 300. With relative distance
            (ie ratio between distance to nearest and distance to second
            nearest), the commonly used value for the threshold is 0.6.
        F (optional): affine fundamental matrix
        epipolar_threshold (optional, default is 10): maximum distance allowed for
            a point to the epipolar line of its match.
        model (optional, default is None): model imposed by RANSAC when
            searching the set of inliers. If None all matches are considered as
            inliers.
        ransac_max_err (float): maximum allowed epipolar error for
            RANSAC inliers. Optional, default is 0.3.

    Returns:
        if any, a numpy 2D array containing the list of inliers matches.
    """
    # compute matches
    matches = keypoints_match_from_nparray(k1, k2, method, sift_thresh,
                                           epipolar_threshold, F)

    # filter matches with ransac
    if model == 'fundamental' and len(matches) >= 7:
        inliers = ransac.find_fundamental_matrix(matches, ntrials=1000,
                                                 max_err=ransac_max_err)[0]
        matches = matches[inliers]

    return matches


def keypoints_match_from_nparray(k1, k2, method, sift_threshold,
                                 epi_threshold=10, F=None):
    """
    Wrapper for the sift keypoints matching function of libsift4ctypes.so.
    """
    # Set expected args and return types
    lib.matching.argtypes = (ndpointer(dtype=ctypes.c_float, shape=k1.shape),
                             ndpointer(dtype=ctypes.c_float, shape=k2.shape),
                             ctypes.c_uint, ctypes.c_uint, ctypes.c_uint,
                             ctypes.c_uint, ctypes.c_float, ctypes.c_float,
                             ndpointer(dtype=ctypes.c_double, shape=(5,)),
                             ctypes.c_bool, ctypes.c_bool,
                             ctypes.POINTER(ctypes.c_uint))
    lib.matching.restype = ctypes.POINTER(ctypes.c_float)

    # Get info of descriptor size
    nb_sift_k1, descr = k1.shape
    sift_offset = 4
    length_descr = descr - sift_offset

    # Transform information of method into boolean
    use_relative_method = (method == 'relative')

    # Format fundamental matrix
    use_fundamental_matrix = False
    coeff_mat = np.zeros(5)
    if F is not None:
        coeff_mat = np.asarray([F[0, 2], F[1, 2], F[2, 0], F[2, 1], F[2, 2]])
        use_fundamental_matrix = True

    # Create variables to be updated by function call
    nb_matches = ctypes.c_uint()

    # Call sift fonction from sift4ctypes.so
    matches_ptr = lib.matching(k1.astype('float32'), k2.astype('float32'),
                               length_descr, sift_offset, len(k1), len(k2),
                               sift_threshold, epi_threshold, coeff_mat,
                               use_fundamental_matrix, use_relative_method,
                               ctypes.byref(nb_matches))

    # Transform result into a numpy array
    matches = np.asarray([matches_ptr[i] for i in range(nb_matches.value * 4)])

    # Delete results to release memory
    lib.delete_buffer.argtypes = ctypes.POINTER(ctypes.c_float),
    lib.delete_buffer(matches_ptr)

    # Reshape keypoints array
    return matches.reshape((nb_matches.value, 4))


def matches_on_rpc_roi(cfg, im1, im2, rpc1, rpc2, x, y, w, h,
                       method, sift_thresh, epipolar_threshold):
    """
    Compute a list of SIFT matches between two images on a given roi.

    The corresponding roi in the second image is determined using the rpc
    functions.

    Args:
        im1, im2: paths to two large tif images
        rpc1, rpc2: two instances of the rpcm.RPCModel class
        x, y, w, h: four integers defining the rectangular ROI in the first
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        method, sift_thresh, epipolar_threshold: see docstring of
            s2p.sift.keypoints_match()

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the full images.
    """
    x2, y2, w2, h2 = rpc_utils.corresponding_roi(cfg, rpc1, rpc2, x, y, w, h)

    # estimate an approximate affine fundamental matrix from the rpcs
    rpc_matches = rpc_utils.matches_from_rpc(cfg, rpc1, rpc2, x, y, w, h, 5)
    F = estimation.affine_fundamental_matrix(rpc_matches)

    # if less than 10 matches, lower thresh_dog. An alternative would be ASIFT
    thresh_dog = 0.0133
    for _ in range(2):
        p1 = image_keypoints(im1, x, y, w, h, thresh_dog=thresh_dog)
        p2 = image_keypoints(im2, x2, y2, w2, h2, thresh_dog=thresh_dog)
        matches = keypoints_match(p1, p2, method, sift_thresh, F,
                                  epipolar_threshold=epipolar_threshold,
                                  model='fundamental')
        if matches is not None and matches.ndim == 2 and matches.shape[0] > 10:
            break
        thresh_dog /= 2.0
    else:
        print("WARNING: sift.matches_on_rpc_roi: found no matches.")
        return None
    return matches
