#!/usr/bin/env python

import numpy as np
import rectification
import rpc_utils
import common

def pointing_accuracy_evaluation(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Measures the maximal pointing error on a Pleiades' pair of images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.

    Returns:
        the highest pointing error, in the direction orthogonal to the epipolar
        lines. This error is measured in pixels.
    """
    # get sift matches
    matches = rectification.matches_from_sift(im1, im2, rpc1, rpc2, x, y, w, h)

    # filter outliers with ransac
    # the binary is from Enric's imscript
    matches_file = common.tmpfile('.txt')
    inliers_file = common.tmpfile('.txt')
    np.savetxt(matches_file, matches)
    common.run("ransac fmn 1000 1 7 %s < %s" % (inliers_file, matches_file))
    inliers = np.loadtxt(inliers_file)
    if not inliers.size:
        raise Exception("no inliers")

    # compute the pointing error for each match
    n = np.shape(inliers)[0]
    e = np.zeros((n, 1))
    x1 = inliers[:, 0]
    y1 = inliers[:, 1]
    x2 = inliers[:, 2]
    y2 = inliers[:, 3]
    for i in range(n):
        e[i] = rpc_utils.compute_height(rpc1, rpc2, x1[i], y1[i], x2[i],
            y2[i])[1]
    print "max, mean, min pointing error, from %d points:" % (n)
    print np.max(e), np.mean(e), np.min(e)

    # return the highest one
    return np.max(np.abs(e))
