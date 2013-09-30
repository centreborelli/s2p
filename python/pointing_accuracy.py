#!/usr/bin/env python

import numpy as np
import rectification
import rpc_utils
import common


def evaluation(im1, im2, rpc1, rpc2, x, y, w, h):
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
    matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, w, h)

    # compute the pointing error for each match
    n = np.shape(matches)[0]
    e = np.zeros((n, 1))
    x1 = matches[:, 0]
    y1 = matches[:, 1]
    x2 = matches[:, 2]
    y2 = matches[:, 3]
    for i in range(n):
        e[i] = rpc_utils.compute_height(rpc1, rpc2, x1[i], y1[i], x2[i],
            y2[i])[1]
    print "max, mean, min pointing error, from %d points:" % (n)
    print np.max(e), np.mean(e), np.min(e)

    # return the highest one
    return np.max(np.abs(e))


def filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the big images.
            If no sift matches are found, then an exception is raised.

    The returned matches are the inliers of an epipolar model found with ransac.
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

    return inliers


def euclidean_transform_matrix(theta, a, b):
    """
    Build the matrix of an euclidean transform in homogeneous coordinates.

    Arguments:
        theta: angle of the rotation
        a, b: coordinates of the vector of the translation

    Returns:
        A numpy 3x3 matrix, representing the euclidean transform (rotation
        followed by translation) in homogeneous coordinates.
    """
    R = np.eye(3)
    R[0, 0] =  np.cos(theta)
    R[0, 1] =  np.sin(theta)
    R[1, 0] = -np.sin(theta)
    R[1, 1] =  np.cos(theta)
    T = np.eye(3)
    T[0, 2] = a
    T[1, 2] = b
    return np.dot(T, R)


def cost_function(v, rpc1, rpc2, matches, alpha=0.01):
    """
    Objective function to minimize in order to correct the pointing error.

    Arguments:
        v: vector of size 3, containing the 3 parameters of the euclidean
            transformation we are looking for.
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is the one of the big images.
        alpha: relative weight of the error terms: e + alpha*(h-h0)^2. See
            paper for more explanations.

    Returns:
        The sum of pointing errors and altitude differences, as written in the
        paper formula (1).
    """
    # compute the altitudes from the matches without correction
    n = np.shape(matches)[0]
    h0 = np.zeros((n, 1))
    x1 = matches[:, 0]
    y1 = matches[:, 1]
    x2 = matches[:, 2]
    y2 = matches[:, 3]
    for i in range(n):
        h0[i] = rpc_utils.compute_height(rpc1, rpc2, x1[i], y1[i], x2[i],
            y2[i])[0]

    # transform the coordinates of points in the second image according to
    # matrix A, built from vector v
    A = euclidean_transform_matrix(v[0], v[1], v[2])
    p2 = common.points_apply_homography(A, matches[:, 2:4])
    x2 = p2[:, 0]
    y2 = p2[:, 1]

    # compute the cost
    cost = 0
    for i in range(n):
        h, e = rpc_utils.compute_height(rpc1, rpc2, x1[i], y1[i], x2[i], y2[i])
        cost += e
        cost += alpha * (h - h0[i])**2
    
    print cost
    return cost


def filtered_sift_matches_full_img(im1, im2, rpc1, rpc2, a=1000):
    """
    Computes a list of sift matches between two full Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        a: length of the squared ROIs used to extract sift points

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the big images.
            If no sift matches are found, then an exception is raised.

    The keypoints are extracted from five zones in the first image. One in the
    center, and four in the corners. The matches are consistent with an
    epipolar model in each of the zones.
    """
    w = rpc1.lastCol
    h = rpc1.lastRow

    # if roi size is too big wrt image size, take a smaller roi size
    if (min(h, w) < 4*a):
        a = round(min(h, w)/4)

    # central zone
    x = round((w-a)/2)
    y = round((h-a)/2)
    matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, a, a)

    # corner zones
    x = round((1*w - 2*a)/4)
    y = round((1*h - 2*a)/4)
    matches = np.vstack((matches, filtered_sift_matches_roi(im1, im2, rpc1,
        rpc2, x, y, a, a)))

    x = round((3*w - 2*a)/4)
    y = round((1*h - 2*a)/4)
    matches = np.vstack((matches, filtered_sift_matches_roi(im1, im2, rpc1,
        rpc2, x, y, a, a)))

    x = round((1*w - 2*a)/4)
    y = round((3*h - 2*a)/4)
    matches = np.vstack((matches, filtered_sift_matches_roi(im1, im2, rpc1,
        rpc2, x, y, a, a)))

    x = round((3*w - 2*a)/4)
    y = round((3*h - 2*a)/4)
    matches = np.vstack((matches, filtered_sift_matches_roi(im1, im2, rpc1,
        rpc2, x, y, a, a)))

    return matches




def optimize_pair(im1, im2, rpc1, rpc2):

    matches = filtered_sift_matches_full_img(im1, im2, rpc1, rpc2)

    from scipy.optimize import fmin_bfgs
    v0 = np.zeros((1, 3))
    v = fmin_bfgs(cost_function, v0, args=(rpc1, rpc2, matches), maxiter=5,
        retall=True)
    return v
