#!/usr/bin/env python

import numpy as np
import rectification
import rpc_utils
import rpc_model
import common
import os


def evaluation(im1, im2, rpc1, rpc2, x, y, w, h, A=None):
    """
    Measures the maximal pointing error on a Pleiades' pair of images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2.

    Returns:
        the highest pointing error, in the direction orthogonal to the epipolar
        lines. This error is measured in pixels.
    """
    matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, w, h)
    p1 = matches[:, 0:2]
    p2 = matches[:, 2:4]

    # apply pointing correction matrix, if available
    if A is not None:
        p2 = common.points_apply_homography(A, p2)

    # compute the pointing error for each match
    x1 = p1[:, 0]
    y1 = p1[:, 1]
    x2 = p2[:, 0]
    y2 = p2[:, 1]
    e = rpc_utils.compute_height(rpc1, rpc2, x1, y1, x2, y2)[1]
    print "max, mean, min pointing error, from %d points:" % (len(matches))
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
    if len(matches) < 7:
        raise Exception("less than 7 matches")
    matches_file = common.tmpfile('.txt')
    np.savetxt(matches_file, matches)

    inliers_file = common.tmpfile('.txt')
    common.run("ransac fmn 1000 1 7 %s < %s" % (inliers_file, matches_file))
    inliers = np.loadtxt(inliers_file)
    if not inliers.size:
        raise Exception("no inliers")

    return inliers


def euclidean_transform_matrix(v):
    """
    Build the matrix of an euclidean transform in homogeneous coordinates.

    Arguments:
        v: numpy 1D array, of length 3 or 4, containing the parameters of the
            planar transform. These parameters are:
            v[0]: angle of the rotation
            v[1]: x coordinate of the vector of the translation
            v[2]: y coordinate of the vector of the translation
            v[3]: horizontal shear parameter

    Returns:
        A numpy 3x3 matrix, representing the euclidean transform (rotation
        followed by translation) in homogeneous coordinates.
    """
    R = np.eye(3)
    R[0, 0] =  np.cos(v[0])
    R[0, 1] =  np.sin(v[0])
    R[1, 0] = -np.sin(v[0])
    R[1, 1] =  np.cos(v[0])
    T = np.eye(3)
    T[0, 2] = v[1]
    T[1, 2] = v[2]
    S = np.eye(3)
    S[0, 1] = v[3]
    return np.dot(np.dot(T, R), S)


def cost_function(v, rpc1, rpc2, matches, alpha=0.01):
    """
    Objective function to minimize in order to correct the pointing error.

    Arguments:
        v: vector of size 4, containing the 4 parameters of the euclidean
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
    x1 = matches[:, 0]
    y1 = matches[:, 1]
    x2 = matches[:, 2]
    y2 = matches[:, 3]
    h0 = rpc_utils.compute_height(rpc1, rpc2, x1, y1, x2, y2)[0]

    # transform the coordinates of points in the second image according to
    # matrix A, built from vector v
    A = euclidean_transform_matrix(v)
    p2 = common.points_apply_homography(A, matches[:, 2:4])
    x2 = p2[:, 0]
    y2 = p2[:, 1]

    # compute the cost
    h, e = rpc_utils.compute_height(rpc1, rpc2, x1, y1, x2, y2)
    cost = np.sum((h - h0)**2)
    cost *= alpha
    cost += np.sum(e)

    #print cost
    return cost


def filtered_sift_matches_full_img(im1, im2, rpc1, rpc2, flag='automatic',
        prev1=None, a=1000):
    """
    Computes a list of sift matches between two full Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        flag: 'automatic', 'interactive' or 'load', to decide if the five zones
            used to search keypoints are queried interactively, chosen
            automatically, or loaded from the file 'pointing_correction_rois.txt'.
        prev1 (optional): path to the jpg preview image of im1 (used in case of
            interactive mode)
        a: length of the squared ROIs used to extract sift points, in the case
            of automatic mode

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
    out = np.zeros(shape = (1, 4))

    if flag == 'automatic':
        # if roi size is too big wrt image size, take a smaller roi size
        if (min(h, w) < 4*a):
            a = round(min(h, w)/4)

        # central zone
        x = round((w-a)/2)
        y = round((h-a)/2)
        try:
            matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the central zone"
            print e

        # corner zones
        x = round((1*w - 2*a)/4)
        y = round((1*h - 2*a)/4)
        try:
            matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 1"
            print e

        x = round((3*w - 2*a)/4)
        y = round((1*h - 2*a)/4)
        try:
            matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 2"
            print e

        x = round((1*w - 2*a)/4)
        y = round((3*h - 2*a)/4)
        try:
            matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 3"
            print e

        x = round((3*w - 2*a)/4)
        y = round((3*h - 2*a)/4)
        try:
            matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 4"
            print e

    if flag == 'interactive':
        for i in range(5):
            x, y, w, h = common.get_roi_coordinates(rpc1, prev1)
            try:
                matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, w, h)
                out = np.vstack((out, matches))
            except Exception as e:
                print "no matches in the selected roi"
                print e

    if flag == 'load':
        im = os.path.dirname(im1)
        fname = os.path.join(im, 'pointing_correction_rois.txt')
        rois = np.loadtxt(fname)
        for i in xrange(len(rois)):
            x, y, w, h = rois[i, :]
            try:
                matches = filtered_sift_matches_roi(im1, im2, rpc1, rpc2, x, y, w, h)
                out = np.vstack((out, matches))
            except Exception as e:
                print "no matches in the selected roi"
                print e


    # return the full list of matches, only if there are enough
    if len(out) < 7:
        raise Exception("not enough matches")
    else:
        return out[1:, :]

def query_rois_save_to_file(im, n=5):
    """
    Save coordinates of the ROIs to use for pointing correction of a dataset.

    Args:
        im: path to the folder containing the Pleiades images. Be careful, the
            path should not end by a '/'
        n: number of ROIs to query

    Returns:
        nothing. It saves the (x, y, w, h) of each ROI in a txt file named
        pointing_correction_rois.txt, one per line
    """
    rpc  = '%s/../rpc/%s/rpc01.xml' % (os.path.dirname(im), os.path.basename(im))
    prev = os.path.join(im, 'prev01.jpg')
    filename = os.path.join(im, 'pointing_correction_rois.txt')

    out = np.zeros((n, 4))
    for i in range(n):
        x, y, w, h = common.get_roi_coordinates(rpc, prev)
        out[i, :] = x, y, w, h

    np.savetxt(filename, out, '%d')

def query_rois_all_datasets(data):
    """
    Run the query_rois_save_to_file function on all the datasets found in data.

    Args:
        data: path to the folder containing all the Pleiades datasets

    Returns:
        nothing. For each dataset, it saves the (x, y, w, h) of each ROI in a
        txt file named dataset/pointing_correction_rois.txt, one per line
    """
    for f in os.listdir(data):
        query_rois_save_to_file(os.path.join(data, f))



def optimize_pair(im1, im2, rpc1, rpc2, prev1=None, matches=None):
    """
    Runs the pointing correction on a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        prev1 (optional): path to the jpg preview image of im1 (used in case of
            interactive mode in the matches computation)
        matches (optional): numpy 4xN array containing a list of matches
            between the two images. If it is not provided the matches are
            computed.

    Returns:
        a 3x3 matrix representing the planar transformation to apply to im2 in
        order to correct the pointing error.
    """

    if matches is None:
        matches = filtered_sift_matches_full_img(im1, im2, rpc1, rpc2,
            'load', prev1)

    from scipy.optimize import fmin_bfgs
    print "running optimization using %d matches" % len(matches)
    v0 = np.zeros((1, 4))
    v = fmin_bfgs(cost_function, v0, args=(rpc1, rpc2, matches), maxiter=30,
        retall=True)[0]

    return euclidean_transform_matrix(v)


def optimize_pair_all_datasets(data):
    """
    Run the optimize_pair function on all the datasets found in data.

    Args:
        data: path to the folder containing all the Pleiades datasets

    Returns:
        nothing. For each dataset, it saves the computed correction matrix in a
        txt file named dataset/pointing_correction_matrix.txt.
    """
    for f in os.listdir(data):
        dataset = os.path.join(data, f)
        im1  = os.path.join(dataset, 'im01.tif')
        im2  = os.path.join(dataset, 'im02.tif')
        rpc1 = '%s/../rpc/%s/rpc01.xml' % (data, f)
        rpc2 = '%s/../rpc/%s/rpc02.xml' % (data, f)
        rpc1 = rpc_model.RPCModel(rpc1)
        rpc2 = rpc_model.RPCModel(rpc2)

        out = os.path.join(dataset, 'pointing_correction_matrix.txt')
        A = optimize_pair(im1, im2, rpc1, rpc2)
        np.savetxt(out, A)
