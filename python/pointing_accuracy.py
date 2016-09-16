# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import os
import numpy as np

import sift
import rectification
import rpc_utils
import rpc_model
import common
import estimation
import evaluation
from config import cfg


def evaluation_iterative(im1, im2, rpc1, rpc2, x, y, w, h, A=None,
                         matches=None):
    """
    Measures the maximal pointing error on a Pleiades' pair of images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers defining the rectangular ROI in the first
            image.  (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2.
        matches (optional): Nx4 numpy array containing a list of matches to use
            to compute the pointing error

    Returns:
        the mean pointing error, in the direction orthogonal to the epipolar
        lines. This error is measured in pixels.
    """
    if not matches:
        matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x, y, w, h)
    p1 = matches[:, 0:2]
    p2 = matches[:, 2:4]
    print '%d sift matches' % len(matches)

    # apply pointing correction matrix, if available
    if A is not None:
        p2 = common.points_apply_homography(A, p2)

    # compute the pointing error for each match
    x1 = p1[:, 0]
    y1 = p1[:, 1]
    x2 = p2[:, 0]
    y2 = p2[:, 1]
    e = rpc_utils.compute_height(rpc1, rpc2, x1, y1, x2, y2)[1]
    # matches = matches[e < 0.1, :]
    # visualisation.plot_matches_pleiades(im1, im2, matches)
    print "max, mean, min pointing error, from %d points:" % (len(matches))
    print np.max(e), np.mean(e), np.min(e)

    # return the mean error
    return np.mean(np.abs(e))


def evaluation_from_estimated_F(im1, im2, rpc1, rpc2, x, y, w, h, A=None,
        matches=None):
    """
    Measures the pointing error on a Pleiades' pair of images, affine approx.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers defining the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2.
        matches (optional): Nx4 numpy array containing a list of matches to use
            to compute the pointing error

    Returns:
        the mean pointing error, in the direction orthogonal to the epipolar
        lines. This error is measured in pixels, and computed from an
        approximated fundamental matrix.
    """
    if not matches:
        matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x, y, w, h)
    p1 = matches[:, 0:2]
    p2 = matches[:, 2:4]
    print '%d sift matches' % len(matches)

    # apply pointing correction matrix, if available
    if A is not None:
        p2 = common.points_apply_homography(A, p2)

    # estimate the fundamental matrix between the two views
    rpc_matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, 5)
    F = estimation.affine_fundamental_matrix(rpc_matches)

    # compute the mean displacement from epipolar lines
    d_sum = 0
    for i in range(len(p1)):
        x  = np.array([p1[i, 0], p1[i, 1], 1])
        xx = np.array([p2[i, 0], p2[i, 1], 1])
        ll  = F.dot(x)
        #d = np.sign(xx.dot(ll)) * evaluation.distance_point_to_line(xx, ll)
        d = evaluation.distance_point_to_line(xx, ll)
        d_sum += d
    return d_sum/len(p1)


def euclidean_transform_matrix(v):
    """
    Build the matrix of an euclidean transform in homogeneous coordinates.

    Arguments:
        v: numpy 1D array, of length 3 or 4, containing the parameters of the
            planar transform. These parameters are:
            v[0]: angle of the rotation, in radians
            v[1]: x coordinate of the vector of the translation
            v[2]: y coordinate of the vector of the translation
            v[3]: horizontal shear parameter

    Returns:
        A numpy 3x3 matrix, representing the euclidean transform (rotation
        followed by translation) in homogeneous coordinates.
    """
    # variables scaling
    v[0] = v[0]/1000000
    v[3] = v[3]/1000000

    # centering of coordinates with respect to the center of the big pleiades
    # image (40000x40000)
    C = np.eye(3)
    C[0, 2] = -20000
    C[1, 2] = -20000

    # matrix construction
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
    #return np.linalg.inv(C).dot(T).dot(R).dot(S).dot(C)
    return np.dot(np.dot(np.dot(np.dot(np.linalg.inv(C), T), R), S), C)


def cost_function(v, *args):
    """
    Objective function to minimize in order to correct the pointing error.

    Arguments:
        v: vector of size 3 or 4, containing the parameters of the euclidean
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
    rpc1, rpc2, matches = args[0], args[1], args[2]
    if len(args) == 4:
        alpha = args[3]
    else:
        alpha = 0.01

    # verify that parameters are in the bounding box
    if (np.abs(v[0]) > 200*np.pi or
        np.abs(v[1]) > 10000 or
        np.abs(v[2]) > 10000):
        print 'warning: cost_function is going too far'
        print v

    if (len(v) > 3):
        if (np.abs(v[3]) > 20000):
            print 'warning: cost_function is going too far'
            print v

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


def cost_function_linear(v, rpc1, rpc2, matches):
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
    print_params(v)

    # verify that parameters are in the bounding box
    if (np.abs(v[0]) > 200*np.pi or
        np.abs(v[1]) > 10000 or
        np.abs(v[2]) > 10000 or
        np.abs(v[3]) > 20000):
        print 'warning: cost_function is going too far'
        print v

    x, y, w, h = common.bounding_box2D(matches[:, 0:2])
    matches_rpc = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, 5)
    F = estimation.fundamental_matrix(matches_rpc)

    # transform the coordinates of points in the second image according to
    # matrix A, built from vector v
    A = euclidean_transform_matrix(v)
    p2 = common.points_apply_homography(A, matches[:, 2:4])

    return evaluation.fundamental_matrix_L1(F, np.hstack([matches[:, 0:2], p2]))


def filtered_sift_matches_full_img(im1, im2, rpc1, rpc2, flag='automatic',
        prev1=None, a=1000, x=None, y=None, w=None, h=None, outfile=None, model='fundamental'):
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
        a: length of the squared tiles used to extract sift points, in the case
            of automatic mode
        x, y, w, h (optional): use a big ROI and extract the five tiles from
            there instead of from the full image.
        outfile (optional): path to a txt where to save the list of matches.
        model (optional, default is 'fundamental'): model imposed by RANSAC
            when searching the set of inliers

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the big images.
            If no sift matches are found, then an exception is raised.

    The keypoints are extracted from five zones in the first image. One in the
    center, and four in the corners. The matches are consistent with an
    epipolar model in each of the zones.
    """
    # if no big ROI is defined, use the full image
    if x is None:
        x = 0
        y = 0
        w = rpc1.lastCol
        h = rpc1.lastRow

    # initialize output array
    out = np.zeros(shape = (1, 4))

    if flag == 'automatic':
        # if roi size is too big wrt image size, take a smaller roi size
        if (min(h, w) < 4*a):
            a = round(min(h, w)/4)

        # central zone
        x0 = round((w-a)/2) + x
        y0 = round((h-a)/2) + y
        try:
            matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x0, y0, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the central zone"
            print e

        # corner zones
        x0 = round((1*w - 2*a)/4) + x
        y0 = round((1*h - 2*a)/4) + y
        try:
            matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x0, y0, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 1"
            print e

        x0 = round((3*w - 2*a)/4) + x
        y0 = round((1*h - 2*a)/4) + y
        try:
            matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x0, y0, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 2"
            print e

        x0 = round((1*w - 2*a)/4) + x
        y0 = round((3*h - 2*a)/4) + y
        try:
            matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x0, y0, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 3"
            print e

        x0 = round((3*w - 2*a)/4) + x
        y0 = round((3*h - 2*a)/4) + y
        try:
            matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x0, y0, a, a)
            out = np.vstack((out, matches))
        except Exception as e:
            print "no matches in the corner 4"
            print e

    if flag == 'interactive':
        for i in range(5):
            x, y, w, h = common.get_roi_coordinates(rpc1, prev1)
            try:
                matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x, y, w, h)
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
                matches = sift.matches_on_rpc_roi(im1, im2, rpc1, rpc2, x, y, w, h)
                out = np.vstack((out, matches))
            except Exception as e:
                print "no matches in the selected roi"
                print e


    # save and return the full list of matches, only if there are enough
    if len(out) < 7:
        raise Exception("not enough matches")
    else:
        if outfile is not None:
            np.savetxt(outfile, out[1:, :])
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


def save_sift_matches_all_datasets(data):
    """
    Run the filtered_sift_matches_full_img function on all the datasets.

    Args:
        data: path to the folder containing all the Pleiades datasets

    Returns:
        nothing. For each dataset, it saves the list of sift matches computed
        between images 1 and 2, according to the ROIs defined in the
        corresponding file.
    """
    for f in os.listdir(data):
        dataset = os.path.join(data, f)
        im1  = os.path.join(dataset, 'im01.tif')
        im2  = os.path.join(dataset, 'im02.tif')
        rpc1 = '%s/../rpc/%s/rpc01.xml' % (data, f)
        rpc2 = '%s/../rpc/%s/rpc02.xml' % (data, f)
        rpc1 = rpc_model.RPCModel(rpc1)
        rpc2 = rpc_model.RPCModel(rpc2)
        fname = os.path.join(dataset, 'sift_matches.txt')
        print fname
        m = filtered_sift_matches_full_img(im1, im2, rpc1, rpc2, 'load', None,
            1000, None, None, None, None, fname)


def print_params(v):
    """
    Display the pointing correction parameters.

    Args:
        v: 1D numpy array of length 4

    This function is called by the fmin_bfgs optimization function at each
    iteration, to display the current values of the parameters.
    """
    print 'rotation: %.3e, translation: (%.3e, %.3e), horizontal shear: %.3e' %(v[0],
        v[1], v[2], v[3])


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
        matches = filtered_sift_matches_full_img(im1, im2, rpc1, rpc2, 'load',
                                                 prev1)

    # Don't use too many matches to keep the evaluation time of 'cost_function'
    # reasonable
    # if len(matches) > 1000:
    #     ind = np.linspace(0, len(matches), 1000, False)
    #     matches = matches[ind.astype(int), :]

    from scipy.optimize import fmin_l_bfgs_b
    print "running optimization using %d matches" % len(matches)
    v0 = np.zeros(4)
    v, min_val, debug = fmin_l_bfgs_b(
            cost_function,
            v0,
            args=(rpc1, rpc2, matches),
            approx_grad=True,
            factr=1,
            bounds=[(-150, 150), (-100, 100), (-100, 100), (-200000, 200000)])
            #maxiter=50,
            #callback=print_params,
            #iprint=0,
            #disp=0)

    # default values are:
    # fmin_l_bfgs_b(func, x0, fprime=None, args=(), approx_grad=0, bounds=None,
    # m=10, factr=10000000.0, pgtol=1e-05, epsilon=1e-08, iprint=-1,
    # maxfun=15000, disp=None)
    print 'theta: %f' % v[0]
    print 'tx: %f' % v[1]
    print 'ty: %f' % v[2]
    print 'shear: %f' % v[3]
    print 'min cost: %f' % min_val
    print debug
    return euclidean_transform_matrix(v)


def optimize_pair_all_datasets(data):
    """
    Run the optimize_pair function on all the datasets found in data.

    Args:
        data: path to the folder containing all the Pleiades datasets

    Returns:
        nothing. For each dataset, it saves the computed correction matrix in a
        txt file named dataset/pointing_correction_matrix.txt.

    It uses precomputed matches, stored in dataset/sift_matches.txt. These
    matches are usually computed with the save_sift_matches_all_datasets
    function.
    """
    for f in os.listdir(data):
        dataset = os.path.join(data, f)
        im1  = os.path.join(dataset, 'im01.tif')
        im2  = os.path.join(dataset, 'im02.tif')
        rpc1 = '%s/../rpc/%s/rpc01.xml' % (data, f)
        rpc2 = '%s/../rpc/%s/rpc02.xml' % (data, f)
        rpc1 = rpc_model.RPCModel(rpc1)
        rpc2 = rpc_model.RPCModel(rpc2)
        matches_file = os.path.join(dataset, 'sift_matches.txt')
        try:
            with open(matches_file, 'r'):
                matches = np.loadtxt(matches_file)
                out = os.path.join(dataset, 'pointing_correction_matrix.txt')
                print f
                A = optimize_pair(im1, im2, rpc1, rpc2, None, matches)
                np.savetxt(out, A)
        except Exception as e:
                print e


def error_vectors(m, F, ind='ref'):
    """
    Computes the error vectors for a list of matches and a given fundamental
    matrix.

    Args:
        m: Nx4 numpy array containing a list of matches, one per line. Each
           match is given by (x, y, x', y') where (x, y) is a point of the
           reference view and (x', y') is the corresponding point in the
           secondary view.
        F: fundamental matrix between the two views
        ind (optional, default is 'ref'): index of the image on which the lines
            are plotted to compute the error vectors. Must be either 'ref' or
            'sec' (reference or secondary image)

    Returns:
        Nx2 numpy array containing a list of planar vectors. Each vector is
        obtained as the difference between x' and its projection on the
        epipolar line Fx.
    """
    # divide keypoints in two lists: x (first image) and xx (second image), and
    # convert them to homogeneous coordinates
    N  = len(m)
    x  = np.ones((N, 3))
    xx = np.ones((N, 3))
    x [:, 0:2] = m[:, 0:2]
    xx[:, 0:2] = m[:, 2:4]

    # epipolar lines: 2D array of size Nx3, one epipolar line per row
    if ind == 'sec':
        l = np.dot(x, F.T)
    elif ind == 'ref':
        l = np.dot(xx, F)
    else:
        print "pointing_accuracy.error_vectors: invalid 'ind' argument"

    # compute the error vectors (going from the projection of x or xx on l to x
    # or xx)
    if ind == 'sec':
        n = np.multiply(xx[:, 0], l[:, 0]) + np.multiply(xx[:, 1], l[:, 1]) + l[:, 2]
    else:
        n = np.multiply(x[:, 0], l[:, 0]) + np.multiply(x[:, 1], l[:, 1]) + l[:, 2]
    d = np.square(l[:, 0]) + np.square(l[:, 1])
    a = np.divide(n, d)
    return np.vstack((np.multiply(a, l[:, 0]), np.multiply(a, l[:, 1]))).T


def local_translation(r1, r2, x, y, w, h, m):
    """
    Estimates the optimal translation to minimise the relative pointing error
    on a given tile.

    Args:
        r1, r2: two instances of the rpc_model.RPCModel class
        x, y, w, h: region of interest in the reference image (r1)
        m: Nx4 numpy array containing a list of matches, one per line. Each
            match is given by (p1, p2, q1, q2) where (p1, p2) is a point of the
            reference view and (q1, q2) is the corresponding point in the
            secondary view.

    Returns:
        3x3 numpy array containing the homogeneous representation of the
        optimal planar translation, to be applied to the secondary image in
        order to correct the pointing error.
    """
    # estimate the affine fundamental matrix between the two views
    n = cfg['n_gcp_per_axis']
    rpc_matches = rpc_utils.matches_from_rpc(r1, r2, x, y, w, h, n)
    F = estimation.affine_fundamental_matrix(rpc_matches)

    # compute the error vectors
    e = error_vectors(m, F, 'sec')

    # compute the median: as the vectors are collinear (because F is affine)
    # computing the median of each component independently is correct
    N = len(e)
    out_x = np.sort(e[:, 0])[int(N/2)]
    out_y = np.sort(e[:, 1])[int(N/2)]

    # the correction to be applied to the second view is the opposite
    A = np.array([[1, 0, -out_x],
                  [0, 1, -out_y],
                  [0, 0, 1]])
    return A


def compute_correction(img1, rpc1, img2, rpc2, x, y, w, h,
                       filter_matches='fundamental'):
    """
    Computes pointing correction matrix for specific ROI

    Args:
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want. If bigger than
            1 Mpix, only five crops will be used to compute sift matches.
        filter_matches (optional, default is 'fundamental'): model imposed by
            RANSAC when searching the set of inliers

    Returns:
        a 3x3 matrix representing the planar transformation to apply to img2 in
        order to correct the pointing error, and the list of sift matches used
        to compute this correction.
    """
    # read rpcs
    r1 = rpc_model.RPCModel(rpc1)
    r2 = rpc_model.RPCModel(rpc2)

    m = sift.matches_on_rpc_roi(img1, img2, r1, r2, x, y, w, h)

    # A = optimize_pair(img1, img2, r1, r2, None, m)
    if m is not None:
        A = local_translation(r1, r2, x, y, w, h, m)
    else:
        A = None

    return A, m


def from_next_tiles(tiles, ntx, nty, col, row):
    """
    Computes the pointing correction of a specific tile using the pointing
    corrections of its neighbors. We assume that all pointing corrections are
    translations.

    Args:
        tiles: list of paths to folders associated to each tile
        ntx: number of tiles in horizontal direction
        nty: number of tiles in vertical direction
        col: horizontal index (from 1 to ntx) of the current tile
        row: horizontal index (from 1 to nty) of the current tile

    Returns:
        the estimated pointing correction for the specified tile
    """
    # TODO
    return None


def global_from_local(tiles):
    """
    Computes the pointing correction of a full roi using local corrections on
    tiles.

    Args:
        tiles: list of paths to folders associated to each tile

    Returns:
        the estimated pointing correction for the specified tile

    In each folder we expect to find the files pointing.txt and center.txt. The
    file pointing.txt contains the local correction (a projective transform
    given in homogeneous coordinates), and the file center.txt contains the
    coordinates of the mean of the keypoints of the secondary image.
    """
    # lists of matching points
    x  = []
    xx = []

    # loop over all the tiles
    for f in tiles:
        center = os.path.join(f, 'center_keypts_sec.txt')
        pointing = os.path.join(f, 'pointing.txt')
        if os.path.isfile(center) and os.path.isfile(pointing):
            A = np.loadtxt(pointing)
            p = np.loadtxt(center)
            if A.shape == (3, 3) and p.shape == (2,):
                q = np.dot(A, np.array([p[0], p[1], 1]))
                x.append(p)
                xx.append(q[0:2])

    if not x:
        return np.eye(3)
    elif len(x) == 1:
        return A
    elif len(x) == 2:
        #TODO: replace translation with similarity
        return estimation.translation(np.array(x), np.array(xx))
    else:
        # estimate an affine transformation transforming x in xx
        return estimation.affine_transformation(np.array(x), np.array(xx))
