# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


from __future__ import print_function
import os
import numpy as np

from s2plib import sift
from s2plib import rectification
from s2plib import rpc_utils
from s2plib import rpc_model
from s2plib import common
from s2plib import estimation
from s2plib import evaluation
from s2plib.config import cfg


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
    print('%d sift matches' % len(matches))

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
    print("max, mean, min pointing error, from %d points:" % (len(matches)))
    print(np.max(e), np.mean(e), np.min(e))

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
    print('%d sift matches' % len(matches))

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
        print('warning: cost_function is going too far')
        print(v)

    if (len(v) > 3):
        if (np.abs(v[3]) > 20000):
            print('warning: cost_function is going too far')
            print(v)

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
        print('warning: cost_function is going too far')
        print(v)

    x, y, w, h = common.bounding_box2D(matches[:, 0:2])
    matches_rpc = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, 5)
    F = estimation.fundamental_matrix(matches_rpc)

    # transform the coordinates of points in the second image according to
    # matrix A, built from vector v
    A = euclidean_transform_matrix(v)
    p2 = common.points_apply_homography(A, matches[:, 2:4])

    return evaluation.fundamental_matrix_L1(F, np.hstack([matches[:, 0:2], p2]))


def print_params(v):
    """
    Display the pointing correction parameters.

    Args:
        v: 1D numpy array of length 4

    This function is called by the fmin_bfgs optimization function at each
    iteration, to display the current values of the parameters.
    """
    print('rotation: %.3e, translation: (%.3e, %.3e), horizontal shear: %.3e' % v)


def optimize_pair(im1, im2, rpc1, rpc2, matches):
    """
    Runs the pointing correction on a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        matches: numpy 4xN array containing a list of matches between the two
            images.

    Returns:
        a 3x3 matrix representing the planar transformation to apply to im2 in
        order to correct the pointing error.
    """
    # Don't use too many matches to keep the evaluation time of 'cost_function'
    # reasonable
    # if len(matches) > 1000:
    #     ind = np.linspace(0, len(matches), 1000, False)
    #     matches = matches[ind.astype(int), :]

    from scipy.optimize import fmin_l_bfgs_b
    print("running optimization using %d matches" % len(matches))
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
    print('theta: %f' % v[0])
    print('tx: %f' % v[1])
    print('ty: %f' % v[2])
    print('shear: %f' % v[3])
    print('min cost: %f' % min_val)
    print(debug)
    return euclidean_transform_matrix(v)


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
        print("pointing_accuracy.error_vectors: invalid 'ind' argument")

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


def compute_correction(img1, rpc1, img2, rpc2, x, y, w, h):
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

    Returns:
        a 3x3 matrix representing the planar transformation to apply to img2 in
        order to correct the pointing error, and the list of sift matches used
        to compute this correction.
    """
    # read rpcs
    r1 = rpc_model.RPCModel(rpc1)
    r2 = rpc_model.RPCModel(rpc2)

    m = sift.matches_on_rpc_roi(img1, img2, r1, r2, x, y, w, h)

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
