# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


from __future__ import print_function
import os
import numpy as np

from s2p import sift
from s2p import rpc_utils
from s2p import rpc_model
from s2p import estimation
from s2p.config import cfg


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


def compute_correction(img1, img2, x, y, w, h):
    """
    Computes pointing correction matrix for specific ROI

    Args:
        img1: path to the reference image.
        img2: path to the secondary image.
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
    r1 = rpc_utils.rpc_from_geotiff(img1)
    r2 = rpc_utils.rpc_from_geotiff(img2)

    m = sift.matches_on_rpc_roi(img1, img2, r1, r2, x, y, w, h)

    if m is not None:
        A = local_translation(r1, r2, x, y, w, h, m)
    else:
        A = None

    return A, m


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
