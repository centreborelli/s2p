# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from __future__ import print_function
import numpy as np

from s2p import common


def fundamental_matrix_cameras(P1, P2):
    """
    Computes the fundamental matrix given the matrices of two cameras.

    Arguments:
        P1, P2: 2D arrays of size 3x4 containing the camera matrices

    Returns:
        the computed fundamental matrix, given by the formula 17.3 (p. 412) in
        Hartley & Zisserman book (2nd ed.).
    """
    X0 = P1[[1, 2], :]
    X1 = P1[[2, 0], :]
    X2 = P1[[0, 1], :]
    Y0 = P2[[1, 2], :]
    Y1 = P2[[2, 0], :]
    Y2 = P2[[0, 1], :]

    F = np.zeros((3, 3))
    F[0, 0] = np.linalg.det(np.vstack([X0, Y0]))
    F[0, 1] = np.linalg.det(np.vstack([X1, Y0]))
    F[0, 2] = np.linalg.det(np.vstack([X2, Y0]))
    F[1, 0] = np.linalg.det(np.vstack([X0, Y1]))
    F[1, 1] = np.linalg.det(np.vstack([X1, Y1]))
    F[1, 2] = np.linalg.det(np.vstack([X2, Y1]))
    F[2, 0] = np.linalg.det(np.vstack([X0, Y2]))
    F[2, 1] = np.linalg.det(np.vstack([X1, Y2]))
    F[2, 2] = np.linalg.det(np.vstack([X2, Y2]))

    return F


def get_angle_from_cos_and_sin(c, s):
    """
    Computes x in ]-pi, pi] such that cos(x) = c and sin(x) = s.
    """
    if s >= 0:
        return np.arccos(c)
    else:
        return -np.arccos(c)


def rectifying_similarities_from_affine_fundamental_matrix(F, debug=False):
    """
    Computes two similarities from an affine fundamental matrix.

    Args:
        F: 3x3 numpy array representing the input fundamental matrix
        debug (optional, default is False): boolean flag to activate verbose
            mode

    Returns:
        S, S': two similarities such that, when used to resample the two images
            related by the fundamental matrix, the resampled images are
            stereo-rectified.
    """
    # check that the input matrix is an affine fundamental matrix
    assert(np.shape(F) == (3, 3))
    assert(np.linalg.matrix_rank(F) == 2)
    assert(F[0, 0] == 0)
    assert(F[0, 1] == 0)
    assert(F[1, 0] == 0)
    assert(F[1, 1] == 0)

    # notations
    a = F[2, 0]
    b = F[2, 1]
    c = F[0, 2]
    d = F[1, 2]
    e = F[2, 2]

    # rotations
    r = np.sqrt(a*a + b*b)
    s = np.sqrt(c*c + d*d)
    R1 = (1.0 / r) * np.array([[b, -a], [a, b]])
    R2 = (1.0 / s) * np.array([[-d, c], [-c, -d]])

    # zoom and translation
    z = np.sqrt(r / s)
    t = 0.5 * e / np.sqrt(r * s)

    if debug:
        theta_1 = get_angle_from_cos_and_sin(b, a)
        print("reference image:")
        print("\trotation: %f deg" % np.rad2deg(theta_1))
        print("\tzoom: %f" % z)
        print("\tvertical translation: %f" % t)
        print()
        theta_2 = get_angle_from_cos_and_sin(-d, -c)
        print("secondary image:")
        print("\trotation: %f deg" % np.rad2deg(theta_2))
        print("\tzoom: %f" % (1.0 / z))
        print("\tvertical translation: %f" % -t)

    # output similarities
    S1 = np.zeros((3, 3))
    S1[0:2, 0:2] = z * R1
    S1[1, 2] = t
    S1[2, 2] = 1

    S2 = np.zeros((3, 3))
    S2[0:2, 0:2] = (1.0 / z) * R2
    S2[1, 2] = -t
    S2[2, 2] = 1

    return S1, S2


def affine_fundamental_matrix(matches):
    """
    Estimates the affine fundamental matrix given a set of point correspondences
    between two images.

    Arguments:
        matches: 2D array of size Nx4 containing a list of pairs of matching
            points. Each line is of the form x1, y1, x2, y2, where (x1, y1) is
            the point in the first view while (x2, y2) is the matching point in
            the second view.

    Returns:
        the estimated affine fundamental matrix, given by the Gold Standard
        algorithm, as described in Hartley & Zisserman book (see chap. 14).
    """
    # revert the order of points to fit H&Z convention (see algo 14.1)
    X = matches[:, [2, 3, 0, 1]]

    # compute the centroid
    N = len(X)
    XX = np.sum(X, axis=0) / N

    # compute the Nx4 matrix A
    A = X - np.tile(XX, (N, 1))

    # the solution is obtained as the singular vector corresponding to the
    # smallest singular value of matrix A. See Hartley and Zissermann for
    # details.
    # It is the last line of matrix V (because np.linalg.svd returns V^T)
    U, S, V = np.linalg.svd(A)
    N = V[-1, :]

    # extract values and build F
    F = np.zeros((3, 3))
    F[0, 2] = N[0]
    F[1, 2] = N[1]
    F[2, 0] = N[2]
    F[2, 1] = N[3]
    F[2, 2] = -np.dot(N, XX)

    return F


def affine_transformation(x, xx):
    """
    Estimates an affine homography from a list of correspondences

    Args:
        x:  Nx2 numpy array, containing a list of points
        xx: Nx2 numpy array, containing the list of corresponding points

    Returns:
        3x3 numpy array, representing (in homogeneous coordinates) an affine
        homography that maps the points of x onto the points of xx.

    This function implements the Gold-Standard algorithm for estimating an
    affine homography, described in Hartley & Zisserman page 130 (second
    edition).
    """
    # check that there are at least 3 points
    if len(x) < 3:
        print("ERROR: estimation.affine_transformation needs 3 correspondences")
        return np.eye(3)

    # translate the input points so that the centroid is at the origin.
    t = -np.mean(x,  0)
    tt = -np.mean(xx, 0)
    x = x + t
    xx = xx + tt

    # compute the Nx4 matrix A
    A = np.hstack((x, xx))

    # two singular vectors corresponding to the two largest singular values of
    # matrix A. See Hartley and Zissermann for details.  These are the first
    # two lines of matrix V (because np.linalg.svd returns V^T)
    U, S, V = np.linalg.svd(A)
    v1 = V[0, :]
    v2 = V[1, :]

    # compute blocks B and C, then H
    tmp = np.vstack((v1, v2)).T
    assert(np.shape(tmp) == (4, 2))
    B = tmp[0:2, :]
    C = tmp[2:4, :]
    H = np.dot(C, np.linalg.inv(B))

    # return A
    A = np.eye(3)
    A[0:2, 0:2] = H
    A[0:2, 2] = np.dot(H, t) - tt
    return A


def translation(x, xx):
    """
    Estimates a planar translation from a list of correspondences

    Args:
        x:  Nx2 numpy array, containing a list of points
        xx: Nx2 numpy array, containing the list of corresponding points

    Returns:
        3x3 numpy array, representing (in homogeneous coordinates) a planar
        translation that maps the points of x onto the points of xx.
    """
    # compute the mean of the displacement vectors
    t = np.mean(xx - x, 0)

    # return A
    A = np.eye(3)
    A[0, 2] = t[0]
    A[1, 2] = t[1]
    return A
