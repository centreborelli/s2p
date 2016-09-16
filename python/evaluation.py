# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from __future__ import print_function
import numpy as np

def distance_point_to_line(x, l):
    """
    Computes the distance between a point and a line expressed in homogeneous
    coordinates.

    Args:
        x: 3-vector, containing the homogeneous coordinates of a point
        l: 3-vector, containing the homogeneous coordinates of a line

    Returns:
        the distance between x and l
        In the case where x is an ideal point or l is the line at infinity, the
        returned distance is infinity
    """
    if (np.abs(x[2]) < np.finfo(float).eps):
        # x is an ideal point, distance is +infty
        return np.finfo(float).max
    if (np.hypot(l[0], l[1]) < np.finfo(float).eps):
        # l is the line at infinity, distance is +infty
        return np.finfo(float).max

    num = np.abs(np.dot(x, l))
    den = np.hypot(l[0], l[1]) * np.abs(x[2])
    return num/den


def max_dist_points_to_lines(x, l):
    """
    Args:
        x: 2D array of size Nx2 containing a list of points in the plane
        l: 2D array of size Nx3 containing a list of lines given by their
            homogeneous coordinates

    Returns:
        the highest pointwise distance between the points and the lines.

    """
    a = np.multiply(x[:, 0], l[:, 0]) + np.multiply(x[:, 1], l[:, 1]) + l[:, 2]
    b = np.hypot(l[:, 0], l[:, 1])
    d = np.abs(np.divide(a, b))
    return np.max(d)


def fundamental_matrix_fast(F, matches):
    """
    Evaluates the precision of a fundamental matrix against a set of point
    correspondences.

    Arguments:
        F: fundamental matrix
        matches: 2D array of size Nx4 containing a list of pairs of matching
            points. Each line is of the form x1, y1, x2, y2, where (x1, y1) is
            the point in the first view while (x2, y2) is the matching point in
            the second view.

    Returns:
        the highest symmetric residual error, ie the maximum over all the
        matches of the following quantity:
                    max( d(x_i, F^Tx'_i), d(x'_i, Fx_i) ),
        where we use the notations of Hartley and Zisserman
    """
    N  = len(matches)
    x  = np.ones((N, 3))
    xx = np.ones((N, 3))

    x[:, 0:2]  = matches[:, 0:2]
    xx[:, 0:2] = matches[:, 2:4]
    l =  np.dot(xx, F)
    ll = np.dot(x, F.T)
    d1 = max_dist_points_to_lines(matches[:, 0:2], l)
    d2 = max_dist_points_to_lines(matches[:, 2:4], ll)
    return max(d1, d2)



def fundamental_matrix(F, matches):
    """
    Evaluates the precision of a fundamental matrix against a set of point
    correspondences.

    Arguments:
        F: fundamental matrix
        matches: 2D array of size Nx4 containing a list of pairs of matching
            points. Each line is of the form x1, y1, x2, y2, where (x1, y1) is
            the point in the first view while (x2, y2) is the matching point in
            the second view.

    Returns:
        the highest symmetric residual error, ie the maximum over all the
        matches of the following quantity:
                    max( d(x_i, F^Tx'_i), d(x'_i, Fx_i) ),
        where we use the notations of Hartley and Zisserman
    """
    d_max = 0
    for i in range(len(matches)):
        x  = np.array([matches[i, 0], matches[i, 1], 1])
        xx = np.array([matches[i, 2], matches[i, 3], 1])
        l = np.dot(F.T, xx)
        ll  = np.dot(F, x)
        d1 = distance_point_to_line(x, l)
        d2 = distance_point_to_line(xx, ll)
        d = max(d1, d2)
#        print(d) # for debug only
        if (d > d_max):
            d_max = d
    return d_max

def fundamental_matrix_L1(F, matches):
    """
    Evaluates the precision of a fundamental matrix against a set of point
    correspondences.

    Arguments:
        F: fundamental matrix
        matches: 2D array of size Nx4 containing a list of pairs of matching
            points. Each line is of the form x1, y1, x2, y2, where (x1, y1) is
            the point in the first view while (x2, y2) is the matching point in
            the second view.

    Returns:
        the sum of symmetric residual error, ie the sum over all the
        matches of the following quantity:
                    max( d(x_i, F^Tx'_i), d(x'_i, Fx_i) ),
        where we use the notations of Hartley and Zisserman
    """
    d_sum = 0
    for i in range(len(matches)):
        x  = np.array([matches[i, 0], matches[i, 1], 1])
        xx = np.array([matches[i, 2], matches[i, 3], 1])
        l = np.dot(F.T, xx)
        ll  = np.dot(F, x)
        d1 = distance_point_to_line(x, l)
        d2 = distance_point_to_line(xx, ll)
        d = max(d1, d2)
        d_sum += d
    return d_sum

def camera_matrix(P, X, x):
    """
    Evaluates the precision of a camera matrix against a set of world to image
    correspondences.

    Arguments:
        P: camera matrix
        X: 2D array of size Nx3 containing the coordinates of the 3-space
            points, one point per line
        x: 2D array of size Nx2 containing the pixel coordinates of the imaged
            points, one point per line
            These two arrays are supposed to have the same number of lines.

    Returns:
        the highest projection error, ie the maximum over all the
        matches of the distance d(x_i, PX_i), computed in the image plane,
        where we use the notations of Hartley and Zisserman.
    """
    d_max = 0
    for i in range(len(x)):
        p_world = np.append(X[i], 1)
        p_proj = np.dot(P, p_world)
        p_proj = p_proj/p_proj[2]
        d = np.linalg.norm(p_proj[0:2] - x[i])
        if (d > d_max):
            d_max = d
    return d_max


def compare_homogeneous(a, b):
    """
    Computes a distance between two homogeneous vectors or matrices.

    The two homogeneous quantities are normalized such that their L2 norm is 1,
    then the norm of their difference is computed.

    Args:
        a, b: two homogeneous quantities, represented as numpy nd-arrays. They
            must have the same shape

    Returns:
        the distance between normalized input data
    """
    if (a.shape != b.shape):
        print("compare_homogeneous: inputs must have the same shape")
        return

    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)
    d = a - b
    return np.linalg.norm(d)
