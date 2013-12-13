# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>

import numpy as np
import common

def normalize_2d_points(pts):
    """
    Translates and scales 2D points.

    The input points are translated and scaled such that the output points are
    centered at origin and the mean distance from the origin is sqrt(2). As
    shown in Hartley (1997), this normalization process typically improves the
    condition number of the linear systems used for solving homographies,
    fundamental matrices, etc.

    References:
        Richard Hartley, PAMI 1997
        Peter Kovesi, MATLAB functions for computer vision and image processing,

    Args:
        pts: 2D array of dimension Nx2 containing the coordinates of the input
            points, one point per line

    Returns:
        new_pts, T: coordinates of the transformed points, together with
            the similarity transform
    """
    # centroid
    cx = np.mean(pts[:, 0])
    cy = np.mean(pts[:, 1])

    # shift origin to centroid
    new_x = pts[:, 0] - cx
    new_y = pts[:, 1] - cy

    # scale such that the average distance from centroid is \sqrt{2}
    mean_dist = np.mean(np.sqrt(new_x**2 + new_y**2))
    s = np.sqrt(2)/mean_dist
    new_x = s*new_x
    new_y = s*new_y

    T = np.eye(3)
    T[0, 0] = s
    T[1, 1] = s         #              s     0   -s*cx
    T[0, 2] = -s*cx     # matrix T  =  0     s   -s*cy
    T[1, 2] = -s*cy     #              0     0     1

    return np.vstack([new_x, new_y]).T, T


def normalize_3d_points(pts):
    """
    Translates and scales 3D points.

    The input points are translated and scaled such that the output points are
    centered at origin and the mean distance from the origin is sqrt(3).

    Args:
        pts: 2D array of dimension Nx3 containing the coordinates of the input
            points, one point per line

    Returns:
        new_pts, U: coordinates of the transformed points, together with
            the similarity transform
    """
    # centroid
    cx = np.mean(pts[:, 0])
    cy = np.mean(pts[:, 1])
    cz = np.mean(pts[:, 2])

    # shift origin to centroid
    new_x = pts[:, 0] - cx
    new_y = pts[:, 1] - cy
    new_z = pts[:, 2] - cz

    # scale such that the average distance from centroid is \sqrt{3}
    mean_dist = np.mean(np.sqrt(new_x**2 + new_y**2 + new_z**2))
    s = np.sqrt(3)/mean_dist
    new_x = s*new_x
    new_y = s*new_y
    new_z = s*new_z

    U = np.eye(4)
    U[0, 0] = s
    U[1, 1] = s         #              s     0      0    -s*cx
    U[2, 2] = s         #              0     s      0    -s*cy
    U[0, 3] = -s*cx     # matrix U  =  0     0      s    -s*cz
    U[1, 3] = -s*cy     #              0     0      0      1
    U[2, 3] = -s*cz

    return np.vstack([new_x, new_y, new_z]).T, U


def camera_matrix(X, x):
    """
    Estimates the camera projection matrix from corresponding 3-space and image
    entities.

    Arguments:
        X: 2D array of size Nx3 containing the coordinates of the 3-space
            points, one point per line
        x: 2D array of size Nx2 containing the pixel coordinates of the imaged
            points, one point per line
            These two arrays are supposed to have the same number of lines.

    Returns:
        the estimated camera projection matrix, given by the Direct Linear
        Transformation algorithm, as described in Hartley & Zisserman book.
    """
    # normalize the input coordinates
    X, U = normalize_3d_points(X)
    x, T = normalize_2d_points(x)

    # make a linear system A*P = 0 from the correspondances, where P is made of
    # the 12 entries of the projection matrix (encoded in a vector P). This
    # system corresponds to the concatenation of correspondance constraints
    # (X_i --> x_i) which can be written as:
    # x_i x P*X_i = 0 (the vectorial product is 0)
    # and lead to 2 independent equations, for each correspondance. The system
    # is thus of size 2n x 12, where n is the number of correspondances. See
    # Zissermann, chapter 7, for more details.

    A = np.zeros((len(x)*2, 12))
    for i in xrange(len(x)):
        A[2*i+0, 4:8]  =       -1*np.array([X[i, 0], X[i, 1], X[i, 2], 1])
        A[2*i+0, 8:12] =  x[i, 1]*np.array([X[i, 0], X[i, 1], X[i, 2], 1])
        A[2*i+1, 0:4]  =          np.array([X[i, 0], X[i, 1], X[i, 2], 1])
        A[2*i+1, 8:12] = -x[i, 0]*np.array([X[i, 0], X[i, 1], X[i, 2], 1])

    # the vector P we are looking for minimizes the norm of A*P, and satisfies
    # the constraint \norm{P}=1 (to avoid the trivial solution P=0). This
    # solution is obtained as the singular vector corresponding to the smallest
    # singular value of matrix A. See Zissermann for details.
    # It is the last line of matrix V (because np.linalg.svd returns V^T)
    W, S, V = np.linalg.svd(A)
    P = V[-1, :].reshape((3, 4))

    # denormalize P
    # implements P = T^-1 * P * U
    P = np.dot(np.dot(np.linalg.inv(T), P), U)
    return P


def fundamental_matrix(matches):
    """
    Estimates the fundamental matrix given a set of point correspondences
    between two images.

    Arguments:
        matches: 2D array of size Nx4 containing a list of pairs of matching
            points. Each line is of the form x1, y1, x2, y2, where (x1, y1) is
            the point in the first view while (x2, y2) is the matching point in
            the second view.

    Returns:
        the estimated fundamental matrix, given by the normalized 8-points
        algorithm, as described in Hartley & Zisserman book.
    """
    # normalize the input coordinates
    pts1, T1 = normalize_2d_points(matches[:, 0:2])
    pts2, T2 = normalize_2d_points(matches[:, 2:4])

    # build the matrix, given by eq 11.3 (p279) in Hartley & Zisserman
    A = np.zeros((len(matches), 9))
    for i in xrange(len(matches)):
        A[i, 0:3] = pts2[i, 0] * np.array([pts1[i, 0], pts1[i, 1], 1])
        A[i, 3:6] = pts2[i, 1] * np.array([pts1[i, 0], pts1[i, 1], 1])
        A[i, 6:9] =              np.array([pts1[i, 0], pts1[i, 1], 1])

    # the vector F we are looking for minimizes the norm of A*F, and satisfies
    # the constraint \norm{F}=1 (to avoid the trivial solution F=0). This
    # solution is obtained as the singular vector corresponding to the smallest
    # singular value of matrix A. See Hartley and Zissermann for details.
    # It is the last line of matrix V (because np.linalg.svd returns V^T)
    U, S, V = np.linalg.svd(A)
    F = V[-1, :].reshape((3, 3))

    # constraint enforcement (F has to be rank 2)
    U, D, V = np.linalg.svd(F)
    D = np.diag([D[0], D[1], 0])
    F = np.dot(np.dot(U, D), V)

    # denormalize F
    F = np.dot(np.dot(T2.T, F), T1)
    return F


def fundamental_matrix_ransac(matches, precision=1.0):
    """
    Estimates the fundamental matrix given a set of point correspondences
    between two images, using ransac.

    Arguments:
        matches: 2D array of size Nx4 containing a list of pair of matching
            points. Each line is of the form x1, y1, x2, y2, where (x1, y1) is
            the point in the first view while (x2, y2) is the matching point in
            the second view.
        precision: optional parameter indicating the maximum error 
            allowed for counting the inliers 

    Returns:
        the estimated fundamental matrix

    The algorithm uses ransac as a search engine.
    """
    # write a file containing the list of correspondences. The
    # expected format is a raw text file with one match per line: x1 y1 x2 y2
    matchfile = common.tmpfile('.txt')
    np.savetxt(matchfile, matches)

    # call ransac binary, from Enric's imscript
    inliers = common.tmpfile('.txt')
    Ffile = common.tmpfile('.txt')
    common.run("""
        ransac fmn 1000 %f 7 %s < %s |
        grep parameters |
        awk \'{ print "[ " $3 " " $4 " " $5 " ; " $6 " " $7 " " $8 " ; " $9 " " $10 " " $11 " ] " }\' |
        tail -1 > %s
        """ % (precision, inliers, matchfile, Ffile) )
    common.matrix_write(Ffile, (common.matrix_read(Ffile, 3, 3)).transpose())
    return common.matrix_read(Ffile, 3, 3)


def loop_zhang(F, w, h):
    """
    Computes rectifying homographies from a fundamental matrix, with Loop-Zhang.

    Args:
        F: 3x3 numpy array containing the fundamental matrix
        w, h: images size. The two images are supposed to have same size

    Returns:
        The two rectifying homographies.

    The rectifying homographies are computed using the Pascal Monasse binary
    named rectify_mindistortion. It uses the Loop-Zhang algorithm.
    """
    Ffile = common.tmpfile('.txt')
    Haf = common.tmpfile('.txt')
    Hbf = common.tmpfile('.txt')
    common.matrix_write(Ffile, F)
    common.run('rectify_mindistortion %s %d %d %s %s > /dev/null' %
                                            (Ffile, w, h, Haf, Hbf))
    Ha = common.matrix_read(Haf, 3, 3)
    Hb = common.matrix_read(Hbf, 3, 3)
    return Ha, Hb
