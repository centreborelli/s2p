# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from __future__ import print_function
import numpy as np

from s2plib import common


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
    T[1, 1] = s         # matrix T           s     0   -s*cx
    T[0, 2] = -s*cx     # is given     T  =  0     s   -s*cy
    T[1, 2] = -s*cy     # by                 0     0     1

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
    U[1, 1] = s         # matrix U             s     0      0    -s*cx
    U[2, 2] = s         # is given             0     s      0    -s*cy
    U[0, 3] = -s*cx     # by this        U  =  0     0      s    -s*cz
    U[1, 3] = -s*cy     # formula              0     0      0      1
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
    for i in range(len(x)):
        A[2*i+0, 4:8] = -1*np.array([X[i, 0], X[i, 1], X[i, 2], 1])
        A[2*i+0, 8:12] = x[i, 1]*np.array([X[i, 0], X[i, 1], X[i, 2], 1])
        A[2*i+1, 0:4] = np.array([X[i, 0], X[i, 1], X[i, 2], 1])
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
    for i in range(len(matches)):
        A[i, 0:3] = np.array([pts1[i, 0], pts1[i, 1], 1]) * pts2[i, 0]
        A[i, 3:6] = np.array([pts1[i, 0], pts1[i, 1], 1]) * pts2[i, 1]
        A[i, 6:9] = np.array([pts1[i, 0], pts1[i, 1], 1])

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


def fundamental_matrix_ransac(matches, precision=1.0, return_inliers=False):
    """
    Estimates the fundamental matrix given a set of point correspondences
    between two images, using ransac.

    Arguments:
        matches: numpy 2D array of size Nx4 containing a list of pair of
            matching points. Each line is of the form x1, y1, x2, y2, where (x1,
            y1) is the point in the first view while (x2, y2) is the matching
            point in the second view.
            It can be the path to a txt file containing such an array.
        precision: optional parameter indicating the maximum error
            allowed for counting the inliers
        return_inliers: optional boolean flag to activate/deactivate inliers
            output

    Returns:
        the estimated fundamental matrix, and optionally the 2D array containing
        the inliers

    The algorithm uses ransac as a search engine.
    """
    if type(matches) is np.ndarray:
        # write a file containing the list of correspondences. The
        # expected format is a text file with one match per line: x1 y1 x2 y2
        matchfile = common.tmpfile('.txt')
        np.savetxt(matchfile, matches)
    else:
        # assume it is a path to a txt file containing the matches
        matchfile = matches

    # call ransac binary, from Enric's imscript
    inliers = common.tmpfile('.txt')
    Ffile = common.tmpfile('.txt')
    awk_command = "awk {\'printf(\"%e %e %e\\n%e %e %e\\n%e %e %e\", $3, $4, $5, $6, $7, $8, $9, $10, $11)\'}"
    common.run("ransac fmn 1000 %f 7 %s < %s | grep param | %s > %s" % (precision, inliers, matchfile, awk_command, Ffile))
    if return_inliers:
        return np.loadtxt(Ffile).transpose(), np.loadtxt(inliers)
    else:
        return np.loadtxt(Ffile).transpose()


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
    common.run('rectify_mindistortion %s %d %d %s %s > /dev/null' % (Ffile, w,
                                                                     h, Haf,
                                                                     Hbf))
    Ha = common.matrix_read(Haf, size=(3, 3))
    Hb = common.matrix_read(Hbf, size=(3, 3))

    # check if both the images are rotated
    a = does_this_homography_change_the_vertical_direction(Ha)
    b = does_this_homography_change_the_vertical_direction(Hb)
    if a and b:
        R = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
        Ha = np.dot(R, Ha)
        Hb = np.dot(R, Hb)
    return Ha, Hb


def does_this_homography_change_the_vertical_direction(H):
    d = H[1,1]
    q = H[1,2]
    s = H[2,1]
    t = H[2,2]
    return ((d+q) / (s+t)) < (q/t)


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
