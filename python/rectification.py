#!/usr/bin/env python

import sys
import numpy as np
import homography_cropper
import rpc_model
import rpc_utils
import estimation
import evaluation
import common


def center_2d_points(pts):
    """
    Translates 2D points.

    The input points are translated such that the output points are centered at
    origin.

    Args:
        pts: 2D array of dimension Nx2 containing the coordinates of the input
            points, one point per line

    Returns:
        new_x, new_y, T: coordinates of the transformed points, together with
            the similarity (translation) matrix. This transformation takes the
            input points on the output points.
    """
    # centroid
    cx = np.mean(pts[:, 0])
    cy = np.mean(pts[:, 1])

    # shift origin to centroid
    new_x = pts[:, 0] - cx
    new_y = pts[:, 1] - cy

    # translation matrix
    T = np.eye(3)     #              1     0    -cx
    T[0, 2] = -cx     # matrix T  =  0     1    -cy
    T[1, 2] = -cy     #              0     0     1

    return np.vstack([new_x, new_y]).T, T


def matches_from_sift(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Computes a list of sift matches between two Pleiades images.

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
    """
    m, M = rpc_utils.altitude_range(rpc1, x, y, w, h)

    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    c = np.array([m, M,   m,   M,   m,   M,   m,   M])

    # corresponding points in im2
    xx, yy = rpc_utils.find_corresponding_point(rpc1, rpc2, a, b, c)[0:2]
    # bounding box in im2
    x2, y2, w2, h2 = common.bounding_box2D(np.vstack([xx, yy]).T)
    x1, y1, w1, h1 = x, y, w, h

    # do crops, to apply sift on reasonably sized images
    crop1 = common.image_crop_LARGE(im1, x1, y1, w1, h1)
    crop2 = common.image_crop_LARGE(im2, x2, y2, w2, h2)
    T1 = common.matrix_translation(x1, y1)
    T2 = common.matrix_translation(x2, y2)

    # apply sift, then transport keypoints coordinates in the big images frame
    kpts1 = common.image_sift_keypoints(crop1, '', 1000)
    kpts2 = common.image_sift_keypoints(crop2, '', 1000)
    matches = common.sift_keypoints_match(kpts1, kpts2, 1, 0.6)
    # 300: distance threshold for sift descriptors
    if matches.size:
        pts1 = common.points_apply_homography(T1, matches[:, 0:2])
        pts2 = common.points_apply_homography(T2, matches[:, 2:4])
        return np.hstack([pts1, pts2])
    else:
        raise Exception("no sift matches")


def filter_matches_epipolar_constraint(F, matches, thresh):
    """
    Discards matches that are not consistent with the epipolar constraint.

    Args:
        F: fundamental matrix
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        thresh: maximum accepted distance between a point and its matched
            epipolar line

    Returns:
        the list of matches that satisfy the constraint. It is a sub-list of
        the input list.
    """
    out = []
    mask = np.zeros((len(matches), 1)) # for debug only
    for i in range(len(matches)):
        x  = np.array([matches[i, 0], matches[i, 1], 1])
        xx = np.array([matches[i, 2], matches[i, 3], 1])
        l  = np.dot(F.T, xx)
        ll = np.dot(F, x)
        d1 = evaluation.distance_point_to_line(x, l)
        d2 = evaluation.distance_point_to_line(xx, ll)
        d = max(d1, d2)
        if (d < thresh):
            out.append(matches[i, :])
            mask[i] = 1 # for debug only

    np.savetxt('/tmp/sift_F_msk', mask, '%d') # for debug only
    return np.array(out)


def register_horizontally(matches, H1, H2, flag='center'):
    """
    Adjust rectifying homographies to modify the disparity range.

    Args:
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        H1, H2: two homographies, stored as numpy 3x3 matrices
        flag: option needed to control how to modify the disparity range:
            'center': move the barycenter of disparities of matches to zero
            'positive': make all the disparities positive
            'negative': make all the disparities negative. Required for
                Hirshmuller stereo (java)

    Returns:
        H2: corrected homography H2
        disp_min, disp_max: horizontal disparity range

    The matches are provided in the original images coordinate system. By
    transforming these coordinates with the provided homographies, we obtain
    matches whose disparity is only along the x-axis. The second homography H2
    is corrected with a horizontal translation to obtain the desired property
    on the disparity range.  The minimum and maximal disparities over the set
    of matches are extracted, with a security margin of 20 percent.
    """
    # transform the matches according to the homographies
    tmp = common.points_apply_homography(H1, matches[:, 0:2])
    x1 = tmp[:, 0]
    y1 = tmp[:, 1]
    tmp = common.points_apply_homography(H2, matches[:, 2:4])
    x2 = tmp[:, 0]
    y2 = tmp[:, 1]

    # compute the disparity offset according to selected option
    if (flag == 'center'):
        t = np.mean(x2 - x1)
    if (flag == 'positive'):
        t = np.min(x2 - x1)
    if (flag == 'negative'):
        t = np.max(x2 - x1)

    # correct H2 with a translation
    H2 = np.dot(common.matrix_translation(-t, 0), H2)
    x2 = x2 - t

    # extract min and max disparities
    dispx_min = np.floor((np.min(x2 - x1)))
    dispx_max = np.ceil((np.max(x2 - x1)))

    # add 20% security margin
    if (dispx_min < 0):
        dispx_min = 1.2 * dispx_min
    else:
        dispx_min = 0.8 * dispx_min
    if (dispx_max > 0):
        dispx_max = 1.2 * dispx_max
    else:
        dispx_max = 0.8 * dispx_max

    # for debug, print the vertical disparities. Should be zero.
    print "Residual vertical disparities: min, max, mean. Should be zero ------"
    print np.min(y2 - y1), np.max(y2 - y1), np.mean(y1 - y2)
    return H2, dispx_min, dispx_max



def compute_rectification_homographies(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Computes rectifying homographies for a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies to be applied
            to the two images.
        disp_min, disp_max: horizontal disparity range, computed on a set of
            sift matches
    """
    # in brief: use 8-pts normalized algo to estimate F, then use loop-zhang to
    # estimate rectifying homographies.

    # the matching points are translated to be centered in 0, in order to deal
    # with coordinates ranging from -1000 to 1000, and decrease imprecision
    # effects of the loop-zhang rectification. These effects may become very
    # important (~ 10 pixels error) when using coordinates around 20000.
    print "step 1: find matches, and center them ------------------------------"
    rpc_matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, 5)
    p1, T1 = center_2d_points(rpc_matches[:, 0:2])
    p2, T2 = center_2d_points(rpc_matches[:, 2:4])

    print "step 2: estimate F (8-points algorithm) ----------------------------"
    F = estimation.fundamental_matrix(np.hstack([p1, p2]))
#    print evaluation.fundamental_matrix(F, np.hstack([p1, p2])

    print "step 3: compute rectifying homographies (loop-zhang algorithm) -----"
    # why does loop-zhang need the size of input images ?
    H1, H2 = estimation.loop_zhang(F, w, h)
    # compose with previous translations to get H1, H2 in the big images frame
    H1 = np.dot(H1, T1)
    H2 = np.dot(H2, T2)

    # for debug
    print "min, max, mean rectification error on rpc matches ------------------"
    tmp = common.points_apply_homography(H1, rpc_matches[:, 0:2])
    y1 = tmp[:, 1]
    tmp = common.points_apply_homography(H2, rpc_matches[:, 2:4])
    y2 = tmp[:, 1]
    err = np.abs(y1 - y2)
    print np.min(err), np.max(err), np.mean(err)

    print "step 4: pull back top-left corner of the ROI in the origin ---------"
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts = common.points_apply_homography(H1, roi)
    x0, y0 = common.bounding_box2D(pts)[0:2]
    T = common.matrix_translation(-x0, -y0)
    H1 = np.dot(T, H1)
    H2 = np.dot(T, H2)

    # add an horizontal translation to H2 to center the disparity range around
    # the origin, if sift matches are available
    print "step 5: horizontal registration ------------------------------------"
    try:
        sift_matches = matches_from_sift(im1, im2, rpc1, rpc2, x, y, w, h)
    except Exception:
        print 'something failed with sift matches'
        return H1, H2

    # filter sift matches with the known fundamental matrix
    F = np.dot(T2.T, np.dot(F, T1)) # convert F for big images coordinate frame
    sift_matches = filter_matches_epipolar_constraint(F, sift_matches, 2.0)
    H2, disp_m, disp_M = register_horizontally(sift_matches, H1, H2, 'negative')

    return H1, H2, disp_m, disp_M



def rectify_pair(im1, im2, rpc1, rpc2, x, y, w, h, out1, out2):
    """
    Rectify a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: paths to the two xml files containing RPC data
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.
        out1, out2: paths to the output crops

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies that
            have been applied to the two (big) images.
        disp_min, disp_max: horizontal disparity range
    """
    # read RPC data
    rpc1 = rpc_model.RPCModel(rpc1)
    rpc2 = rpc_model.RPCModel(rpc2)

    # compute rectifying homographies
    H1, H2, disp_min, disp_max = compute_rectification_homographies(im1, im2,
                                                        rpc1, rpc2, x, y, w, h)

    # compute output images size
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts1 = common.points_apply_homography(H1, roi)
    x0, y0, w0, h0 = common.bounding_box2D(pts1)
    # check that the first homography maps the ROI in the positive quadrant
    assert (x0 == 0)
    assert (y0 == 0)

    # apply homographies and do the crops
    homography_cropper.crop_and_apply_homography(out1, im1, H1, w0, h0)
    homography_cropper.crop_and_apply_homography(out2, im2, H2, w0, h0)

    return H1, H2, disp_min, disp_max


def main():

    if len(sys.argv) > 12:
      im1  = sys.argv[1]
      im2  = sys.argv[2]
      rpc1 = sys.argv[3]
      rpc2 = sys.argv[4]
      x    = int(sys.argv[5])
      y    = int(sys.argv[6])
      w    = int(sys.argv[7])
      h    = int(sys.argv[8])
      H1f  = sys.argv[9]
      H2f  = sys.argv[10]
      out1  = sys.argv[11]
      out2  = sys.argv[12]
    else:
      print """
      Incorrect syntax, use:
        > %s im1 im2 rpc1 rpc2 x y w h H1 H2 out1 out2
        Computes rectification homographies for two pleiades images, using
        a region of interest ROI (x, y, w, h) defined over the first image.
        The area in the second image corresponding to the ROI is determined
        from the RPC metadata (height range), and, if available, from SRTM
        data. The rectified crops are computed and saved in files out1 out2
        Uses 8-points algorithm for F estimation, and Loop-Zhang for
        rectification.
        im1/2:        paths to the two Pleiades images (usually jp2 or tif)
        rpc1/2:       RPCs xml files associated to the two images
        x y w h:      ROI of first image used to compute the rectification
        H1/H2:        output rectification homographies in the
                      coordinate systems of the two BIG images
        out1/2:       paths to output crops
      """ % sys.argv[0]
      sys.exit(1)

    H1, H2, dm, dM = rectify_pair(im1, im2, rpc1, rpc2, x, y, w, h, out1, out2)
    np.savetxt(H1f, H1)
    np.savetxt(H2f, H2)

    return

# main call
if __name__ == '__main__': main()
