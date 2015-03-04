# Copyright (C) 2013, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2013, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2013, Julien Michel <julien.michel@cnes.fr>

import sys
import numpy as np
import homography_cropper
import rpc_model
import rpc_utils
import estimation
import evaluation
import common
from config import cfg


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


def matches_from_sift(im1, im2, asift_only=False):
    """
    Computes a list of sift matches between two images.

    This function uses the parameter subsampling_factor_registration from the
    config module. If factor > 1 then the registration is performed over
    subsampled images, but the resulting keypoints are then scaled back to
    cancel the subsampling.

    Two implementations of sift are used: first the one by Pascal Monasse, and
    if we get less than 10 matches we use the one from ipol, which usually
    detects more keypoints.

    Args:
        im1, im2: paths to the two images
        asift_only: set to true to use ASIFT on the first try. If False, SIFT
            will be used first, and ASIFT will be used only if not enough
            matches are found

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
    """
    # convert to gray
    if common.image_pix_dim(im1) == 4:
        im1 = common.pansharpened_to_panchro(im1)
    if common.image_pix_dim(im2) == 4:
        im2 = common.pansharpened_to_panchro(im2)

    # zoom out
    zoom = cfg['subsampling_factor_registration']
    if zoom != 1:
        im1 = common.image_safe_zoom_fft(im1, zoom)
        im2 = common.image_safe_zoom_fft(im2, zoom)

    # rescale on 8 bits
    im1_8b = common.image_qauto(im1)
    im2_8b = common.image_qauto(im2)

    if not asift_only:
        p1 = common.image_sift_keypoints(im1_8b, max_nb=2000)
        p2 = common.image_sift_keypoints(im2_8b, max_nb=2000)
        matches = common.sift_keypoints_match(p1, p2, 'relative',
                                              cfg['sift_match_thresh'])

    if asift_only or matches.shape[0] < 10:
        ver = common.tmpfile('.png')
        hor = common.tmpfile('.png')
        match_f = common.tmpfile('.txt')
        common.run('demo_ASIFT %s %s %s %s %s /dev/null /dev/null' % (im1_8b,
                                                                      im2_8b,
                                                                      ver, hor,
                                                                      match_f))
        matches = np.loadtxt(match_f, skiprows=1)

    # Below is an alternative to ASIFT: lower the thresh_dog for the sift calls.
    # Default value for thresh_dog is 0.0133
    if False:
        thresh_dog = 0.0133
        nb_sift_tries = 1
        while (matches.shape[0] < 10 and nb_sift_tries < 6):
            nb_sift_tries += 1
            thresh_dog /= 2.0
            p1 = common.image_sift_keypoints(im1_8b, None, None,
                                             '-thresh_dog %f' % thresh_dog)
            p2 = common.image_sift_keypoints(im2_8b, None, None,
                                             '-thresh_dog %f' % thresh_dog)
            matches = common.sift_keypoints_match(p1, p2, 'relative',
                                                  cfg['sift_match_thresh'])

    # compensate coordinates for the zoom
    return matches * zoom


def matches_from_sift_rpc_roi(im1, im2, rpc1, rpc2, x, y, w, h, asift_only=False):
    """
    Computes a list of sift matches between two Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers defining the rectangular ROI in the first
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.

        This function uses the parameter subsampling_factor_registration
        from the config module. If factor > 1 then the registration
        is performed over subsampled images, but the resulting keypoints
        are then scaled back to conceal the subsampling

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the big images.
    """
    x1, y1, w1, h1 = x, y, w, h
    x2, y2, w2, h2 = rpc_utils.corresponding_roi(rpc1, rpc2, x, y, w, h)

    # do crops, to apply sift on reasonably sized images
    crop1 = common.image_crop_LARGE(im1, x1, y1, w1, h1)
    crop2 = common.image_crop_LARGE(im2, x2, y2, w2, h2)
    T1 = common.matrix_translation(x1, y1)
    T2 = common.matrix_translation(x2, y2)

    # compute sift matches for the images
    matches = matches_from_sift(crop1, crop2, asift_only)

    if matches.size:
        # compensate coordinates for the crop and the zoom
        pts1 = common.points_apply_homography(T1, matches[:, 0:2])
        pts2 = common.points_apply_homography(T2, matches[:, 2:4])
        return np.hstack([pts1, pts2])
    else:
        return np.array([[]])


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
    if not matches.size:
        return np.array(out)
    for i in range(len(matches)):
        x = np.array([matches[i, 0], matches[i, 1], 1])
        xx = np.array([matches[i, 2], matches[i, 3], 1])
        l = np.dot(F.T, xx)
        ll = np.dot(F, x)
        d1 = evaluation.distance_point_to_line(x, l)
        d2 = evaluation.distance_point_to_line(xx, ll)
        d = max(d1, d2)
        if (d < thresh):
            out.append(matches[i, :])

    return np.array(out)


def register_horizontally(matches, H1, H2, do_shear=True, flag='center'):
    """
    Adjust rectifying homographies to modify the disparity range.

    Args:
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        H1, H2: two homographies, stored as numpy 3x3 matrices
        do_shear: boolean flag indicating wheter to minimize the shear on im2
            or not.
        flag: option needed to control how to modify the disparity range:
            'center': move the barycenter of disparities of matches to zero
            'positive': make all the disparities positive
            'negative': make all the disparities negative. Required for
                Hirshmuller stereo (java)

    Returns:
        H2: corrected homography H2

    The matches are provided in the original images coordinate system. By
    transforming these coordinates with the provided homographies, we obtain
    matches whose disparity is only along the x-axis. The second homography H2
    is corrected with a horizontal translation to obtain the desired property
    on the disparity range.
    """
    # transform the matches according to the homographies
    pt1 = common.points_apply_homography(H1, matches[:, 0:2])
    x1 = pt1[:, 0]
    y1 = pt1[:, 1]
    pt2 = common.points_apply_homography(H2, matches[:, 2:4])
    x2 = pt2[:, 0]
    y2 = pt2[:, 1]

    # for debug, print the vertical disparities. Should be zero.
    print "Residual vertical disparities: max, min, mean. Should be zero ------"
    print np.max(y2 - y1), np.min(y2 - y1), np.mean(y2 - y1)

    # shear correction
    # we search the (s, b) vector that minimises \sum (x1 - (x2+s*y2+b))^2
    # it is a least squares minimisation problem
    if do_shear:
        A = np.vstack((y2, y2*0+1)).T
        B = x1 - x2
        z = np.linalg.lstsq(A, B)[0]
        s = z[0]
        b = z[1]
        H2 = np.dot(np.array([[1, s, b], [0, 1, 0], [0, 0, 1]]), H2)
        x2 = x2 + s*y2 + b

    # compute the disparity offset according to selected option
    if (flag == 'center'):
        t = np.mean(x2 - x1)
    if (flag == 'positive'):
        t = np.min(x2 - x1)
    if (flag == 'negative'):
        t = np.max(x2 - x1)
    if (flag == 'none'):
        t = 0

    # correct H2 with a translation
    H2 = np.dot(common.matrix_translation(-t, 0), H2)
    return H2


def update_disp_range(matches, H1, H2, w, h):
    """
    Update the disparity range considering the extrapolation of the affine
    registration estimated from the matches. Extrapolate on the whole region
    of the region of interest.

    Args:
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        H1, H2: two rectifying homographies, stored as numpy 3x3 matrices
        w, h: width and height of the region of interest

    Returns:
        disp_min, disp_max: horizontal disparity range
    """
    # transform the matches according to the homographies
    pt1 = common.points_apply_homography(H1, matches[:, 0:2])
    x1 = pt1[:, 0]
    pt2 = common.points_apply_homography(H2, matches[:, 2:4])
    x2 = pt2[:, 0]
    y2 = pt2[:, 1]

    # estimate an affine transformation (tilt, shear and bias)
    # that maps pt1 on pt2
    A = np.vstack((x2, y2, y2*0+1)).T
    B = x1
    z = np.linalg.lstsq(A, B)[0]
    t, s, b = z[0:3]

    # corners of ROI
    xx2 = np.array([0, w, 0, w])
    yy2 = np.array([0, 0, h, h])

    # compute the max and min disparity values according to the estimated
    # model. The min/max disp values are necessarily obtained at the ROI
    # corners
    roi_disparities_by_the_affine_model = (xx2*t + yy2*s + b) - xx2
    max_roi = np.max(roi_disparities_by_the_affine_model)
    min_roi = np.min(roi_disparities_by_the_affine_model)

    # min/max disparities according to the keypoints
    max_kpt = np.max(x2 - x1)
    min_kpt = np.min(x2 - x1)

    # compute the range with the extracted min and max disparities
    dispx_min = np.floor(min(min_roi, min_kpt))
    dispx_max = np.floor(max(max_roi, max_kpt))

    # add a security margin to the disp range
    d = cfg['disp_range_extra_margin']
    if (dispx_min < 0):
        dispx_min = (1 + d) * dispx_min
    else:
        dispx_min = (1 - d) * dispx_min
    if (dispx_max > 0):
        dispx_max = (1 + d) * dispx_max
    else:
        dispx_max = (1 - d) * dispx_max

    return dispx_min, dispx_max


def compute_rectification_homographies(im1, im2, rpc1, rpc2, x, y, w, h, A=None,
                                       m=None):
    """
    Computes rectifying homographies for a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2. This matrix is usually estimated with the pointing_accuracy
            module.
        m (optional): Nx4 numpy array containing a list of matches.

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies to be
            applied to the two images.
        disp_min, disp_max: horizontal disparity range, computed on a set of
            sift matches
    """
    # in brief: use 8-pts normalized algo to estimate F, then use loop-zhang to
    # estimate rectifying homographies.

    print "step 1: find virtual matches, and center them ----------------------"
    n = cfg['n_gcp_per_axis']
    rpc_matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, n)
    p1 = rpc_matches[:, 0:2]
    p2 = rpc_matches[:, 2:4]

    if A is not None:
        print "applying pointing error correction"
        # correct coordinates of points in im2, according to A
        p2 = common.points_apply_homography(np.linalg.inv(A), p2)

    # the matching points are translated to be centered in 0, in order to deal
    # with coordinates ranging from -1000 to 1000, and decrease imprecision
    # effects of the loop-zhang rectification. These effects may become very
    # important (~ 10 pixels error) when using coordinates around 20000.
    pp1, T1 = center_2d_points(p1)
    pp2, T2 = center_2d_points(p2)

    print "step 2: estimate F (Gold standard algorithm) -----------------------"
    F = estimation.affine_fundamental_matrix(np.hstack([pp1, pp2]))

    print "step 3: compute rectifying homographies (loop-zhang algorithm) -----"
    H1, H2 = estimation.loop_zhang(F, w, h)
    S1, S2 = estimation.rectifying_similarities_from_affine_fundamental_matrix(
        F, True)
    print "F\n", F, "\n"
    print "H1\n", H1, "\n"
    print "S1\n", S1, "\n"
    print "H2\n", H2, "\n"
    print "S2\n", S2, "\n"
    # compose with previous translations to get H1, H2 in the big images frame
    H1 = np.dot(H1, T1)
    H2 = np.dot(H2, T2)

    # for debug
    print "max, min, mean rectification error on rpc matches ------------------"
    tmp = common.points_apply_homography(H1, p1)
    y1 = tmp[:, 1]
    tmp = common.points_apply_homography(H2, p2)
    y2 = tmp[:, 1]
    err = np.abs(y1 - y2)
    print np.max(err), np.min(err), np.mean(err)

    print "step 4: pull back top-left corner of the ROI in the origin ---------"
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts = common.points_apply_homography(H1, roi)
    x0, y0 = common.bounding_box2D(pts)[0:2]
    T = common.matrix_translation(-x0, -y0)
    H1 = np.dot(T, H1)
    H2 = np.dot(T, H2)

    # add an horizontal translation to H2 to center the disparity range around
    # the origin, if sift matches are available
    if m is not None:
        print "step 5: horizontal registration --------------------------------"

        # filter sift matches with the known fundamental matrix
        # but first convert F for big images coordinate frame
        F = np.dot(T2.T, np.dot(F, T1))
        print '%d sift matches before epipolar constraint filering', len(m)
        m = filter_matches_epipolar_constraint(F, m, cfg['epipolar_thresh'])
        print 'd sift matches after epipolar constraint filering', len(m)
        if len(m) < 2:
            # 0 or 1 sift match
            print 'rectification.compute_rectification_homographies: less than'
            print '2 sift matches after filtering by the epipolar constraint.'
            print 'This may be due to the pointing error, or to strong'
            print 'illumination changes between the input images.'
            print 'No registration will be performed.'
        else:
            H2 = register_horizontally(m, H1, H2)
            disp_m, disp_M = update_disp_range(m, H1, H2, w, h)

    # expand disparity range with srtm according to cfg params
    if (cfg['disp_range_method'] is "srtm") or (m is None) or (len(m) < 2):
        disp_m, disp_M = rpc_utils.srtm_disp_range_estimation(
            rpc1, rpc2, x, y, w, h, H1, H2, A,
            cfg['disp_range_srtm_high_margin'],
            cfg['disp_range_srtm_low_margin'])
    if ((cfg['disp_range_method'] is "wider_sift_srtm") and (m is not None) and
            (len(m) >= 2)):
        d_m, d_M = rpc_utils.srtm_disp_range_estimation(
            rpc1, rpc2, x, y, w, h, H1, H2, A,
            cfg['disp_range_srtm_high_margin'],
            cfg['disp_range_srtm_low_margin'])
        disp_m = min(disp_m, d_m)
        disp_M = max(disp_M, d_M)

    print "disparity range:  [%s, %s]" % (disp_m, disp_M)
    return H1, H2, disp_m, disp_M


def rectify_pair(im1, im2, rpc1, rpc2, x, y, w, h, out1, out2, A=None, m=None,
                 flag='rpc'):
    """
    Rectify a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: paths to the two xml files containing RPC data
        x, y, w, h: four integers defining the rectangular ROI in the first
            image.  (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        out1, out2: paths to the output crops
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2. This matrix is usually estimated with the pointing_accuracy
            module.
        m (optional): Nx4 numpy array containing a list of sift matches, in the
            full image coordinates frame
        flag (default: 'rpc'): option to decide wether to use rpc of sift
            matches for the fundamental matrix estimation.

        This function uses the parameter subsampling_factor from the
        config module.  If the factor z > 1 then the output images will
        be subsampled by a factor z.  The output matrices H1, H2, and the
        ranges are also updated accordingly:
        Hi = Z*Hi   with Z = diag(1/z,1/z,1)   and
        disp_min = disp_min/z  (resp _max)

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies that
            have been applied to the two (big) images.
        disp_min, disp_max: horizontal disparity range
    """
    # read RPC data
    rpc1 = rpc_model.RPCModel(rpc1)
    rpc2 = rpc_model.RPCModel(rpc2)

    # compute rectifying homographies
    if flag == 'rpc':
        H1, H2, disp_min, disp_max = compute_rectification_homographies(
            im1, im2, rpc1, rpc2, x, y, w, h, A, m)
    else:
        H1, H2, disp_min, disp_max = compute_rectification_homographies_sift(
            im1, im2, rpc1, rpc2, x, y, w, h)

    # compute output images size
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts1 = common.points_apply_homography(H1, roi)
    x0, y0, w0, h0 = common.bounding_box2D(pts1)
    # check that the first homography maps the ROI in the positive quadrant
    assert(round(x0) == 0)
    assert(round(y0) == 0)

    # apply homographies and do the crops
    homography_cropper.crop_and_apply_homography(out1, im1, H1, w0, h0,
                                                 cfg['subsampling_factor'],
                                                 True)
    homography_cropper.crop_and_apply_homography(out2, im2, H2, w0, h0,
                                                 cfg['subsampling_factor'],
                                                 True)

    #  If subsampling_factor'] the homographies are altered to reflect the zoom
    if cfg['subsampling_factor'] != 1:
        from math import floor, ceil
        # update the H1 and H2 to reflect the zoom
        Z = np.eye(3)
        Z[0, 0] = Z[1, 1] = 1.0 / cfg['subsampling_factor']

        H1 = np.dot(Z, H1)
        H2 = np.dot(Z, H2)
        disp_min = floor(disp_min / cfg['subsampling_factor'])
        disp_max = ceil(disp_max / cfg['subsampling_factor'])

    return H1, H2, disp_min, disp_max


def compute_rectification_homographies_sift(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Computes rectifying homographies for a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies to be
            applied to the two images.
        disp_min, disp_max: horizontal disparity range, computed on a set of
            sift matches
    """
    # in brief: use ransac to estimate F from a set of sift matches, then use
    # loop-zhang to estimate rectifying homographies.

    matches = matches_from_sift(im1, im2, rpc1, rpc2, x, y, w, h)
    p1 = matches[:, 0:2]
    p2 = matches[:, 2:4]

    # the matching points are translated to be centered in 0, in order to deal
    # with coordinates ranging from -1000 to 1000, and decrease imprecision
    # effects of the loop-zhang rectification. These effects may become very
    # important (~ 10 pixels error) when using coordinates around 20000.
    pp1, T1 = center_2d_points(p1)
    pp2, T2 = center_2d_points(p2)

    F = estimation.fundamental_matrix_ransac(np.hstack([pp1, pp2]))
    H1, H2 = estimation.loop_zhang(F, w, h)

    # compose with previous translations to get H1, H2 in the big images frame
    H1 = np.dot(H1, T1)
    H2 = np.dot(H2, T2)

    # for debug
    print "max, min, mean rectification error on sift matches ----------------"
    tmp = common.points_apply_homography(H1, p1)
    y1 = tmp[:, 1]
    tmp = common.points_apply_homography(H2, p2)
    y2 = tmp[:, 1]
    err = np.abs(y1 - y2)
    print np.max(err), np.min(err), np.mean(err)

    # pull back top-left corner of the ROI in the origin
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts = common.points_apply_homography(H1, roi)
    x0, y0 = common.bounding_box2D(pts)[0:2]
    T = common.matrix_translation(-x0, -y0)
    H1 = np.dot(T, H1)
    H2 = np.dot(T, H2)

    # add an horizontal translation to H2 to center the disparity range around
    H2 = register_horizontally(matches, H1, H2)
    disp_m, disp_M = update_disp_range(matches, H1, H2, w, h)

    return H1, H2, disp_m, disp_M


def main():

    if len(sys.argv) > 12:
        im1 = sys.argv[1]
        im2 = sys.argv[2]
        rpc1 = sys.argv[3]
        rpc2 = sys.argv[4]
        x = int(sys.argv[5])
        y = int(sys.argv[6])
        w = int(sys.argv[7])
        h = int(sys.argv[8])
        H1f = sys.argv[9]
        H2f = sys.argv[10]
        out1 = sys.argv[11]
        out2 = sys.argv[12]
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
if __name__ == '__main__':
    main()
