#!/usr/bin/env python

import sys
import numpy as np
import homography_cropper
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



def matches_from_sift(im1, im2):
    """
    Computes a list of sift matches between two images.

    Args:
        im1, im2: paths to the two images (usually jp2 or tif)

        This function uses the parameter subsampling_factor_registration
        from the global_params module. If factor > 1 then the registration
        is performed over subsampled images, but the resulting keypoints
        are then scaled back to conceal the subsampling

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the big images.
            If no sift matches are found, then an exception is raised.
    """
    ## Try to import the global parameters module
    #  it permits to pass values between different modules
    try:
        from global_params import subsampling_factor_registration
        if subsampling_factor_registration != 1:
            im1 = common.image_safe_zoom_fft(im1, subsampling_factor_registration)
            im2 = common.image_safe_zoom_fft(im2, subsampling_factor_registration)
    except ImportError:
        subsampling_factor_registration = 1


    # apply sift, then transport keypoints coordinates in the big images frame
    kpts1 = common.image_sift_keypoints(im1, '')
    kpts2 = common.image_sift_keypoints(im2, '')
    try:
        from global_params import sift_match_thresh
    except ImportError:
        sift_match_thresh = 0.6
    matches = common.sift_keypoints_match(kpts1, kpts2, 1, sift_match_thresh)
    # compensate coordinates for the crop and the zoom
    return matches*subsampling_factor_registration


def matches_from_projection_matrices_roi(im1, im2, rpc1, rpc2, x, y, w, h):
    """
    Computes a list of sift matches between two Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.

        This function uses the parameter subsampling_factor_registration
        from the global_params module. If factor > 1 then the registration
        is performed over subsampled images, but the resulting keypoints
        are then scaled back to conceal the subsampling

    Returns:
        matches: 2D numpy array containing a list of matches. Each line
            contains one pair of points, ordered as x1 y1 x2 y2.
            The coordinate system is that of the big images.
            If no sift matches are found, then an exception is raised.
    """
    #m, M = rpc_utils.altitude_range(rpc1, x, y, w, h)
    m=5
    M=20

    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    # also include the midpoints because the 8 corners of the frustum alone don't seem to work
    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w,x+w/2,x+w/2,x+w/2,x+w/2,x+w/2,x+w/2,x    ,x    ,x+w  ,x+w  ])
    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h,y    ,y    ,y+h/2,y+h/2,y+h  ,y+h  ,y+h/2,y+h/2,y+h/2,y+h/2])
    c = np.array([m, M,   m,   M,   m,   M,   m,   M,m    ,M    ,m    ,M    ,m    ,M    ,m    ,M    ,m    ,M    ])

    xx = np.zeros(len(a))
    yy = np.zeros(len(a))

    # corresponding points in im2
    P1 = np.loadtxt(rpc1)
    P2 = np.loadtxt(rpc2)
    
    M  = P1[:,:3]
    p4 = P1[:,3]
    m3 = M[2,:]

    inv_M = np.linalg.inv(M)

    v = np.vstack((a,b,c*0+1))

    for i in range(len(a)):
       v = np.array([a[i],b[i],1])
       mu = c[i] / np.sign ( np.linalg.det(M) ) 

       X3D = inv_M.dot (mu * v - p4 )
   
       # backproject 
       newpoints = P2.dot(np.hstack([X3D,1]))
       xx[i] = newpoints[0]  / newpoints[2]
       yy[i] = newpoints[1]  / newpoints[2]


    print xx
    print yy

    matches = np.vstack([a, b,xx,yy]).T
    return matches

   ##### xx, yy = rpc_utils.find_corresponding_point(rpc1, rpc2, a, b, c)[0:2]

    
    # bounding box in im2
    x2, y2, w2, h2 = common.bounding_box2D(np.vstack([xx, yy]).T) ## GF NOT USED
    x1, y1, w1, h1 = x, y, w, h
    x2, y2, w2, h2 = x, y, w, h

    # do crops, to apply sift on reasonably sized images
    crop1 = common.image_crop_LARGE(im1, x1, y1, w1, h1)
    crop2 = common.image_crop_LARGE(im2, x2, y2, w2, h2)
    T1 = common.matrix_translation(x1, y1)
    T2 = common.matrix_translation(x2, y2)

    # call sift matches for the images
    matches = matches_from_sift(crop1, crop2)

    if matches.size:
        # compensate coordinates for the crop and the zoom
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


def register_horizontally(matches, H1, H2, do_shear=True, do_scale_horizontally=False , flag='center'):
    """
    Adjust rectifying homographies to modify the disparity range.

    Args:
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        H1, H2: two homographies, stored as numpy 3x3 matrices
        do_shear: boolean flag indicating wheter to minimize the shear on im2
            or not.
        do_scale_horizontally: boolean flag indicating wheter to minimize 
            with respect to the horizontal scaling on im2 or not.
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
    pt1 = common.points_apply_homography(H1, matches[:, 0:2])
    x1 = pt1[:, 0]
    y1 = pt1[:, 1]
    pt2 = common.points_apply_homography(H2, matches[:, 2:4])
    x2 = pt2[:, 0]
    y2 = pt2[:, 1]

    # shear correction
    # we search the (s, b) vector that minimises \sum (x1 - (x2+s*y2+b))^2
    # it is a least squares minimisation problem
    if do_shear:
      # horizontal scale correction
      if do_scale_horizontally: # | x1 -  (s*x2 + t*y2 +d) |^2
          A = np.vstack((x2, y2, y2*0+1)).T
          b = x1 
          z = np.linalg.lstsq(A, b)[0]
          s = z[0]
          t = z[1]
          d = z[2]
          H2 = np.dot(np.array([[s, t, d], [0, 1, 0], [0, 0, 1]]), H2)
          x2 = s*x2  + t*y2 + d
      else:
          A = np.vstack((y2, y2*0+1)).T
          b = x1 - x2
          z = np.linalg.lstsq(A, b)[0]
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
    x2 = x2 - t

    # extract min and max disparities
    dispx_min = np.floor((np.min(x2 - x1)))
    dispx_max = np.ceil((np.max(x2 - x1)))

    # add a security margin to the disp range
    try:
        from python.global_params import disp_range_extra_margin as d
    except ImportError:
        d = 0.2
    if (dispx_min < 0):
        dispx_min = (1+d) * dispx_min
    else:
        dispx_min = (1-d) * dispx_min
    if (dispx_max > 0):
        dispx_max = (1+d) * dispx_max
    else:
        dispx_max = (1-d) * dispx_max

    # for debug, print the vertical disparities. Should be zero.
    print "Residual vertical disparities: min, max, mean. Should be zero ------"
    print np.min(y2 - y1), np.max(y2 - y1), np.mean(y1 - y2)
    return H2, dispx_min, dispx_max


def update_minmax_range_extrapolating_registration_affinity(matches, H1, H2,w_roi,h_roi):
    """
    Update the disparity range considering the extrapolation of the affine registration 
    transformation. Extrapolate until the boundary of the region of interest

    Args:
        matches: list of pairs of 2D points, stored as a Nx4 numpy array
        H1, H2: two homographies, stored as numpy 3x3 matrices
        roi_w/h: width and height of the region of interest

    Returns:
        disp_min, disp_max: horizontal disparity range
    """
    # transform the matches according to the homographies
    pt1 = common.points_apply_homography(H1, matches[:, 0:2])
    x1 = pt1[:, 0]
    y1 = pt1[:, 1]
    pt2 = common.points_apply_homography(H2, matches[:, 2:4])
    x2 = pt2[:, 0]
    y2 = pt2[:, 1]

    # estimate an affine transformation (tilt, shear and bias)
    # from the matched keypoints 
    A = np.vstack((x2, y2, y2*0+1)).T
#    A = x2[:, np.newaxis]
    b = x1
    z = np.linalg.lstsq(A, b)[0]
    t,s,dx = z[0:3]

    # corners of ROI
    xx2 = np.array([0,w_roi,0,w_roi])
    yy2 = np.array([0,0,h_roi,h_roi])

    # compute the max and min disparity values (according to
    # the estimated model) at the ROI corners
    roi_disparities_by_the_affine_model = (xx2*t + yy2*s + dx) - xx2
    maxb = np.max(roi_disparities_by_the_affine_model)
    minb = np.min(roi_disparities_by_the_affine_model)
    #print minb,maxb

    # compute the rage with the extract min and max disparities
    dispx_min = np.floor(minb + np.min(x2 - x1))
    dispx_max = np.ceil(maxb + np.max(x2 - x1))

    # add 20% security margin
    if (dispx_min < 0):
        dispx_min = 1.2 * dispx_min
    else:
        dispx_min = 0.8 * dispx_min
    if (dispx_max > 0):
        dispx_max = 1.2 * dispx_max
    else:
        dispx_max = 0.8 * dispx_max

    return dispx_min, dispx_max


def compute_rectification_homographies(im1, im2, rpc1, rpc2, x, y, w, h, A=None):
    """
    Computes rectifying homographies for a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2. This matrix is usually estimated with the pointing_accuracy
            module.

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies to be applied
            to the two images.
        disp_min, disp_max: horizontal disparity range, computed on a set of
            sift matches
    """
    # in brief: use 8-pts normalized algo to estimate F, then use loop-zhang to
    # estimate rectifying homographies.

    print "step 1: find matches, and center them ------------------------------"
    sift_matches = matches_from_projection_matrices_roi(im1, im2, rpc1, rpc2, x+w/4, y+h/4, w*2/4, h*2/4)
    #sift_matches2 = matches_from_sift(im1, im2)
    #sift_matches = sift_matches2
#    import visualisation
#    print visualisation.plot_matches(im1,im2,sift_matches)

    p1 = sift_matches[:, 0:2]
    p2 = sift_matches[:, 2:4]


    # the matching points are translated to be centered in 0, in order to deal
    # with coordinates ranging from -1000 to 1000, and decrease imprecision
    # effects of the loop-zhang rectification. These effects may become very
    # important (~ 10 pixels error) when using coordinates around 20000.
    pp1, T1 = center_2d_points(p1)
    pp2, T2 = center_2d_points(p2)

    print "step 2: estimate F (8-points algorithm) ----------------------------"
    F = estimation.fundamental_matrix(np.hstack([pp1, pp2]))
    F = np.dot(T2.T, np.dot(F, T1)) # convert F for big images coordinate frame

    print "step 3: compute rectifying homographies (loop-zhang algorithm) -----"
    H1, H2 = estimation.loop_zhang(F, w, h)
    #### ATTENTION: LOOP-ZHANG IMPLICITLY ASSUMES THAT F IS IN THE FINAL (CROPPED)
    # IMAGE GEOMETRY. THUS 0,0 IS THE UPPER LEFT CORNER OF THE IMAGE AND W,H ARE 
    # USED TO ESTIMATE THE DISTORTION WITHIN THE REGION. BY CENTERING THE COORDINATES
    # OF THE PIXELS WE ARE CONSTRUCTING A RECTIFICATION DOES NOT TAKE INTO ACCOUNT THE 
    # CORRECT IMAGE PORTION.
    # compose with previous translations to get H1, H2 in the big images frame
    #H1 = np.dot(H1, T1)
    #H2 = np.dot(H2, T2)

    # for debug
    print "min, max, mean rectification error on rpc matches ------------------"
    tmp = common.points_apply_homography(H1, p1)
    y1 = tmp[:, 1]
    tmp = common.points_apply_homography(H2, p2)
    y2 = tmp[:, 1]
    err = np.abs(y1 - y2)
    print np.min(err), np.max(err), np.mean(err)

#    print "step 4: pull back top-left corner of the ROI in the origin ---------"
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts = common.points_apply_homography(H1, roi)
    x0, y0 = common.bounding_box2D(pts)[0:2]
    T = common.matrix_translation(-x0, -y0)
    H1 = np.dot(T, H1)
    H2 = np.dot(T, H2)

    # add an horizontal translation to H2 to center the disparity range around
    # the origin, if sift matches are available
    print "step 5: horizontal registration ------------------------------------"
    sift_matches2 = matches_from_sift(im1, im2)

    # filter sift matches with the known fundamental matrix
    try:
        from python.global_params import epipolar_thresh
    except ImportError:
        epipolar_thresh = 2.0
    epipolar_thresh = 2.0
    sift_matches2 = filter_matches_epipolar_constraint(F, sift_matches2,
        epipolar_thresh)
    if not len(sift_matches2):
        print """all the sift matches have been discarded by the epipolar
        constraint. This is probably due to the pointing error. Try with a
        bigger value for epipolar_thresh."""
        sys.exit()

    H2, disp_m, disp_M = register_horizontally(sift_matches2, H1, H2, do_scale_horizontally=True)
    disp_m, disp_M = update_minmax_range_extrapolating_registration_affinity(sift_matches2,
        H1, H2, w, h)

    return H1, H2, disp_m, disp_M


def rectify_pair(im1, im2, rpc1, rpc2, x, y, w, h, out1, out2, A=None):
    """
    Rectify a ROI in a pair of Pleiades images.

    Args:
        im1, im2: paths to the two Pleiades images (usually jp2 or tif)
        rpc1, rpc2: paths to the two xml files containing RPC data
        x, y, w, h: four integers defining the rectangular ROI in the first image.
            (x, y) is the top-left corner, and (w, h) are the dimensions of the
            rectangle.
        out1, out2: paths to the output crops
        A (optional): 3x3 numpy array containing the pointing error correction
            for im2. This matrix is usually estimated with the pointing_accuracy
            module.

        This function uses the parameter subsampling_factor from the global_params module.
        If the factor z > 1 then the output images will be subsampled by a factor z.
        The output matrices H1, H2, and the ranges are also updated accordingly:
        Hi = Z*Hi   with Z = diag(1/z,1/z,1)   and
        disp_min = disp_min/z  (resp _max)

    Returns:
        H1, H2: Two 3x3 matrices representing the rectifying homographies that
            have been applied to the two (big) images.
        disp_min, disp_max: horizontal disparity range
    """

    # compute rectifying homographies
    H1, H2, disp_min, disp_max = compute_rectification_homographies(im1, im2,
        rpc1, rpc2, x, y, w, h, A)

    ## compute output images size
    roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
    pts1 = common.points_apply_homography(H1, roi)
    x0, y0, w0, h0 = common.bounding_box2D(pts1)
    #x0,y0,w0,h0 = x,y,w,h
    
    # check that the first homography maps the ROI in the positive quadrant
    assert (round(x0) == 0)
    assert (round(y0) == 0)

    ## Try to import the global parameters module
    #  it permits to pass values between different modules
    try:
        from python.global_params import subsampling_factor
    except ImportError:
        subsampling_factor = 1

    # apply homographies and do the crops 
    # THIS STEP IS HERE TO PRODUCE THE MASKS WHERE THE IMAGE IS KNOWN
    # SURE THIS IS A CRAPPY WAY TO DO THIS, WE SHOULD DEFINITIVELY DO IT 
    # SIMULTANEOUSLY WITH THE HOMOGRAPHIC TRANSFORMATION
    msk1 = common.tmpfile('.png')
    msk2 = common.tmpfile('.png')
    common.run('plambda %s "x 255" | iion - %s'%(im1, msk1))
    common.run('plambda %s "x 255" | iion - %s'%(im2, msk2))
    homography_cropper.crop_and_apply_homography(msk1, msk1, H1, w0, h0, subsampling_factor)
    homography_cropper.crop_and_apply_homography(msk2, msk2, H2, w0, h0, subsampling_factor)
    # FINALLY : apply homographies and do the crops of the images
    homography_cropper.crop_and_apply_homography(out1, im1, H1, w0, h0, subsampling_factor)
    homography_cropper.crop_and_apply_homography(out2, im2, H2, w0, h0, subsampling_factor)
    # COMBINE THE MASK TO REMOVE THE POINTS THAT FALL OUTSIDE THE IMAGE
    common.run('plambda %s %s "x 200 > y nan if" | iion - %s'%(msk1, out1, out1))
    common.run('plambda %s %s "x 200 > y nan if" | iion - %s'%(msk2, out2, out2))

#    This also does the job but when subsampling_factor != 1 it fails (segfault: homography)
#    TODO: FIX homography, maybe code a new one
#    common.image_apply_homography(out1, im1, H1, w0, h0)
#    common.image_apply_homography(out2, im2, H2, w0, h0)

    #  If subsampling_factor the homographies are altered to reflect the zoom
    if subsampling_factor != 1:
        from math import floor, ceil
        # update the H1 and H2 to reflect the zoom
        Z = np.eye(3);
        Z[0,0] = Z[1,1] = 1.0/subsampling_factor

        H1 = np.dot(Z, H1)
        H2 = np.dot(Z, H2)
        disp_min = floor(disp_min/subsampling_factor)
        disp_max = ceil(disp_max/subsampling_factor)
        w0 = w0/subsampling_factor
        h0 = h0/subsampling_factor


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
