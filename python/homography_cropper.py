#!/usr/bin/env python

import numpy as np
import common

def image_apply_pleiades_unsharpening_filter(im):
    """
    Returns the image convolved by the unsharpening MTF idata_0009_MTF_89x89.tif
    This filter specifically undoes the sharpening applied to the sensor perfect
    Pleiades images
    """
    unsharpening_mtf_small = common.image_pleiades_unsharpening_mtf()
    tmp_mtf_large = common.image_zeropadding_from_image_with_target_size(
            unsharpening_mtf_small, im)
    return common.image_fftconvolve(im, tmp_mtf_large)


def crop_and_apply_homography(im_out, im_in, H, w, h, subsampling_factor=1):
    """
    Warps a piece of a Pleiades (panchro or ms) image with a homography.

    Args:
        im_out: path to the output image
        im_in: path to the input (jp2 or tif) full Pleiades image
        H: numpy array containing the 3x3 homography matrix
        w, h: size of the output image
        subsampling_factor (optional, default=1): when set to z>1
           will result in the application of the homography Z*H
           where Z = diag(1/z,1/z,1)
           So the output will be zoomed out by a factor z
           The output image will be (w/z, h/z)

    Returns:
        nothing

    The homography has to be used as: coord_out = H coord_in. The produced
    output image corresponds to coord_out in [0, w] x [0, h]. The warp is made
    by Pascal Monasse's binary named 'homography'.
    """

    # crop a piece of the big input image, to which the homography will be
    # applied
    pts = [[0, 0], [w, 0], [w, h], [0, h]]
    inv_H_pts = common.points_apply_homography(np.linalg.inv(H), pts)
    x0, y0, w0, h0 = common.bounding_box2D(inv_H_pts)
    tmp = common.image_crop_LARGE(im_in, x0, y0, w0, h0)

    # compensate the homography with the translation induced by the preliminary
    # crop, then apply the homography and crop.
    H = np.dot(H, common.matrix_translation(x0, y0))

    # This filter is needed (for panchro images) because the original PLEAIDES
    # SENSOR PERFECT images are aliased
    if (common.image_pix_dim(tmp) == 1 and subsampling_factor == 1):
        tmp = image_apply_pleiades_unsharpening_filter(tmp)

    # the output image is zoomed out by subsampling_factor so
    # H, w, and h are updated accordingly
    assert(subsampling_factor >= 1)
    Z = np.eye(3);
    Z[0,0] = Z[1,1] = 1.0/subsampling_factor
    H = np.dot(Z, H)

    w = int(w/subsampling_factor)
    h = int(h/subsampling_factor)


    # Since the objective is to commpute a zoomed out homographic application
    # to save computations we zoom out the image before applying the homography
    # and then update H accordingly
    if subsampling_factor != 1:
        # the DCT zoom is NOT SAFE, when the input image size is not a multiple
        # of the zoom factor
        tmpw, tmph = common.image_size(tmp)
        tmpw, tmph = int(tmpw/subsampling_factor), int(tmph/subsampling_factor)
        tmp  = common.image_crop(tmp, 0, 0, tmpw*subsampling_factor, tmph*subsampling_factor)
        # zoom out the input image
        tmp = common.image_safe_zoom_fft(tmp, subsampling_factor)

        # TODO DCT IS STILL NOT SAFE THE position 0,0 is translated half pixel !

        # update H
        Z = np.eye(3)
        Z[0, 0] = Z[1, 1] = 1.0*subsampling_factor
        H = np.dot(H, Z)

    common.image_apply_homography(im_out, tmp, H, w, h)
    return
