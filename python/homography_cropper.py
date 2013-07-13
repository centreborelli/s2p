#!/usr/bin/env python

import numpy as np
import cv2
import piio
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


def crop_and_apply_homography_panchro(im_out, im_in, H, w, h):
    """
    Warps a piece of a Pleiades panchromatic (ie gray) image with a homography.

    Args:
        im_out: path to the output image
        im_in: path to the input (jp2 or tif) full Pleiades panchro image
        H: numpy array containing the 3x3 homography matrix
        w, h: size of the output image

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

    # preliminary crop
    tmp = common.image_crop_LARGE(im_in, x0, y0, w0, h0)

    # This filter is needed because the original PLEAIDES SENSOR PERFECT images
    # are aliased
    tmp = image_apply_pleiades_unsharpening_filter(tmp)

    # compensate the homography with the translation induced by the preliminary
    # crop, then apply the homography and crop.
    H = np.dot(H, common.matrix_translation(x0, y0))
    common.image_apply_homography(im_out, tmp, H, w, h)
    return


def crop_and_apply_homography_ms(im_out, im_in, H, w, h):
    """
    Warps a piece of a Pleiades multispectral (ie BGRI) image with a homography.

    Args:
        im_out: path to the output image
        im_in: path to the input (jp2 or tif) full Pleiades ms image
        H: numpy array containing the 3x3 homography matrix
        w, h: size of the output image

    Returns:
        nothing

    The homography has to be used as: coord_out = H coord_in. The produced
    output image corresponds to coord_out in [0, w] x [0, h]. The warp is made
    by the openCV 'warpPerspective' function.
    """

    # crop a piece of the big input image, to which the homography will be
    # applied
    pts = [[0, 0], [w, 0], [w, h], [0, h]]
    inv_H_pts = common.points_apply_homography(np.linalg.inv(H), pts)
    x0, y0, w0, h0 = common.bounding_box2D(inv_H_pts)
    tmp = common.image_crop_LARGE(im_in, x0, y0, w0, h0)

    # compensate the homography with the translation induced by the preliminary
    # crop
    H = np.dot(H, common.matrix_translation(x0, y0))

    # read the tmp image in a float numpy array, and convert H data_type to
    # float (needed by cv2.warpPerspective)
    src = piio.read(tmp)
    H = H.astype(np.float32)

    # do the warp, and save the image
    dest = cv2.warpPerspective(src, H, (w, h), flags = cv2.INTER_CUBIC)
    piio.write(im_out, dest)
    return
