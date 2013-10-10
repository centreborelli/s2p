#!/usr/bin/env python

import numpy as np
import common
import homography_cropper

def compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height, rpc_err):
    """
    Computes a height map from a disparity map, using rpc.

    Args:
        rpc1, rpc2: paths to the xml files
        H1, H2: paths to the files containing the homography matrices
        disp, mask: paths to the diparity and mask maps
        height: path to the output height map
        rpc_err: path to the output rpc_error of triangulation
    """
    common.run("disp_to_h %s %s %s %s %s %s %s %s" % (rpc1, rpc2, H1, H2,
        disp, mask, height, rpc_err))
    return


def transfer_height_map(height, msk, H, rpc, x, y, w, h, zoom, out_height,
    out_msk):
    """
    Transfer the heights computed on the rectified grid to the original
    Pleiades image grid.

    Args:
        height: path to the input height map (on the rectified grid)
        msk: path to the associated mask
        H: path to the file containing the rectifying homography matrix
        rpc: path to the xml file
        x, y, w, h: four integers defining the rectangular ROI in the original
            image.  (x, y) is the top-left corner, and (w, h) are the dimensions of
            the rectangle.
        zoom: zoom factor (usually 1, 2 or 4) used to produce the input height
            map
        out_height: path to the output height map
        out_msk: path to the output mask
    """
    # write the matrix associated to crop
    A = common.matrix_translation(-x, -y)
    f = 1.0/zoom
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    H_crop = common.tmpfile('.txt')
    np.savetxt(H_crop, A)

    # run the height_rpc_move binary
    tmp_h = common.tmpfile('.tif')
    common.run("height_rpc_move %s %s %s %s %s %s %s %s %d %d" % (rpc, H,
        height, msk, rpc, H_crop, tmp_h, out_msk, w*f, h*f))

    # replace the -inf with nan
    # implements: if isinf(x) then nan, else x
    common.run('plambda %s "x isinf nan x if" > %s' % (tmp_h, out_height))
    return


def colorize(crop_panchro, im_color, H, out_colorized):
    """
    Colorizes a Pleiades gray crop using low-resolution color information.

    Args:
        crop_panchro: path to the panchro (ie gray) rectified crop
        im_color: path to the full color image (tiff or jp2)
        H: path to the file containing the coefficients of the rectifying
            homography, that was used to generate crop_panchro
        out_colorized: path to the output file
    """
    # 1. Get a rectified and zoomed crop from the color image. It has to be
    # sampled on exactly the same grid as the panchro rectified crop. To do
    # that we compose the rectifying homography with a 4x zoom (because color
    # pleiades images have 4x lower resolution).
    # There is also a small horizontal translation (4 pixels at the panchro
    # resolution)
    H = np.loadtxt(H)
    H_zoom = np.array([[4, 0, -4], [0, 4, 0], [0, 0, 1]])
    H = np.dot(H, H_zoom)
    w, h = common.image_size(crop_panchro)
    crop_ms = common.tmpfile('.tif')
    homography_cropper.crop_and_apply_homography(crop_ms, im_color, H, w, h)

    # convert rgbi to rgb and requantify between 0 and 255
    crop_rgb = common.rgbi_to_rgb(crop_ms)
    #rgb      = common.image_qeasy(crop_rgb, 300, 3000)
    #panchro  = common.image_qeasy(crop_panchro, 300, 3000)
    rgb      = common.image_qauto(crop_rgb)
    panchro  = common.image_qauto(crop_panchro)

    # 2. Combine linearly the intensity and the color to obtain the result
    common.run('plambda %s %s "dup split + + / *" | qeasy 0 85 - %s' %
        (panchro, rgb, out_colorized))
    return



def compute_point_cloud(crop_colorized, heights, rpc, H, cloud):
    """
    Computes a color point cloud from a height map.

    Args:
        crop_colorized: path to the colorized rectified crop
        heights: height map. Its size is the same as the crop_color image
        rpc: path to xml file containing RPC data for the current Pleiade image
        H: path to the file containing the coefficients of the rectifying
            homography
        cloud: path to the output points cloud (ply format)
    """
    common.run("colormesh %s %s %s %s %s" % (crop_colorized, heights, rpc, H,
        cloud))
    return
