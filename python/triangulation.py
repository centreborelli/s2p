#!/usr/bin/env python

import numpy as np
import common
import homography_cropper

def compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height, rpc_err, A=None):
    """
    Computes a height map from a disparity map, using rpc.

    Args:
        rpc1, rpc2: paths to the xml files
        H1, H2: two 3x3 numpy arrays defining the rectifying homographies
        disp, mask: paths to the diparity and mask maps
        height: path to the output height map
        rpc_err: path to the output rpc_error of triangulation
        A (optional): pointing correction matrix for im2
    """
    # save homographies to files
    hom1 = common.tmpfile('.txt')
    hom2 = common.tmpfile('.txt')
    np.savetxt(hom1, H1)
    if A is not None:
        np.savetxt(hom2, H2.dot(np.linalg.inv(A)))
    else:
        np.savetxt(hom2, H2)

    common.run("disp_to_h %s %s %s %s %s %s %s %s" % (rpc1, rpc2, hom1, hom2,
        disp, mask, height, rpc_err))
    return


def transfer_map(in_map, ref_crop, H, x, y, zoom, out_map):
    """
    Transfer the heights computed on the rectified grid to the original
    Pleiades image grid.

    Args:
        in_map: path to the input map, usually a height map or a mask, sampled
            on the rectified grid
        ref_crop: path to the reference crop sampled on the original grid
        H: numpy 3x3 array containing the rectifying homography
        x, y: two integers defining the top-left corner of the rectangular ROI
            in the original image.
        zoom: zoom factor (usually 1, 2 or 4) used to produce the input height
            map
        out_map: path to the output map
    """
    # write the inverse of the resampling transform matrix. In brief it is:
    # homography * translation * zoom
    # This matrix transports the coordinates of the original croopped and
    # zoomed grid (the one desired for out_height) to the rectified cropped and
    # zoomed grid (the one we have for height)
    Z = np.diag([zoom, zoom, 1])
    A = common.matrix_translation(x, y)
    HH = np.dot(H, np.dot(A, Z))

    # apply the homography
    # write the 9 coefficients of the homography to a string, then call synflow
    # to produce the flow, then backflow to apply it
    hij = ' '.join(['%r' % num for num in HH.flatten()])
    common.run('synflow hom "%s" %s /dev/null - | BILINEAR=1 backflow - %s %s' % (
        hij, ref_crop, in_map, out_map))

    # replace the -inf with nan in the heights map, because colormesh filter
    # out nans but not infs
    # implements: if isinf(x) then nan, else x
    #common.run('plambda %s "x isinf nan x if" > %s' % (tmp_h, out_height))



def colorize(crop_panchro, im_color, x, y, zoom, out_colorized):
    """
    Colorizes a Pleiades gray crop using low-resolution color information.

    Args:
        crop_panchro: path to the panchro (ie gray) crop
        im_color: path to the full color image (tiff or jp2)
        x, y: coordinates of the top-left corner of crop_panchro, in the full
            Pleiade image frame.
        zoom: subsampling zoom-factor that was used to generate crop_panchro
        out_colorized: path to the output file
    """
    # 1. Get a translated and zoomed crop from the color image. It has to be
    # sampled on exactly the same grid as the panchro crop.
    # To do that we compose the translation + zoom transformation with a 4x
    # zoom (because color pleiades images have 4x lower resolution).  There is
    # also a small horizontal translation (4 pixels at the panchro resolution)
    # The resulting transformation is the composition of:
    #   translation (-1 - x/4, -y/4)
    #   zoom 4/z
    w, h = common.image_size(crop_panchro)
    ww = np.ceil(w * zoom / 4.0)
    hh = np.ceil(h * zoom / 4.0)
    xx = np.floor(1 + x / 4.0)
    yy = np.floor(y / 4.0)
    crop_ms = common.image_crop_TIFF(im_color, xx, yy, ww, hh)
    crop_ms = common.image_safe_zoom_fft(crop_ms, zoom/4.0)

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
