#!/usr/bin/env python

import numpy as np
import common
import homography_cropper

def compute_height_map(rpc1, H1, rpc2, H2, disp, mask, im1, cloud, height):
    common.run("disp_to_h %s %s %s %s %s %s %s %s %s" % (rpc1, H1, rpc2, H2,
                                            disp, mask, im1, cloud, height))
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
    # pleiades images have 4x lower resolution)
    H = np.loadtxt(H)
    H_zoom = np.array([[4, 0, 0], [0, 4, 0], [0, 0, 1]])
    H = np.dot(H, H_zoom)
    w, h = common.image_size(crop_panchro)
    crop_ms = common.tmpfile('.tif')
    homography_cropper.crop_and_apply_homography_ms(crop_ms, im_color, H, w, h)

    # convert rgbi to rgb and requantify between 0 and 255
    crop_rgb = common.rgbi_to_rgb(crop_ms)
    rgb      = common.image_qeasy(crop_rgb, 300, 3000)
    panchro  = common.image_qeasy(crop_panchro, 300, 3000)

    # 2. Combine linearly the intensity and the color to obtain the result
    common.run('plambda %s %s "dup split + + / *" | qeasy 0 85 - %s' % (panchro,
                                                            rgb, out_colorized))
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
    common.run("colormeshh %s %s %s %s > %s" % (crop_colorized, heights, rpc,
                                                                    H, cloud))
    return
