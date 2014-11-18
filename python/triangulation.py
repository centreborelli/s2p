# Copyright (C) 2013, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2013, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2013, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import numpy as np

import common
from config import cfg


def compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height, rpc_err, A=None):
    """
    Computes a height map from a disparity map, using rpc.

    Args:
        rpc1, rpc2: paths to the xml files
        H1, H2: path to txt files containing two 3x3 numpy arrays defining
            the rectifying homographies
        disp, mask: paths to the diparity and mask maps
        height: path to the output height map
        rpc_err: path to the output rpc_error of triangulation
        A (optional): pointing correction matrix for im2
    """
    if A is not None:
        HH2 = common.tmpfile('.txt')
        np.savetxt(HH2, np.dot(np.loadtxt(H2), np.linalg.inv(A)))
    else:
        HH2 = H2

    common.run("disp_to_h %s %s %s %s %s %s %s %s" % (rpc1, rpc2, H1, HH2, disp,
                                                      mask, height, rpc_err))
    return


def transfer_map(in_map, H, x, y, w, h, zoom, out_map):
    """
    Transfer the heights computed on the rectified grid to the original
    Pleiades image grid.

    Args:
        in_map: path to the input map, usually a height map or a mask, sampled
            on the rectified grid
        H: path to txt file containing a numpy 3x3 array representing the
            rectifying homography
        x, y, w, h: four integers defining the rectangular ROI in the original
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        zoom: zoom factor (usually 1, 2 or 4) used to produce the input height
            map
        out_map: path to the output map
    """
    # write the inverse of the resampling transform matrix. In brief it is:
    # homography * translation * zoom
    # This matrix transports the coordinates of the original cropped and
    # zoomed grid (the one desired for out_height) to the rectified cropped and
    # zoomed grid (the one we have for height)
    Z = np.diag([zoom, zoom, 1])
    A = common.matrix_translation(x, y)
    HH = np.dot(np.loadtxt(H), np.dot(A, Z))

    # apply the homography
    # write the 9 coefficients of the homography to a string, then call synflow
    # to produce the flow, then backflow to apply it
    # zero:256x256 is the iio way to create a 256x256 image filled with zeros
    hij = ' '.join(['%r' % num for num in HH.flatten()])
    common.run('synflow hom "%s" zero:%dx%d /dev/null - | BILINEAR=1 backflow - %s %s' % (
        hij, w/zoom, h/zoom, in_map, out_map))

    # w and h, provided by the s2p.process_pair_single_tile function, are
    # always multiples of zoom.

    # replace the -inf with nan in the heights map, because colormesh filter
    # out nans but not infs
    # implements: if isinf(x) then nan, else x
    # common.run('plambda %s "x isinf nan x if" > %s' % (tmp_h, out_height))


def compute_dem(out, x, y, w, h, z, rpc1, rpc2, H1, H2, disp, mask, rpc_err,
                A=None):
    """
    Computes an altitude map, on the grid of the original reference image, from
    a disparity map given on the grid of the rectified reference image.

    Args:
        out: path to the output file
        x, y, w, h: four integers defining the rectangular ROI in the original
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        z: zoom factor (usually 1, 2 or 4) used to produce the input disparity
            map
        rpc1, rpc2: paths to the xml files
        H1, H2: path to txt files containing two 3x3 numpy arrays defining
            the rectifying homographies
        disp, mask: paths to the diparity and mask maps
        rpc_err: path to the output rpc_error of triangulation
        A (optional): pointing correction matrix for im2

    Returns:
        nothing
    """
    out_dir = os.path.dirname(out)

    # redirect stdout and stderr to log file, in append mode
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'a', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    tmp = common.tmpfile('.tif')
    compute_height_map(rpc1, rpc2, H1, H2, disp, mask, tmp, rpc_err, A)
    transfer_map(tmp, H1, x, y, w, h, z, out)

    # close logs
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()


def compute_ply(out, rpc1, rpc2, H1, H2, disp, mask, img, A=None):
    """
    Computes a 3D point cloud from a disparity map.

    Args:
        out: path to the output ply file
        rpc1, rpc2: paths to the xml files
        H1, H2: path to txt files containing two 3x3 numpy arrays defining
            the rectifying homographies
        disp, mask: paths to the diparity and mask maps
        img: path to the png image containing the colors
        A (optional): pointing correction matrix for im2
    """
    # redirect stdout and stderr to log file
    if not cfg['debug']:
        log_file = '%s/stdout.log' % os.path.dirname(out)
        fout = open(log_file, 'a', 0)  # 'a' for append, 0 for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # apply correction matrix
    if A is not None:
        HH2 = '%s/H_sec_corrected.txt' % os.path.dirname(out)
        np.savetxt(HH2, np.dot(np.loadtxt(H2), np.linalg.inv(A)))
    else:
        HH2 = H2

    # do the job
    common.run("disp2ply %s %s %s %s %s %s %s %s" % (out, disp,  mask, H1, HH2,
                                                     rpc1, rpc2, img))
    # close logs
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return


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
    w, h = common.image_size_tiffinfo(crop_panchro)
    xx = np.floor(x / 4.0) + 1
    yy = np.floor(y / 4.0)
    ww = np.ceil((x + w * zoom) / 4.0) - xx + 1
    hh = np.ceil((y + h * zoom) / 4.0) - yy
    crop_ms = common.image_crop_TIFF(im_color, xx, yy, ww, hh)
    crop_ms = common.image_zoom_gdal(crop_ms, zoom/4.0)
    # crop_ms = common.image_safe_zoom_fft(crop_ms, zoom/4.0)

    # crop the crop_ms image to remove the extra-pixels due to the integer crop
    # followed by zoom
    x0 = x - 4*xx
    y0 = y - 4*yy
    crop_ms = common.image_crop_TIFF(crop_ms, x0, y0, w, h)
    assert(common.image_size_gdal(crop_panchro) ==
           common.image_size_gdal(crop_ms))

    # convert rgbi to rgb and requantify between 0 and 255
    crop_rgb = common.rgbi_to_rgb(crop_ms)
    # rgb      = common.image_qeasy(crop_rgb, 300, 3000)
    # panchro  = common.image_qeasy(crop_panchro, 300, 3000)
    rgb = common.image_qauto_gdal(crop_rgb)
    panchro = common.image_qauto_gdal(crop_panchro)

    # 2. Combine linearly the intensity and the color to obtain the result
    common.run('plambda %s %s "dup split + + / *" | qeasy 0 85 - %s' % (
        panchro, rgb, out_colorized))
    return


def compute_point_cloud(crop_colorized, heights, rpc, H, cloud, off_x=None,
                        off_y=None, ascii_ply=False, with_normals=False):
    """
    Computes a color point cloud from a height map.

    Args:
        crop_colorized: path to a colorized crop of a Pleiades image
        heights: height map, sampled on the same grid as the crop_colorized
            image. In particular, its size is the same as crop_colorized.
        rpc: path to xml file containing RPC data for the current Pleiade image
        H: path to the file containing the coefficients of the homography
            transforming the coordinates system of the original full size image
            into the coordinates system of the crop we are dealing with.
        cloud: path to the output points cloud (ply format)
        off_{x,y} (optional, default None): coordinates of the point we want to
            use as origin in the local coordinate system of the computed cloud
        ascii_ply (optional, default false): boolean flag to tell if the output
            ply file should be encoded in plain text (ascii).
    """
    hom = np.loadtxt(H)
    hij = ' '.join(['%f' % x for x in hom.flatten()])
    asc = "--ascii" if ascii_ply else ""
    nrm = "--with-normals" if with_normals else ""
    command = "colormesh %s %s %s %s -h \"%s\" %s %s" % (cloud, heights, rpc,
                                                         crop_colorized, hij,
                                                         asc, nrm)
    if off_x:
        command += " --offset_x %d" % off_x
    if off_y:
        command += " --offset_y %d" % off_y
    common.run(command)

    # if LidarViewer is installed, convert the point cloud to its format
    # this is useful for huge point clouds
    if common.which('LidarPreprocessor'):
        cloud_lidar_viewer = "%s.lidar_viewer" % os.path.splitext(cloud)[0]
        common.run("LidarPreprocessor %s -o %s" % (cloud, cloud_lidar_viewer))
    return
