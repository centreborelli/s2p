# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import os
import numpy as np

import common
from config import cfg


def height_map_rectified(rpc1, rpc2, H1, H2, disp, mask, height, rpc_err, A=None):
    """
    Computes a height map from a disparity map, using rpc.

    Args:
        rpc1, rpc2: paths to the xml files
        H1, H2: path to txt files containing two 3x3 numpy arrays defining
            the rectifying homographies
        disp, mask: paths to the diparity and mask maps
        height: path to the output height map
        rpc_err: path to the output rpc_error of triangulation
        A (optional): path to txt file containing the pointing correction matrix
            for im2
    """
    if A is not None:
        HH2 = common.tmpfile('.txt')
        np.savetxt(HH2, np.dot(np.loadtxt(H2), np.linalg.inv(np.loadtxt(A))))
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


def height_map(out, x, y, w, h, z, rpc1, rpc2, H1, H2, disp, mask, rpc_err,
               out_filt, A=None):
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
        A (optional): path to txt file containing the pointing correction matrix
            for im2
    """
    tmp = common.tmpfile('.tif')
    height_map_rectified(rpc1, rpc2, H1, H2, disp, mask, tmp, rpc_err, A)
    transfer_map(tmp, H1, x, y, w, h, z, out)

    # apply output filter
    common.run('plambda {0} {1} "x 0 > y nan if" -o {1}'.format(out_filt, out))


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
    # apply correction matrix
    if A is not None:
        HH2 = '%s/H_sec_corrected.txt' % os.path.dirname(out)
        np.savetxt(HH2, np.dot(np.loadtxt(H2), np.linalg.inv(A)))
    else:
        HH2 = H2

    # do the job
    common.run("disp2ply %s %s %s %s %s %s %s %s" % (out, disp,  mask, H1, HH2,
                                                     rpc1, rpc2, img))


def compute_point_cloud(cloud, heights, rpc, H=None, crop_colorized='',
                        off_x=None, off_y=None, ascii_ply=False,
                        with_normals=False, utm_zone=None, llbbx=None):
    """
    Computes a color point cloud from a height map.

    Args:
        cloud: path to the output points cloud (ply format)
        heights: height map, sampled on the same grid as the crop_colorized
            image. In particular, its size is the same as crop_colorized.
        rpc: path to xml file containing RPC data for the current Pleiade image
        H (optional, default None): numpy array of size 3x3 defining the
            homography transforming the coordinates system of the original full
            size image into the coordinates system of the crop we are dealing
            with.
        crop_colorized (optional, default ''): path to a colorized crop of a
            Pleiades image
        off_{x,y} (optional, default None): coordinates of the point we want to
            use as origin in the local coordinate system of the computed cloud
        ascii_ply (optional, default false): boolean flag to tell if the output
            ply file should be encoded in plain text (ascii).
        utm_zone (optional, default None):
    """
    if not os.path.exists(crop_colorized):
        crop_colorized = ''
    hij = " ".join([str(x) for x in H.flatten()]) if H is not None else ""
    asc = "--ascii" if ascii_ply else ""
    nrm = "--with-normals" if with_normals else ""
    utm = "--utm-zone %s" % utm_zone if utm_zone else ""
    lbb = "--lon-m %s --lon-M %s --lat-m %s --lat-M %s" % llbbx if llbbx else ""
    if not (os.path.isfile(cloud) and cfg['skip_existing']):
        command = "colormesh %s %s %s %s -h \"%s\" %s %s %s %s" % (cloud,
                          heights, rpc, crop_colorized, hij, asc, nrm, utm, lbb)
        if off_x:
            command += " --offset_x %d" % off_x
        if off_y:
            command += " --offset_y %d" % off_y
        common.run(command)
