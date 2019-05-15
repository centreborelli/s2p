# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import os
import numpy as np

from s2p import common
from s2p.config import cfg


def height_map_rectified(rpc1, rpc2, H1, H2, disp, mask, height, rpc_err, A=None):
    """
    Computes a height map from a disparity map, using rpc.

    Args:
        rpc1, rpc2: instances of the rpc_model.RPCModel class
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

    # write rpc coefficients to txt files
    rpcfile1 = common.tmpfile('.txt')
    rpcfile2 = common.tmpfile('.txt')
    rpc1.write_to_file(rpcfile1)
    rpc2.write_to_file(rpcfile2)

    common.run("disp_to_h %s %s %s %s %s %s %s %s" % (rpcfile1, rpcfile2, H1, HH2, disp,
                                                      mask, height, rpc_err))


def transfer_map(in_map, H, x, y, w, h, out_map):
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
        out_map: path to the output map
    """
    # write the inverse of the resampling transform matrix. In brief it is:
    # homography * translation
    # This matrix transports the coordinates of the original cropped and
    # grid (the one desired for out_height) to the rectified cropped and
    # grid (the one we have for height)
    HH = np.dot(np.loadtxt(H), common.matrix_translation(x, y))

    # apply the homography
    # write the 9 coefficients of the homography to a string, then call synflow
    # to produce the flow, then backflow to apply it
    # zero:256x256 is the iio way to create a 256x256 image filled with zeros
    hij = ' '.join(['%r' % num for num in HH.flatten()])
    common.run('synflow hom "%s" zero:%dx%d /dev/null - | BILINEAR=1 backflow - %s %s' % (
        hij, w, h, in_map, out_map))

    # replace the -inf with nan in the heights map, because colormesh filter
    # out nans but not infs
    # implements: if isinf(x) then nan, else x
    # common.run('plambda %s "x isinf nan x if" > %s' % (tmp_h, out_height))


def height_map(out, x, y, w, h, rpc1, rpc2, H1, H2, disp, mask, rpc_err,
               out_filt, A=None):
    """
    Computes an altitude map, on the grid of the original reference image, from
    a disparity map given on the grid of the rectified reference image.

    Args:
        out: path to the output file
        x, y, w, h: four integers defining the rectangular ROI in the original
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        rpc1, rpc2: instances of the rpc_model.RPCModel class
        H1, H2: path to txt files containing two 3x3 numpy arrays defining
            the rectifying homographies
        disp, mask: paths to the diparity and mask maps
        rpc_err: path to the output rpc_error of triangulation
        A (optional): path to txt file containing the pointing correction matrix
            for im2
    """
    tmp = common.tmpfile('.tif')
    height_map_rectified(rpc1, rpc2, H1, H2, disp, mask, tmp, rpc_err, A)
    transfer_map(tmp, H1, x, y, w, h, out)

    # apply output filter
    common.run('plambda {0} {1} "x 0 > y nan if" -o {1}'.format(out_filt, out))


def disp_map_to_point_cloud(out, disp, mask, rpc1, rpc2, H1, H2, A, colors, extra='',
                            utm_zone=None, llbbx=None, xybbx=None, xymsk=None):
    """
    Computes a 3D point cloud from a disparity map.

    Args:
        out: path to the output ply file
        disp, mask: paths to the diparity and mask maps
        rpc1, rpc2: instances of the rpc_model.RPCModel class
        H1, H2: path to txt files containing two 3x3 numpy arrays defining
            the rectifying homographies
        A: path to txt file containing the pointing correction matrix
            for im2
        colors: path to the png image containing the colors
    """
    href = " ".join(str(x) for x in np.loadtxt(H1).flatten())
    hsec = " ".join(str(x) for x in np.dot(np.loadtxt(H2),
                                           np.linalg.inv(np.loadtxt(A))).flatten())
    utm = "--utm-zone %s" % utm_zone if utm_zone else ""
    lbb = "--lon-m %s --lon-M %s --lat-m %s --lat-M %s" % llbbx if llbbx else ""
    xbb = "--col-m %s --col-M %s --row-m %s --row-M %s" % xybbx if xybbx else ""
    msk = "--mask-orig %s" % xymsk if xymsk else ""

    # write rpc coefficients to txt files
    rpcfile1 = common.tmpfile('.txt')
    rpcfile2 = common.tmpfile('.txt')
    rpc1.write_to_file(rpcfile1)
    rpc2.write_to_file(rpcfile2)

    # run disp2ply
    command = 'disp2ply {} {} {} {} {}'.format(out, disp, mask, rpcfile1, rpcfile2)
    # extra: is an optinonal extra data channel in the ply its default value '' ignores it
    command += ' {} {} -href "{}" -hsec "{}"'.format(colors, extra, href, hsec)
    command += ' {} {} {} {}'.format(utm, lbb, xbb, msk)
    common.run(command)


def multidisp_map_to_point_cloud(out, disp_list, rpc_ref, rpc_list, colors,
                                 utm_zone=None, llbbx=None, xybbx=None):
    """
    Computes a 3D point cloud from N disparity maps.

    Args:
        out: path to the output ply file
        disp_list: paths to the diparity maps
        rpc_ref, rpc_list: paths to the xml files
        colors: path to the png image containing the colors
    """

    disp_command = ['--disp%d %s' % (i+1, disp) for i, disp in enumerate(disp_list)]
    rpc_command = ['--rpc_sec%d %s' % (i+1, rpc) for i, rpc in enumerate(rpc_list)]

    utm = "--utm-zone %s" % utm_zone if utm_zone else ""
    lbb = "--lon-m %s --lon-M %s --lat-m %s --lat-M %s" % llbbx if llbbx else ""
    xbb = "--col-m %s --col-M %s --row-m %s --row-M %s" % xybbx if xybbx else ""

    command = 'multidisp2ply {} {} {} {} {}'.format(out, len(disp_list),
                                                    " ".join(disp_command),
                                                    "--rpc_ref %s" % rpc_ref,
                                                    " ".join(rpc_command))
    command += ' --color {}'.format(colors)
    command += ' {} {} {}'.format(utm, lbb, xbb)
    common.run(command)

def height_map_to_point_cloud(cloud, heights, rpc, H=None, crop_colorized='',
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
    hij = " ".join(str(x) for x in H.flatten()) if H is not None else ""
    asc = "--ascii" if ascii_ply else ""
    nrm = "--with-normals" if with_normals else ""
    utm = "--utm-zone %s" % utm_zone if utm_zone else ""
    lbb = "--lon-m %s --lon-M %s --lat-m %s --lat-M %s" % llbbx if llbbx else ""
    command = "colormesh %s %s %s %s -h \"%s\" %s %s %s %s" % (cloud, heights,
                                                               rpc,
                                                               crop_colorized,
                                                               hij, asc, nrm,
                                                               utm, lbb)
    if off_x:
        command += " --offset_x %d" % off_x
    if off_y:
        command += " --offset_y %d" % off_y
    common.run(command)
