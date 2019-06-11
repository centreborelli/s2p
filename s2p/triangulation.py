# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import os
import ctypes
from ctypes import c_int, c_float, c_double, byref, POINTER
from numpy.ctypeslib import ndpointer
import numpy as np

from s2p import common
from s2p.config import cfg


here = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(os.path.dirname(here), 'lib', 'disp_to_h.so')
lib = ctypes.CDLL(lib_path)


class RPCStruct(ctypes.Structure):
    """
    ctypes version of the RPC C struct defined in rpc.h.
    """
    _fields_ = [("numx", c_double * 20),
                ("denx", c_double * 20),
                ("numy", c_double * 20),
                ("deny", c_double * 20),
                ("scale", c_double * 3),
                ("offset", c_double * 3),
                ("inumx", c_double * 20),
                ("idenx", c_double * 20),
                ("inumy", c_double * 20),
                ("ideny", c_double * 20),
                ("iscale", c_double * 3),
                ("ioffset", c_double * 3),
                ("dmval", c_double * 4),
                ("imval", c_double * 4)]

    def __init__(self, rpc):
        """
        Args:
            rpc (rpc_model.RPCModel): rpc model
        """
        self.offset[0] = rpc.col_offset
        self.offset[1] = rpc.row_offset
        self.offset[2] = rpc.alt_offset
        self.ioffset[0] = rpc.lon_offset
        self.ioffset[1] = rpc.lat_offset
        self.ioffset[2] = rpc.alt_offset

        self.scale[0] = rpc.col_scale
        self.scale[1] = rpc.row_scale
        self.scale[2] = rpc.alt_scale
        self.iscale[0] = rpc.lon_scale
        self.iscale[1] = rpc.lat_scale
        self.iscale[2] = rpc.alt_scale

        for i in range(20):
            self.inumx[i] = rpc.col_num[i]
            self.idenx[i] = rpc.col_den[i]
            self.inumy[i] = rpc.row_num[i]
            self.ideny[i] = rpc.row_den[i]

        if hasattr(rpc, 'lat_num'):
            for i in range(20):
                self.numx[i] = rpc.lon_num[i]
                self.denx[i] = rpc.lon_den[i]
                self.numy[i] = rpc.lat_num[i]
                self.deny[i] = rpc.lat_den[i]
        else:
            for i in range(20):
                self.numx[i] = np.nan
                self.denx[i] = np.nan
                self.numy[i] = np.nan
                self.deny[i] = np.nan


def disp_to_xyz(rpc1, rpc2, H1, H2, disp, mask, utm_zone, A=None):
    """
    Compute a height map from a disparity map, using RPC camera models.

    Args:
        rpc1, rpc2 (rpc_model.RPCModel): camera models
        H1, H2 (arrays): 3x3 numpy arrays defining the rectifying homographies
        disp, mask (array): 2D arrays of shape (h, w) representing the diparity
            and mask maps
        utm_zone (int): desired UTM zone number (between 1 and 60) for the
            output xyz map
        A (array): 3x3 array with the pointing correction matrix for im2

    Returns:
        xyz: array of shape (h, w, 3) where each pixel contains the UTM
            easting, northing, and altitude of the triangulated point.
        err: array of shape (h, w) where each pixel contains the triangulation
            error
    """
    if A is not None:  # apply pointing correction
        H2 = np.dot(H2, np.linalg.inv(A))

    # copy rpc coefficients to an RPCStruct object
    rpc1_c_struct = RPCStruct(rpc1)
    rpc2_c_struct = RPCStruct(rpc2)

    # define the argument types of the disp_to_xyz function disp_to_h.so
    h, w = disp.shape
    lib.disp_to_xyz.argtypes = (ndpointer(dtype=c_float, shape=(h, w, 3)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                c_int, c_int,
                                ndpointer(dtype=c_double, shape=(3, 3)),
                                ndpointer(dtype=c_double, shape=(3, 3)),
                                POINTER(RPCStruct), POINTER(RPCStruct), c_int)


    # call the disp_to_xyz fonction from disp_to_h.so
    xyz =  np.zeros((h, w, 3), dtype='float32')
    err =  np.zeros((h, w), dtype='float32')
    dispx = disp.astype('float32')
    dispy = np.zeros((h, w), dtype='float32')
    msk = mask.astype('float32')
    lib.disp_to_xyz(xyz, err, dispx, dispy, msk, w, h, H1, H2,
                    byref(rpc1_c_struct), byref(rpc2_c_struct), utm_zone)

    return xyz, err


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
    xyz, err = disp_to_xyz(rpc1, rpc2, H1, H2, disp, mask, A)
    common.image_apply_homography()
    transfer_map(xyz, H1, x, y, w, h, out)

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
    # extra is an optional additional channel in the ply. Its default value '' ignores it
    command += ' {} {} -href "{}" -hsec "{}"'.format(colors, extra, href, hsec)
    command += ' {} {} {} {}'.format(utm, lbb, xbb, msk)
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
        rpc: instances of the rpc_model.RPCModel class
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
    # write rpc coefficients to txt file
    rpcfile = common.tmpfile('.txt')
    rpc.write_to_file(rpcfile)

    if not os.path.exists(crop_colorized):
        crop_colorized = ''
    hij = " ".join(str(x) for x in H.flatten()) if H is not None else ""
    asc = "--ascii" if ascii_ply else ""
    nrm = "--with-normals" if with_normals else ""
    utm = "--utm-zone %s" % utm_zone if utm_zone else ""
    lbb = "--lon-m %s --lon-M %s --lat-m %s --lat-M %s" % llbbx if llbbx else ""
    command = "colormesh %s %s %s %s -h \"%s\" %s %s %s %s" % (cloud, heights,
                                                               rpcfile,
                                                               crop_colorized,
                                                               hij, asc, nrm,
                                                               utm, lbb)
    if off_x:
        command += " --offset_x %d" % off_x
    if off_y:
        command += " --offset_y %d" % off_y
    common.run(command)
