# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import os
import ctypes
from ctypes import c_int, c_float, c_double, byref, POINTER
from numpy.ctypeslib import ndpointer
import numpy as np
from scipy import ndimage

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
            rpc (rpcm.RPCModel): rpc model
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


def disp_to_xyz(rpc1, rpc2, H1, H2, disp, mask, utm_zone, img_bbx=None, A=None):
    """
    Compute a height map from a disparity map, using RPC camera models.

    Args:
        rpc1, rpc2 (rpcm.RPCModel): camera models
        H1, H2 (arrays): 3x3 numpy arrays defining the rectifying homographies
        disp, mask (array): 2D arrays of shape (h, w) representing the diparity
            and mask maps
        utm_zone (int): desired UTM zone number (between 1 and 60) for the
            output xyz map
        img_bbx (4-tuple): col_min, col_max, row_min, row_max defining the
            unrectified image domain to process.
        A (array): 3x3 array with the pointing correction matrix for im2

    Returns:
        xyz: array of shape (h, w, 3) where each pixel contains the UTM
            easting, northing, and altitude of the triangulated point.
        err: array of shape (h, w) where each pixel contains the triangulation
            error
    """
    # copy rpc coefficients to an RPCStruct object
    rpc1_c_struct = RPCStruct(rpc1)
    rpc2_c_struct = RPCStruct(rpc2)

    # handle optional arguments
    if A is not None:  # apply pointing correction
        H2 = np.dot(H2, np.linalg.inv(A))

    if img_bbx is None:
        img_bbx = (-np.inf, np.inf, -np.inf, np.inf)

    # define the argument types of the disp_to_xyz function from disp_to_h.so
    h, w = disp.shape
    lib.disp_to_xyz.argtypes = (ndpointer(dtype=c_float, shape=(h, w, 3)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                c_int, c_int,
                                ndpointer(dtype=c_double, shape=(9,)),
                                ndpointer(dtype=c_double, shape=(9,)),
                                POINTER(RPCStruct), POINTER(RPCStruct), c_int,
                                ndpointer(dtype=c_float, shape=(4,)))


    # call the disp_to_xyz function from disp_to_h.so
    xyz =  np.zeros((h, w, 3), dtype='float32')
    err =  np.zeros((h, w), dtype='float32')
    dispx = disp.astype('float32')
    dispy = np.zeros((h, w), dtype='float32')
    msk = mask.astype('float32')
    lib.disp_to_xyz(xyz, err, dispx, dispy, msk, w, h,
                    H1.flatten(), H2.flatten(),
                    byref(rpc1_c_struct), byref(rpc2_c_struct), utm_zone,
                    np.asarray(img_bbx, dtype='float32'))

    return xyz, err


def count_3d_neighbors(xyz, r, p):
    """
    Count 3D neighbors of a gridded set of 3D points.

    Args:
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            UTM easting, northing, and altitude of a 3D point.
        r (float): filtering radius, in meters
        p (int): the filering window has size 2p + 1, in pixels

    Returns:
        array of shape (h, w) with the count of the number of 3D points located
        less than r meters from the current 3D point
    """
    h, w, d = xyz.shape
    assert(d == 3)

    # define the argument types of the count_3d_neighbors function from disp_to_h.so
    lib.count_3d_neighbors.argtypes = (ndpointer(dtype=c_int, shape=(h, w)),
                                       ndpointer(dtype=c_float, shape=(h, w, 3)),
                                       c_int, c_int, c_float, c_int)

    # call the count_3d_neighbors function from disp_to_h.so
    out = np.zeros((h, w), dtype='int32')
    lib.count_3d_neighbors(out, np.ascontiguousarray(xyz), w, h, r, p)

    return out


def remove_isolated_3d_points(xyz, r, p, n, q=1):
    """
    Discard (in place) isolated (groups of) points in a gridded set of 3D points

    Discarded points satisfy the following conditions:
    - they have less than n 3D neighbors in a ball of radius r meters;
    - all their neighboring points of the grid in a square window of size 2q+1
      that are closer than r meters are also discarded.

    Args:
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            UTM easting, northing, and altitude of a 3D point.
        r (float): filtering radius, in meters
        p (int): filering window radius, in pixels (square window of size 2p+1)
        n (int): filtering threshold, in number of points
        q (int): 2nd filtering window radius, in pixels (square of size 2q+1)
    """
    h, w, d = xyz.shape
    assert d == 3, 'expecting a 3-channels image with shape (h, w, 3)'

    lib.remove_isolated_3d_points.argtypes = (
        ndpointer(dtype=c_float, shape=(h, w, 3)),
        c_int, c_int, c_float, c_int, c_int, c_int)

    lib.remove_isolated_3d_points(np.ascontiguousarray(xyz), w, h, r, p, n, q)


def filter_xyz(xyz, r, n, img_gsd):
    """
    Discard (in place) points that have less than n points closer than r meters.

    Args:
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            UTM easting, northing, and altitude of a 3D point.
        r (float): filtering radius, in meters
        n (int): filtering threshold, in number of points
        img_gsd (float): ground sampling distance, in meters / pix
    """
    p = np.ceil(r / img_gsd).astype(int)
    remove_isolated_3d_points(xyz, r, p, n)


def height_map(x, y, w, h, rpc1, rpc2, H1, H2, disp, mask, utm_zone, A=None):
    """
    Computes an altitude map, on the grid of the original reference image, from
    a disparity map given on the grid of the rectified reference image.

    Args:
        x, y, w, h (ints): rectangular AOI in the original image. (x, y) is the
            top-left corner, and (w, h) are the dimensions of the rectangle.
        rpc1, rpc2 (rpcm.RPCModel): camera models
        H1, H2 (arrays): 3x3 numpy arrays defining the rectifying homographies
        disp, mask (array): 2D arrays of shape (h, w) representing the diparity
            and mask maps
        utm_zone (int): desired UTM zone number (between 1 and 60) for the
            output xyz map
        A (array): 3x3 array with the pointing correction matrix for im2

    Returns:
        array of shape (h, w) with the height map
    """
    xyz, err = disp_to_xyz(rpc1, rpc2, H1, H2, disp, mask, utm_zone, A=A)
    height_map = xyz[:, :, 2].squeeze()

    # transfer the rectified height map onto an unrectified height map
    H = np.dot(H1, common.matrix_translation(x, y))
    out = ndimage.affine_transform(np.nan_to_num(height_map).T, H,
                                   output_shape=(w, h), order=1).T

    # nearest-neighbor interpolation of nan locations in the resampled image
    if np.isnan(height_map).any():
        i = ndimage.affine_transform(np.isnan(height_map).T, H,
                                     output_shape=(w, h), order=0).T
        i = ndimage.binary_dilation(i, structure=np.ones((3, 3)))
        out[i] = np.nan  # put nans back in the resampled image

    return out


def height_map_to_point_cloud(cloud, heights, rpc, H=None, crop_colorized='',
                              off_x=None, off_y=None, ascii_ply=False,
                              with_normals=False, utm_zone=None, llbbx=None):
    """
    Computes a color point cloud from a height map.

    Args:
        cloud: path to the output points cloud (ply format)
        heights: height map, sampled on the same grid as the crop_colorized
            image. In particular, its size is the same as crop_colorized.
        rpc: instances of the rpcm.RPCModel class
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
