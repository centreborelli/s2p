# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import os
import ctypes
from ctypes import c_int, c_float, c_double, byref, POINTER
from numpy.ctypeslib import ndpointer
import numpy as np
from scipy import ndimage
import rasterio

from s2p import common
from s2p.config import cfg
from s2p import ply
from s2p import geographiclib

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


def disp_to_xyz(rpc1, rpc2, H1, H2, disp, mask, out_crs=None, img_bbx=None, A=None):
    """
    Compute a height map from a disparity map, using RPC camera models.

    Args:
        rpc1, rpc2 (rpcm.RPCModel): camera models
        H1, H2 (arrays): 3x3 numpy arrays defining the rectifying homographies
        disp, mask (array): 2D arrays of shape (h, w) representing the diparity
            and mask maps
        out_crs (pyproj.crs.CRS): object defining the desired coordinate reference system for the
            output xyz map
        img_bbx (4-tuple): col_min, col_max, row_min, row_max defining the
            unrectified image domain to process.
        A (array): 3x3 array with the pointing correction matrix for im2

    Returns:
        lonlatalt: array of shape (h, w, 3) where each pixel contains the
            longitude, latitude and altitude of the triangulated point.
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

    # define the argument types of the disp_to_lonlatalt function from disp_to_h.so
    h, w = disp.shape
    lib.disp_to_lonlatalt.argtypes = (ndpointer(dtype=c_double, shape=(h, w, 3)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                ndpointer(dtype=c_float, shape=(h, w)),
                                c_int, c_int,
                                ndpointer(dtype=c_double, shape=(9,)),
                                ndpointer(dtype=c_double, shape=(9,)),
                                POINTER(RPCStruct), POINTER(RPCStruct),
                                ndpointer(dtype=c_float, shape=(4,)))


    # call the disp_to_lonlatalt function from disp_to_h.so
    lonlatalt =  np.zeros((h, w, 3), dtype='float64')
    err =  np.zeros((h, w), dtype='float32')
    dispx = disp.astype('float32')
    dispy = np.zeros((h, w), dtype='float32')
    msk = mask.astype('float32')
    lib.disp_to_lonlatalt(lonlatalt, err, dispx, dispy, msk, w, h,
                    H1.flatten(), H2.flatten(),
                    byref(rpc1_c_struct), byref(rpc2_c_struct),
                    np.asarray(img_bbx, dtype='float32'))

    # output CRS conversion
    in_crs = geographiclib.pyproj_crs("epsg:4979")

    if out_crs and out_crs != in_crs:

        # reshape the xyz array into a 3-column 2D-array
        lonlatalt_shape = lonlatalt.shape
        lonlatalt = lonlatalt.reshape(-1, 3)
        xyz_array = np.empty(lonlatalt_shape, dtype=np.float32)

        x, y, z = geographiclib.pyproj_transform(lonlatalt[:,0], lonlatalt[:,1],
                                             in_crs, out_crs, lonlatalt[:,2])

        xyz_array[:,:,0] = x.reshape(*lonlatalt_shape[:2])
        xyz_array[:,:,1] = y.reshape(*lonlatalt_shape[:2])
        xyz_array[:,:,2] = z.reshape(*lonlatalt_shape[:2])
    else:
        xyz_array = lonlatalt

    return xyz_array, err


def count_3d_neighbors(xyz, r, p):
    """
    Count 3D neighbors of a gridded set of 3D points.

    Args:
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            UTM easting, northing, and altitude of a 3D point.
        r (float): filtering radius, in the unit of the CRS (ex: meters)
        p (int): the filering window has size 2p + 1, in pixels

    Returns:
        array of shape (h, w) with the count of the number of 3D points located
        less than r units from the current 3D point
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
    - they have less than n 3D neighbors in a ball of radius r units (ex: meters);
    - all their neighboring points of the grid in a square window of size 2q+1
      that are closer than r units are also discarded.

    Args:
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            UTM easting, northing, and altitude of a 3D point.
        r (float): filtering radius, in the unit of the CRS (ex: meters)
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
    Discard (in place) points that have less than n points closer than r units (ex: meters).

    Args:
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            UTM easting, northing, and altitude of a 3D point.
        r (float): filtering radius, in the unit of the CRS (ex: meters)
        n (int): filtering threshold, in number of points
        img_gsd (float): ground sampling distance, in units of the CRS (ex: meters) / pix
    """
    p = np.ceil(r / img_gsd).astype(int)
    remove_isolated_3d_points(xyz, r, p, n)


def height_map(x, y, w, h, rpc1, rpc2, H1, H2, disp, mask, A=None):
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
        A (array): 3x3 array with the pointing correction matrix for im2

    Returns:
        array of shape (h, w) with the height map
    """
    xyz, err = disp_to_xyz(rpc1, rpc2, H1, H2, disp, mask, A=A)
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


def filter_xyz_and_write_to_ply(path_to_ply_file, xyz, r, n, img_gsd, colors='', proj_com=''):
    """
    Filter points that have less than n points closer than r units (ex: meters) and write them in a .ply file

    Args:
        path_to_ply_file (str): path to a .ply file
        xyz (array): 3D array of shape (h, w, 3) where each pixel contains the
            x, y, and z  coordinates of a 3D point.
        r (float): filtering radius, in the unit of the CRS (ex: meters)
        n (int): filtering threshold, in number of points
        img_gsd (float): ground sampling distance, in units of the CRS (ex: meters) / pix
        colors (optional, default ''): path to a colorized image
        proj_com (str): projection comment in the .ply file
    """
    # 3D filtering
    if r and n:
        filter_xyz(xyz, r, n, img_gsd)

    # flatten the xyz array into a list and remove nan points
    xyz_list = xyz.reshape(-1, 3)
    valid = np.all(np.isfinite(xyz_list), axis=1)

    # write the point cloud to a ply file
    if colors:
        with rasterio.open(colors, 'r') as f:
            img = f.read()
        colors_list = img.transpose(1, 2, 0).reshape(-1, img.shape[0])[valid]
    else:
        colors_list = None

    ply.write_3d_point_cloud_to_ply(path_to_ply_file, xyz_list[valid],
                                    colors=colors_list,
                                    extra_properties=None,
                                    extra_properties_names=None,
                                    comments=["created by S2P",
                                              "projection: {}".format(proj_com)])


def height_map_to_point_cloud(cloud, heights, rpc, off_x=None, off_y=None, crop_colorized=''):
    """
    Computes a color point cloud from a height map.

    Args:
        cloud: path to the output points cloud (ply format)
        heights: height map, sampled on the same grid as the crop_colorized
            image. In particular, its size is the same as crop_colorized.
        rpc: instances of the rpcm.RPCModel class
        off_{x,y} (optional, default None): coordinates of the origin of the crop
            we are dealing with in the pixel coordinates of the original full
            size image
        crop_colorized (optional, default ''): path to a colorized crop of a
            Pleiades image
    """
    with rasterio.open(heights) as src:
        h_map = src.read(1)
        map_shape = h_map.shape

    heights = h_map.ravel()
    indices = np.indices(map_shape)

    non_nan_ind = np.where(~np.isnan(heights))[0]

    alts = heights[non_nan_ind]
    cols = indices[1].ravel()[non_nan_ind]
    rows = indices[0].ravel()[non_nan_ind]

    if off_x or off_y:
        cols = cols + (off_x or 0)
        rows = rows + (off_y or 0)

    # localize pixels
    lons = np.empty_like(heights, dtype=np.float64)
    lats = np.empty_like(heights, dtype=np.float64)
    lons[non_nan_ind], lats[non_nan_ind] = rpc.localization(cols, rows, alts)

    # output CRS conversion
    in_crs = geographiclib.pyproj_crs("epsg:4979")
    pyproj_out_crs = geographiclib.pyproj_crs(cfg['out_crs'])
    proj_com = "CRS {}".format(cfg['out_crs'])

    if pyproj_out_crs != in_crs:
        x, y, z = geographiclib.pyproj_transform(lons, lats,
                                             in_crs, pyproj_out_crs, heights)
    else:
        x, y, z = lons, lats, heights

    xyz_array = np.zeros((*map_shape, 3))
    xyz_array[:,:,0] = x.reshape(map_shape)
    xyz_array[:,:,1] = y.reshape(map_shape)
    xyz_array[:,:,2] = z.reshape(map_shape)

    filter_xyz_and_write_to_ply(cloud, xyz_array,
                                              cfg['3d_filtering_r'], cfg['3d_filtering_n'],
                                              cfg['gsd'], crop_colorized, proj_com)
