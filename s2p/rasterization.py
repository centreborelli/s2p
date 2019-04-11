# Copyright (C) 2019, David Youssefi (CNES) <david.youssefi@cnes.fr>

from __future__ import print_function

import os
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer

from plyfile import PlyData, PlyElement
import affine
import pyproj

# TODO: This is kind of ugly. Cleaner way to do this is to update
# LD_LIBRARY_PATH, which we should do once we have a proper config file
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lib = ctypes.CDLL(os.path.join(parent_dir, 'lib', 'libplyflatten.so'))

def plyflatten(cloud,
               xoff, yoff,
               resolution,
               xsize, ysize,
               radius, sigma):
    """
    Projects a points cloud into the raster band(s) of a raster image

    Args:
        cloud: A nb_points x (2+nb_extra_columns) numpy array:
            | x0 y0 [z0 r0 g0 b0 ...] |
            | x1 y1 [z1 r1 g1 b1 ...] |
            | ...                     |
            | xN yN [zN rN gN bN ...] |
            x, y give positions of the points into the final raster, the "extra
            columns" give the values
        xoff, yoff: offset position (upper left corner) considering the georeferenced image
        resolution: resolution of the output georeferenced image
        xsize, ysize: size of the georeferenced image
        radius: controls the spread of the blob from each point
        sigma: radius of influence for each point (unit: pixel)

    Returns;
        A numpy array of shape (ysize, xsize, nb_extra_columns)
    """
    # Set expected args and return types
    lib.rasterize_cloud.argtypes = (ndpointer(dtype=ctypes.c_double,
                                              shape=np.shape(cloud)),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_double, ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_int, ctypes.c_int,
                                    ctypes.c_int, ctypes.c_float)
    lib.rasterize_cloud.restype = ctypes.POINTER(ctypes.c_float)

    # Call rasterize_cloud function from libplyflatten.so
    nb_points, nb_extra_columns = cloud.shape[0], cloud.shape[1] - 2
    raster_ptr = lib.rasterize_cloud(cloud.astype(np.float64),
                                     nb_points,
                                     nb_extra_columns,
                                     xoff, yoff,
                                     resolution,
                                     xsize, ysize,
                                     radius, sigma)

    # Transform result into a numpy array
    raster = np.asarray([raster_ptr[i] for i in range(nb_extra_columns*xsize*ysize)])
    raster = raster.reshape((ysize, xsize, nb_extra_columns))

    # Delete results to release memory
    lib.delete_buffer.argtypes = (ctypes.POINTER(ctypes.c_float)),
    lib.delete_buffer(raster_ptr)

    return raster


def plyflatten_from_plyfiles_list(clouds_list, resolution, radius=0, roi=None, sigma=None):
    """
    Projects a points cloud into the raster band(s) of a raster image (points clouds as files)

    Args:
        clouds_list: list of cloud.ply files
        resolution: resolution of the georeferenced output raster file
        roi: region of interest: (xoff, yoff, xsize, ysize), compute plyextrema if None

    Returns:
        raster: georeferenced raster
        profile: profile for rasterio
    """
    # read points clouds
    full_cloud = list()
    for cloud in clouds_list:
        plydata = PlyData.read(cloud)
        cloud_data = np.array(plydata.elements[0].data)
        proj = "projection:"
        utm_zone = [comment.split(proj)[-1] for comment in plydata.comments \
                    if proj in comment][0].split()[-1]

        # nb_extra_columns: z, r, g, b (all columns except x, y)
        nb_extra_columns = len(cloud_data.dtype) - 2
        full_cloud += [np.array([cloud_data[el] for el in cloud_data.dtype.names]).astype(np.float64).T]

    full_cloud = np.concatenate(full_cloud)
    nb_points = np.shape(full_cloud)[0]

    # region of interest (compute plyextrema if roi is None)
    if roi is not None:
        xoff, yoff, xsize, ysize = roi
    else:
        xx = full_cloud[:, 0]
        yy = full_cloud[:, 1]
        xmin = np.amin(xx)
        xmax = np.amax(xx)
        ymin = np.amin(yy)
        ymax = np.amax(yy)

        xsize = int(1 + np.floor((xmax - xmin) / resolution))
        ysize = int(1 + np.floor((ymax - ymin) / resolution))
        xoff = (xmax + xmin - resolution * xsize) / 2
        yoff = (ymax + ymin + resolution * ysize) / 2

    # The copy() method will reorder to C-contiguous order by default:
    full_cloud = full_cloud.copy()
    sigma = float("inf")  if sigma is None else sigma
    raster = plyflatten(full_cloud, xoff, yoff, resolution,
                        xsize, ysize,
                        radius, sigma)

    utm = pyproj.Proj(proj='utm', zone=utm_zone[:-1], ellps='WGS84', datum='WGS84',
                      south=(utm_zone[-1]=='S'))

    # construct profile dict
    profile = dict()
    profile['crs'] = utm.srs
    profile['transform'] = affine.Affine(resolution, 0.0, xoff,
                                         0.0, -resolution, yoff)

    return raster, profile
