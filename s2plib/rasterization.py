# Copyright (C) 2019, David Youssefi (CNES) <david.youssefi@cnes.fr>

from __future__ import print_function

import os
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer

# TODO: This is kind of ugly. Cleaner way to do this is to update
# LD_LIBRARY_PATH, which we should do once we have a proper config file
plyflatten_library = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), '../lib/libplyflatten.so')
ctypes.CDLL(plyflatten_library)

def plyflatten(cloud,
               nb_points,
               nb_extra_columns,
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
        nb_points: number of points to project
        nb_extra_columns: number of columns (all columns except x, y): x, y give
        positions of the points into the final raster, the "extra columns" give the
        values
        xoff, yoff: offset position (upper left corner) considering the georeferenced image
        resolution: resolution of the output georeferenced image
        xsize, ysize: size of the georeferenced image
        radius: controls the spread of the blob from each point
        sigma: radius of influence for each point (unit: pixel)

    Returns;
        A numpy array of shape (ysize, xsize, nb_extra_columns)
    """

    # load shared library
    lib = ctypes.CDLL(plyflatten_library)
    
    # Set expected args and return types
    lib.rasterize_cloud.argtypes = (ndpointer(dtype=ctypes.c_double,
                                              shape=np.shape(cloud)),
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_float, ctypes.c_float,
                                    ctypes.c_float,
                                    ctypes.c_int, ctypes.c_int,
                                    ctypes.c_int, ctypes.c_float)
    lib.rasterize_cloud.restype = ctypes.POINTER(ctypes.c_float)

    # Call rasterize_cloud function from libplyflatten.so
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
