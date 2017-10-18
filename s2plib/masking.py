# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import subprocess
import numpy as np
from osgeo import gdal

from s2plib import common
from s2plib.config import cfg


def cloud_water_image_domain(x, y, w, h, rpc, roi_gml=None, cld_gml=None,
                             wat_msk=None, use_srtm_for_water=False):
    """
    Compute a mask for pixels masked by clouds, water, or out of image domain.

    Args:
        x, y, w, h: coordinates of the ROI
        rpc: path to the xml file containing the rpc coefficients of the image
            RPC model is used with SRTM data to derive the water mask
        roi_gml (optional): path to a gml file containing a mask
            defining the area contained in the full image
        cld_gml (optional): path to a gml file containing a mask
            defining the areas covered by clouds
        wat_msk (optional): path to an image file containing a water mask

    Returns:
        2D array containing the output binary mask. 0 indicate masked pixels, 1
        visible pixels.
    """
    # coefficients of the transformation associated to the crop
    H = common.matrix_translation(-x, -y)
    hij = ' '.join([str(el) for el in H.flatten()])
    
    mask = np.ones((h, w), dtype=np.bool)

    if roi_gml is not None:  # image domain mask (polygons)
        tmp = common.tmpfile('.png')
        subprocess.check_call('cldmask %d %d -h "%s" %s %s' % (w, h, hij,
                                                               roi_gml, tmp),
                              shell=True)

        f = gdal.Open(tmp)
        mask = np.logical_and(mask, f.ReadAsArray())
        f = None  # this is the gdal way of closing files

    if not mask.any():
        return mask

    if cld_gml is not None:  # cloud mask (polygons)
        tmp = common.tmpfile('.png')
        subprocess.check_call('cldmask %d %d -h "%s" %s %s' % (w, h, hij,
                                                               cld_gml, tmp),
                              shell=True)
        f = gdal.Open(tmp)
        mask = np.logical_and(mask, ~f.ReadAsArray().astype(bool))
        f = None  # this is the gdal way of closing files

    if not mask.any():
        return mask

    if wat_msk is not None:  # water mask (raster)
        f = gdal.Open(wat_msk)
        mask = np.logical_and(mask, f.ReadAsArray(x, y, w, h))
        f = None  # this is the gdal way of closing files

    elif use_srtm_for_water:  # water mask (srtm)
        tmp = common.tmpfile('.png')
        env = os.environ.copy()
        env['SRTM4_CACHE'] = cfg['srtm_dir']
        subprocess.check_call('watermask %d %d -h "%s" %s %s' % (w, h, hij, rpc,
                                                                 tmp),
                              shell=True, env=env)
        f = gdal.Open(tmp)
        mask = np.logical_and(mask, f.ReadAsArray())
        f = None  # this is the gdal way of closing files

    return mask


def erosion(out, msk, radius):
    """
    Erodes the accepted regions (ie eliminates more pixels)

    Args:
        out: path to the ouput mask image file
        msk: path to the input mask image file
        radius (in pixels): size of the disk used for the erosion
    """
    if radius >= 2:
        common.run('morsi disk%d erosion %s %s' % (int(radius), msk, out))
