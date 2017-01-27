# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

from __future__ import print_function
import subprocess
import zipfile
import os

from s2plib import common
from s2plib import rpc_model
from s2plib import rpc_utils
from s2plib.config import cfg


def list_srtm_tiles(rpcfile, x, y, w, h):
    """
    Tells which srtm tiles are needed to cover a given region of interest.

    Args:
        rpcfile: path to the xml file describing the rpc model.
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are
            the dimensions of the rectangle.

    Returns:
        the set of srtm tiles corresponding to the input ROI.
    """
    rpc = rpc_model.RPCModel(rpcfile)
    lon_min, lon_max, lat_min, lat_max = rpc_utils.geodesic_bounding_box(rpc,
                                                                         x, y,
                                                                         w, h)
    out = []
    for lon in [lon_min, lon_max]:
        for lat in [lat_min, lat_max]:
            p = subprocess.Popen(['srtm4_which_tile', str(lon), str(lat)],
                                 stdout=subprocess.PIPE)
            out.append(p.stdout.readline().split()[0].decode())
    out = set(out)
    print("Needed srtm tiles: ", out)
    return out


def get_srtm_tile(srtm_tile, out_dir):
    """
    Downloads and extract an srtm tile from the internet.

    Args:
        srtm_tile: string following the pattern 'srtm_%02d_%02d', identifying
            the desired strm tile
        out_dir: directory where to store and extract the srtm tiles
    """
    # check if the tile is already there
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if os.path.exists(os.path.join(out_dir, '%s.tif' % srtm_tile)):
        return

    # download the zip file
    srtm_tile_url = '%s/%s.zip' % (cfg['srtm_url'], srtm_tile)
    zip_path = os.path.join(out_dir, '%s.zip' % srtm_tile)
    common.download(zip_path, srtm_tile_url)

    # extract the tif file
    if zipfile.is_zipfile(zip_path):
        z = zipfile.ZipFile(zip_path, 'r')
        z.extract('%s.tif' % srtm_tile, out_dir)
    else:
        print("%s not available" % srtm_tile)

    # remove the zip file
    os.remove(zip_path)


def srtm4(lon, lat):
    """
    Gives the SRTM height of a point. It is a wrapper to the srtm4 binary.

    Args:
        lon, lat: longitude and latitude

    Returns:
        the height, in meters above the WGS84 geoid (not ellipsoid)
    """

    new_env = os.environ.copy()
    new_env['SRTM4_CACHE'] = cfg['srtm_dir']
    p = subprocess.Popen(['srtm4', str(lon), str(lat)],
                         stdout=subprocess.PIPE, env=new_env)
    return float(p.stdout.readline().split()[0])
