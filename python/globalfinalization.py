# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import numpy as np
from config import cfg
import os
import shutil

from python import tile_composer
from python import common


def write_vrt_files(tiles_full_info):
    """
    Merges pieces of data into single VRT files : height map comprising the N pairs, height map for each signle pair, and err_rpc 

    Args:
         tiles_full_info: a list of tile_info dictionaries
    """
    # VRT file : height map (N pairs)

    #-------------- tileComposerInfo --------------
    # tileComposerInfo : a simple list that gives the size of the tiles (after having removed the overlapping areas),
    # regarding their positions inside the ROI --> ULw , ULh, Lw , Lh, Uw ,
    # Uh, Mw , Mh
    x0, y0, tw, th = tiles_full_info[0]['coordinates']
    x, y, w, h = tiles_full_info[0]['roi_coordinates']
    ov = tiles_full_info[0]['overlap']
    nb_pairs = tiles_full_info[0]['number_of_pairs']

    ULw, ULh = tw - ov / 2, th - ov / 2  # Size of Corner tiles (UL UR BL BR)
    Lw, Lh = tw - ov / 2, th - ov  # Size of Left/Right tile (L R)
    Uw, Uh = tw - ov, th - ov / 2  # Size of Upper/Bottom tile (U B)
    Mw, Mh = tw - ov, th - ov  # Size of Middle tile

    rangey = np.arange(y, y + h - ov, th - ov)
    rangex = np.arange(x, x + w - ov, tw - ov)
    rowmin, rowmax = rangey[0], rangey[-1]
    colmin, colmax = rangex[0], rangex[-1]

    # Tile composer info (bis)
    imax = len(rangey) - 1
    jmax = len(rangex) - 1
    # height of the final height map after removing margins
    fh = 2 * ULh + (imax - 1) * Mh
    # width of the final height map after removing margins
    fw = 2 * ULw + (jmax - 1) * Mh
    #-------------- tileComposerInfo (end) --------------

    tileSizesAndPositions = {}
    for tile_info in tiles_full_info:
        col, row, tw, th = tile_info['coordinates']
        pos = tile_info['position_type']
        i, j = tile_info['index_in_roi']

        dicoPos = {}
        dicoPos['M'] = [ULw + (j - 1) * Mw, ULh + (i - 1) * Mh, Mw, Mh]
        dicoPos['UL'] = [0, 0, ULw, ULh]
        dicoPos['UR'] = [fw - ULw, 0, ULw, ULh]
        dicoPos['BR'] = [fw - ULw, fh - ULh, ULw, ULh]
        dicoPos['BL'] = [0, fh - ULh, ULw, ULh]
        dicoPos['L'] = [0, ULh + (i - 1) * Lh, Lw, Lh]
        dicoPos['R'] = [fw - ULw, ULh + (i - 1) * Lh, Lw, Lh]
        dicoPos['U'] = [ULw + (j - 1) * Uw, 0, Uw, Uh]
        dicoPos['B'] = [ULw + (j - 1) * Uw, fh - ULh, Uw, Uh]
        dicoPos['Single'] = [0, 0, tw, th]

        tile_reldir = 'tile_%d_%d_row_%d/col_%d/' % (tw, th, row, col)

        tileSizesAndPositions[tile_reldir] = dicoPos[pos]

    z = cfg['subsampling_factor']
    tile_composer.mosaic_gdal2(cfg['out_dir'] + '/heightMap_N_pairs.vrt',
                               tileSizesAndPositions, 'local_merged_height_map_crop.tif', fw, fh, z)

    # VRT file : height map (for each single pair)
    # VRT file : rpc_err (for each single pair)
    for i in range(0, nb_pairs):

        pair_id = i + 1

        pairSizesAndPositions = {}
        for tile_reldir in tileSizesAndPositions:
            pair_reldir = tile_reldir + '/pair_%d' % (pair_id)
            pairSizesAndPositions[
                pair_reldir] = tileSizesAndPositions[tile_reldir]

        tile_composer.mosaic_gdal2(cfg['out_dir'] + '/heightMap_pair_%d.vrt' % (
            pair_id), pairSizesAndPositions, 'height_map_crop.tif', fw, fh, z)
        tile_composer.mosaic_gdal2(cfg['out_dir'] + '/rpc_err_pair_%d.vrt' % (
            pair_id), pairSizesAndPositions, 'rpc_err_crop.tif', fw, fh, z)


def write_dsm(tiles_full_info, n=5):
    """
    Writes the DSM, from the ply files given by each tile.

    Args :
         tiles_full_info: a list of tile_info dictionaries
    """
    clouds_dir = os.path.join(cfg['out_dir'], 'clouds')
    if (os.path.exists(clouds_dir)):
        shutil.rmtree(clouds_dir)
    os.mkdir(clouds_dir)

    for tile_info in tiles_full_info:
        tile_dir = tile_info['directory']
        x, y, w, h = tile_info['coordinates']
        cloud = os.path.join(os.path.abspath(tile_dir), 'cloud.ply')
        cloud_link_name = os.path.join(clouds_dir,
                                       'cloud_%d_%d_row_%d_col_%d.ply' % (w, h,
                                                                          x, y))
        if (os.path.exists(cloud)):
            common.run('ln -s %s %s' % (cloud, cloud_link_name))
    out_dsm = os.path.join(cfg['out_dir'], 'dsm.tif')

    common.run("ls %s | plyflatten %f %s" % (os.path.join(clouds_dir, 'cloud*'),
                                             cfg['dsm_resolution'], out_dsm))
    #common.run("gdalbuildvrt %s %s" %
    #           (cfg['out_dir'] + '/dsm.vrt', out_dsm_dir + '/dsm*'))
