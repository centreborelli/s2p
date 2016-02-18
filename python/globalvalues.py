# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import os 
import numpy as np

from config import cfg
from python import pointing_accuracy


def global_pointing_correction(tiles_full_info):
    """
    Computes the global pointing correction. For that purpose, global_pointing_correction needs all the tile related to one pair (given by pairedTilesPerPairId),
    and the total number of pairs as well.
    """

    # Build pairedTilesPerPairId : a dictionary that provides the list of all
    # the tile for a given pair_id (key)
    pairedTilesPerPairId = {}
    for tile_info in tiles_full_info:
        col, row, tw, th, ov, i, j, pos, x, y, w, h, images, NbPairs, cld_msk, roi_msk, tile_dir = tile_info

        for i in range(0, NbPairs):
            pair_id = i + 1
            pair_dir = tile_dir + 'pair_%d' % pair_id
            if pairedTilesPerPairId.has_key(pair_dir):
                pairedTilesPerPairId[pair_id].append(pair_dir)
            else:
                pairedTilesPerPairId[pair_id] = []
     # Build pairedTilesPerPairId (end)

    for pair_id in pairedTilesPerPairId:
        tiles = pairedTilesPerPairId[pair_id]
        A_globalMat = pointing_accuracy.global_from_local(tiles)
        np.savetxt('%s/global_pointing_pair_%d.txt' %
                   (cfg['out_dir'], pair_id), A_globalMat)


def global_minmax_intensities(tiles_full_info):
    """
    Computes the min and max intensities from the tiles that will be processed.
    This will allow to re-code colors by using 8-bits instead of 12-bits or more, and to better vizualise the ply files.

    Args:
         tiles_full_info: a list of tile_info dictionaries
    """
    minlist = []
    maxlist = []
    for tile_info in tiles_full_info:
        minmax = np.loadtxt(os.path.join(tile_info[-1], 'local_minmax.txt'))
        minlist.append(minmax[0])
        maxlist.append(minmax[1])

    global_minmax = [min(minlist), max(maxlist)]

    np.savetxt(os.path.join(cfg['out_dir'], 'global_minmax.txt'), global_minmax)
