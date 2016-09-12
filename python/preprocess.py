# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import multiprocessing
import numpy as np

from config import cfg
from python import common
from python import masking
from python import pointing_accuracy
from python import visualisation


def minmax_color_on_tile(tile_info):
    """
    Compute min and max intensities on a given tile and save them to a txt file.

    Args:
        tile_info: dictionary containing all the information needed to process
            the tile
    """
    # read info
    img = cfg['images'][0]['img']
    coords = tile_info['coordinates']
    tile_dir = tile_info['directory']
    z = cfg['subsampling_factor']

    # output files
    crop_ref = os.path.join(tile_dir, 'roi_ref.tif')
    local_minmax = os.path.join(tile_dir, 'local_minmax.txt')

    # do the job
    if os.path.isfile(crop_ref) and cfg['skip_existing']:
        print 'roi_ref.tif for tile %s already generated, skip' % tile_dir
    else:
        common.cropImage(img, crop_ref, *coords, zoom=z)
    if os.path.isfile(os.path.join(tile_dir, 'local_minmax.txt')) and cfg['skip_existing']:
        print 'extrema intensities on tile %s already computed, skip' % tile_dir
    else:
        common.image_getminmax(crop_ref, local_minmax)


def pointing_correction(tile_info):
    """
    Compute the translations that correct the pointing errors on a tile.

    There is one correction per pair of images.

    Args:
        tile_info: dictionary containing all the information needed to process
            the tile
    """
    tile_dir = tile_info['directory']
    x, y, w, h = tile_info['coordinates']

    for i in xrange(1, len(cfg['images'])):
        paired_tile_dir = os.path.join(tile_dir, 'pair_%d' % i)
        img1, rpc1 = cfg['images'][0]['img'], cfg['images'][0]['rpc']
        img2, rpc2 = cfg['images'][i]['img'], cfg['images'][i]['rpc']
        roi_msk = cfg['images'][0]['roi']
        cld_msk = cfg['images'][0]['cld']
        wat_msk = cfg['images'][0]['wat']

        # output files
        pointing = '%s/pointing.txt' % paired_tile_dir
        center = '%s/center_keypts_sec.txt' % paired_tile_dir
        sift_matches = '%s/sift_matches.txt' % paired_tile_dir

        # check if the tile is already done
        if os.path.isfile('%s/pointing.txt' % paired_tile_dir) and cfg['skip_existing']:
            print "pointing correction on tile %d %d (pair %d) already done, skip" % (x, y, i)
        else:
            # correct pointing error
            # A is the correction matrix and m is the list of sift matches
            A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2,
                                                        x, y, w, h)
            if A is not None:
                np.savetxt(pointing, A, fmt='%6.3f')
            if m is not None:
                np.savetxt(sift_matches, m, fmt='%9.3f')
                np.savetxt(center, np.mean(m[:, 2:4], 0), fmt='%9.3f')
                if cfg['debug']:
                    png = '%s/sift_matches_plot.png' % paired_tile_dir
                    visualisation.plot_matches_pleiades(img1, img2, rpc1, rpc2,
                                                        m, x, y, w, h, png)
