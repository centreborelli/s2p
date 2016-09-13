# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import os
import numpy as np

from config import cfg
from python import pointing_accuracy
from python import visualisation


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
