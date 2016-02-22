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
    img1 = cfg['images'][0]['img']
    coords = tile_info['coordinates']
    tile_dir = tile_info['directory']
    z = cfg['subsampling_factor']

    # output files
    crop_ref = os.path.join(tile_dir, 'roi_ref.tif')
    local_minmax = os.path.join(tile_dir, 'local_minmax.txt')

    # do the job
    common.cropImage(img1, crop_ref, *coords, zoom=z)
    if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
        print 'tile %s is masked, skip' % tile_dir
    elif os.path.isfile(os.path.join(tile_dir, 'local_minmax.txt')) and cfg['skip_existing']:
        print 'extrema intensities on tile %s already computed, skip' % tile_dir
    else:
        common.image_getminmax(crop_ref, local_minmax)


def pointing_correction(tile_info):
    """
    Computes pointing corrections

    Args:
        tile_info: dictionary containing all the information needed to process
            the tile
    """
    tile_dir = tile_info['directory']
    x, y, w, h = tile_info['coordinates']

    for i in range(1, tile_info['number_of_pairs'] + 1):
        paired_tile_dir = os.path.join(tile_dir, 'pair_%d' % i)
        img1, rpc1 = cfg['images'][0]['img'], cfg['images'][0]['rpc']
        img2, rpc2 = cfg['images'][i]['img'], cfg['images'][i]['rpc']
        roi_msk = cfg['images'][0]['roi']
        cld_msk = cfg['images'][0]['cld']

        # create a directory for the experiment
        if not os.path.exists(paired_tile_dir):
            os.makedirs(paired_tile_dir)

        # redirect stdout and stderr to log file
        if not cfg['debug']:
            # the last arg '0' is for no buffering
            fout = open(os.path.join(paired_tile_dir, 'stdout.log'), 'w', 0)
            sys.stdout = fout
            sys.stderr = fout

        # debug print
        print 'tile %d %d running on process %s' % (x, y, multiprocessing.current_process())

        # output files
        cwid_msk = '%s/cloud_water_image_domain_mask.png' % (paired_tile_dir)
        pointing = '%s/pointing.txt' % paired_tile_dir
        center = '%s/center_keypts_sec.txt' % paired_tile_dir
        sift_matches = '%s/sift_matches.txt' % paired_tile_dir
        sift_matches_plot = '%s/sift_matches_plot.png' % paired_tile_dir

        # check if the tile is already done
        if os.path.isfile('%s/pointing.txt' % paired_tile_dir) and cfg['skip_existing']:
            print "pointing correction on tile %d %d (pair %d) already done, skip" % (x, y, i)
        else:

            # check if the ROI is completely masked (water, or outside the
            # image domain)
            H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
            if masking.cloud_water_image_domain(cwid_msk, w, h, H, rpc1, roi_msk, cld_msk):
                print "Tile masked by water or outside definition domain, skip"
                open("%s/this_tile_is_masked.txt" % paired_tile_dir, 'a').close()
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                if not cfg['debug']:
                    fout.close()
            else:

                # correct pointing error
                # A is the correction matrix and m is the list of sift matches
                A, m = pointing_accuracy.compute_correction(img1, rpc1, img2,
                                                            rpc2, x, y, w, h)
                if A is not None:
                    np.savetxt(pointing, A)
                if m is not None:
                    np.savetxt(sift_matches, m)
                    np.savetxt(center, np.mean(m[:, 2:4], 0))
                    if cfg['debug']:
                        visualisation.plot_matches_pleiades(img1, img2, rpc1,
                                                            rpc2, m, x, y, w, h,
                                                            sift_matches_plot)

        # close logs
        common.garbage_cleanup()
        if not cfg['debug']:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            fout.close()

    return
