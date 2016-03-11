#!/usr/bin/env python

# s2p - Satellite Stereo Pipeline
# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@polytechnique.org>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import time
import shutil
import os.path
import datetime
import numpy as np
import multiprocessing

from python.config import cfg
from python import common
from python import initialization
from python import preprocess
from python import globalvalues
from python import process
from python import globalfinalization


def show_progress(a):
    """
    Print the number of tiles that have been processed.

    Args:
        a: useless argument, but since this function is used as a callback by
            apply_async, it has to take one argument.
    """
    show_progress.counter += 1
    print 'done %d / %d tiles' % (show_progress.counter, show_progress.total)


def preprocess_tile(tile_info):
    """
    Compute pointing corrections and extrema intensities for a single tile.

    Args:
        tile_info: list containing all the informations needed to process a
            tile.
    """
    preprocess.pointing_correction(tile_info)
    preprocess.minmax_color_on_tile(tile_info)


def global_values(tiles_full_info):
    """
    Compute the global pointing correction and extrema intensities for the ROI.
    """
    globalvalues.pointing_correction(tiles_full_info)
    globalvalues.minmax_intensities(tiles_full_info)


def process_tile_pair(tile_info, pair_id):
    """
    Process a pair of images on a given tile.

    Processing includes rectification, disparity estimation and triangulation.

    Args:
        tile_info: list containing all the informations needed to process a
            tile.
        pair_id: index of the pair to process
    """
    # read all the information
    tile_dir = tile_info['directory']
    col, row, tw, th = tile_info['coordinates']
    images = cfg['images']
    cld_msk = cfg['images'][0]['cld']
    roi_msk = cfg['images'][0]['roi']

    img1, rpc1 = images[0]['img'], images[0]['rpc']
    img2, rpc2 = images[pair_id]['img'], images[pair_id]['rpc']

    out_dir = os.path.join(tile_dir, 'pair_%d' % pair_id)

    A_global = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_%d.txt' % pair_id)

    print 'processing tile %d %d...' % (col, row)
    # rectification
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'rectified_ref.tif')) and
        os.path.isfile(os.path.join(out_dir, 'rectified_sec.tif'))):
        print '\trectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\trectifying tile %d %d (pair %d)...' % (col, row, pair_id)
        process.rectify(out_dir, np.loadtxt(A_global), img1, rpc1,
                        img2, rpc2, col, row, tw, th, None, cld_msk,
                        roi_msk)

    # disparity estimation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'rectified_disp.tif'))):
        print '\tdisparity estimation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\testimating disparity on tile %d %d (pair %d)...' % (col, row, pair_id)
        process.disparity(out_dir, img1, rpc1, img2, rpc2, col, row,
                          tw, th, None, cld_msk, roi_msk)

    # triangulation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'height_map.tif'))):
        print '\ttriangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\ttriangulating tile %d %d (pair %d)...' % (col, row, pair_id)
        process.triangulate(out_dir, img1, rpc1, img2, rpc2, col,
                            row, tw, th, None, cld_msk, roi_msk,
                            np.loadtxt(A_global))


def process_tile(tile_info):
    """
    Process a tile by merging the height maps computed for each image pair.

    Args:
        tile_info: a dictionary that provides all you need to process a tile
    """
    tile_dir = tile_info['directory']

    # check that the tile is not masked
    if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
        print 'tile %s already masked, skip' % tile_dir
        return

    # process each pair to get a height map
    nb_pairs = tile_info['number_of_pairs']
    for pair_id in range(1, nb_pairs + 1):
        process_tile_pair(tile_info, pair_id)

    # finalization
    height_maps = [os.path.join(tile_dir, 'pair_%d' % i, 'height_map.tif') for i
                   in range(1, nb_pairs + 1)]
    process.finalize_tile(tile_info, height_maps)


def global_finalization(tiles_full_info):
    """
    Produce single height map and DSM for the whole region of interest.

    The height maps associated to each pair, as well as the height map obtained
    by merging all the pairs, are stored as VRT files. The final DSM is obtained
    by projecting the 3D points from the ply files obtained on each tile.

    Args:
        tiles_full_info: dictionary providing all the information about the
            processed tiles
    """
    globalfinalization.write_vrt_files(tiles_full_info)
    globalfinalization.write_dsm(tiles_full_info)

    # crop area corresponding to the ROI in the secondary images
    if not cfg['full_img']:
        common.crop_corresponding_areas(cfg['out_dir'], cfg['images'],
                                        cfg['roi'])

    # copy RPC xml files in the output directory
    for img in cfg['images']:
        shutil.copy2(img['rpc'], cfg['out_dir'])


def launch_parallel_calls(fun, list_of_args, nb_workers):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel. It must have a
            unique input argument.
        list_of_args: list of arguments passed to fun, one per call
        nb_workers: number of calls run simultaneously
    """
    results = []
    show_progress.counter = 0
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        p = pool.apply_async(fun, args=(x,), callback=show_progress)
        results.append(p)

    for r in results:
        try:
            r.get(3600)  # wait at most one hour per call
        except multiprocessing.TimeoutError:
            print "Timeout while running %s" % str(r)
        except common.RunFailure as e:
            print "FAILED call: ", e.args[0]["command"]
            print "\toutput: ", e.args[0]["output"]
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(1)


def main(config_file):
    """
    Launch the entire s2p pipeline with the parameters given by a json file.

    It is a succession of five steps:
        initialization
        preprocessing
        global_values
        processing
        global_finalization

    Args:
        config_file: path to a json configuration file
    """
    t0 = time.time()

    # initialization
    initialization.init_dirs_srtm_roi(config_file)
    tiles_full_info = initialization.init_tiles_full_info(config_file)
    show_progress.total = len(tiles_full_info)

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_nb_threads']:
        nb_workers = min(nb_workers, cfg['max_nb_threads'])

    # omp_num_threads: should not exceed nb_workers when multiplied by the
    # number of tiles
    cfg['omp_num_threads'] = max(1, int(nb_workers / len(tiles_full_info)))

    # do the job
    print '\npreprocessing tiles...'
    launch_parallel_calls(preprocess_tile, tiles_full_info, nb_workers)
    print "Elapsed time:", datetime.timedelta(seconds=int(time.time() - t0))

    print '\ncomputing global values...'
    global_values(tiles_full_info)
    print "Elapsed time:", datetime.timedelta(seconds=int(time.time() - t0))

    print '\nprocessing tiles...'
    launch_parallel_calls(process_tile, tiles_full_info, nb_workers)
    print "Elapsed time:", datetime.timedelta(seconds=int(time.time() - t0))

    print '\nglobal finalization...'
    global_finalization(tiles_full_info)
    print "Total runtime:", datetime.timedelta(seconds=int(time.time() - t0))

    # cleanup
    common.garbage_cleanup()


if __name__ == '__main__':

    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print """
        Incorrect syntax, use:
          > %s config.json

          Launches the s2p pipeline. All the parameters, paths to input and
          output files, are defined in the json configuration file.
        """ % sys.argv[0]
        sys.exit(1)
