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
    show_progress.counter += 1
    if show_progress.counter > 1:
        print "Processed %d tiles" % show_progress.counter
    else:
        print "Processed 1 tile"


def preprocess_tile(tile_info):
    """
    1) Computes pointing corrections,
    2) Crops ref image into tiles and get min/max intensity values for each one (useful for colorizing ply files with 8-bits colors)
    """
    preprocess.pointing_correction(tile_info)
    preprocess.getMinMaxFromExtract(tile_info)


def global_values(tiles_full_info):
    """
    1) Computes the global pointing correction
    2) Computes the global min and max intensities from the tiles that will be processed
    """
    globalvalues.global_pointing_correction(tiles_full_info)
    globalvalues.global_minmax_intensities(tiles_full_info)


def process_tile_pair(tile_dir, pair_id, tiles_full_info):
    """
    Process a pair of images on a given tile.

    Processing includes rectification, disparity estimation and triangulation.

    Args:
        tile_dir: directory of the tile
        pair_id: index of the pair to process
        tiles_full_info: a dictionary providing all you need to process a tile
    """
    # read all the information
    col, row, tw, th = tiles_full_info[tile_dir][:4]
    images = tiles_full_info[tile_dir][12]
    cld_msk = tiles_full_info[tile_dir][14]
    roi_msk = tiles_full_info[tile_dir][15]

    img1, rpc1 = images[0]['img'], images[0]['rpc']
    img2, rpc2 = images[pair_id]['img'], images[pair_id]['rpc']

    paired_tile_dir = os.path.join(tile_dir, 'pair_%d' % pair_id)

    A_global = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_%d.txt' % pair_id)

    if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
        print 'tile %s already masked, skip' % paired_tile_dir
        return

    else:
        print 'processing tile %d %d...' % (col, row)

        # rectification
        if cfg['skip_existing'] and os.path.isfile(os.path.join(paired_tile_dir, 'rectified_ref.tif')) and os.path.isfile(os.path.join(paired_tile_dir, 'rectified_sec.tif')):
            print '\trectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print '\trectifying tile %d %d (pair %d)...' % (col, row, pair_id)
            process.rectify(paired_tile_dir, np.loadtxt(A_global), img1, rpc1,
                            img2, rpc2, col, row, tw, th, None, cld_msk,
                            roi_msk)

        # disparity estimation
        if cfg['skip_existing'] and os.path.isfile(os.path.join(paired_tile_dir,
                                                                'rectified_disp.tif')):
            print '\tdisparity estimation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print '\testimating disparity on tile %d %d (pair %d)...' % (col, row, pair_id)
            process.disparity(paired_tile_dir, img1, rpc1, img2, rpc2, col, row,
                              tw, th, None, cld_msk, roi_msk)

        # triangulation
        if cfg['skip_existing'] and os.path.isfile(os.path.join(paired_tile_dir,
                                                                'height_map.tif')):
            print '\ttriangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print '\ttriangulating tile %d %d (pair %d)...' % (col, row, pair_id)
            process.triangulate(paired_tile_dir, img1, rpc1, img2, rpc2, col,
                                row, tw, th, None, cld_msk, roi_msk,
                                np.loadtxt(A_global))

        return os.path.join(paired_tile_dir, 'height_map.tif')


def process_tile(tile_dir, tiles_full_info):
    """
    Process a tile by merging the height maps computed for each image pair.

    Args:
        tile_dir: directory of the tile to be processed
        tiles_full_info: a dictionary providing all you need to process a tile
    """
    nb_pairs = tiles_full_info[tile_dir][13]
    height_maps = []

    # process each pair to get a height map
    for pair_id in range(1, nb_pairs + 1):
        height_map = process_tile_pair(tile_dir, pair_id, tiles_full_info)
        if height_map:
            height_maps.append(height_map)

    # finalization
    process.finalize_tile(tile_dir, height_maps, tiles_full_info)


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


def map_processing(config_file):
    """
    Runs the entire s2p pipeline.

    It is a succession of five steps:
        initialization
        preprocessing
        global_values
        processing
        global_finalization

    Args:
        config_file: path to a json config file
    """
    try:
        # initialization
        initialization.init_dirs_srtm_roi(config_file)
        tiles_full_info = initialization.init_tiles_full_info(config_file)

        if cfg['debug']:  # monoprocessing

            print 'preprocess_tile...'
            for tile_info in tiles_full_info.values():
                preprocess_tile(tile_info)

            print 'global values...'
            global_values(tiles_full_info)

            print 'process_tile...'
            for tile_dir in tiles_full_info:
                process_tile(tile_dir, tiles_full_info)

            print 'global finalization...'
            global_finalization(tiles_full_info)

        else:  # multiprocessing

            nb_workers = multiprocessing.cpu_count()
            if cfg['max_nb_threads']:
                # use less workers than available cores
                nb_workers = min(nb_workers, cfg['max_nb_threads'])

            print 'preprocess_tile...'
            results = []
            show_progress.counter = 0
            pool = multiprocessing.Pool(nb_workers)
            for tile_info in tiles_full_info.values():
                p = pool.apply_async(preprocess_tile, args=(tile_info,),
                                     callback=show_progress)
                results.append(p)

            for r in results:
                try:
                    r.get(3600)  # wait at most one hour per tile
                except multiprocessing.TimeoutError:
                    print "Timeout while computing tile " + str(r)

            print 'global values...'
            global_values(tiles_full_info)

            print 'process_tile...'
            results = []
            show_progress.counter = 0
            pool = multiprocessing.Pool(nb_workers)
            for tile_dir in tiles_full_info:
                p = pool.apply_async(process_tile, args=(tile_dir,
                                                         tiles_full_info),
                                     callback=show_progress)
                results.append(p)

            for r in results:
                try:
                    r.get(3600)  # wait at most one hour per tile
                except multiprocessing.TimeoutError:
                    print "Timeout while computing tile " + str(r)

            print 'global finalization...'
            global_finalization(tiles_full_info)

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]


def main(config_file):
    """
    Launches s2p with the parameters given by a json file.

    Args:
        config_file: path to the config json file
    """
    t0 = time.time()

    # run the pipeline
    map_processing(config_file)

    # measure total runtime
    t = int(time.time() - t0)
    h = t / 3600
    m = (t / 60) % 60
    s = t % 60
    print "Total runtime: %dh:%dm:%ds" % (h, m, s)

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
