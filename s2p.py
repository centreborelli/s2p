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
import glob
import shutil
import os.path
import datetime
import traceback
import numpy as np
import multiprocessing
from osgeo import gdal

gdal.UseExceptions()

from python.config import cfg
from python import common
from python import initialization
from python import preprocess
from python import globalvalues
from python import process
from python import fusion
from python import triangulation


def show_progress(a):
    """
    Print the number of tiles that have been processed.

    Args:
        a: useless argument, but since this function is used as a callback by
            apply_async, it has to take one argument.
    """
    show_progress.counter += 1
    status = "done {:{fill}{width}} / {} tiles".format(show_progress.counter,
                                                       show_progress.total,
                                                       fill='',
                                                       width=len(str(show_progress.total)))
    if show_progress.counter < show_progress.total:
        status += chr(8) * len(status)
    else:
        status += '\n'
    sys.stdout.write(status)
    sys.stdout.flush()


def print_elapsed_time(since_first_call=False):
    """
    Print the elapsed time since the last call or since the first call.

    Args:
        since_first_call:
    """
    t2 = datetime.datetime.now()
    if since_first_call:
        print "Total elapsed time:", t2 - print_elapsed_time.t0
    else:
        try:
            print "Elapsed time:", t2 - print_elapsed_time.t1
        except AttributeError:
            print t2 - print_elapsed_time.t0
    print_elapsed_time.t1 = t2


def preprocess_tile(tile):
    """
    Compute pointing corrections and extrema intensities for a single tile.

    Args:
        tile: dictionary containing all the information needed to process a
            tile.
    """
    # redirect stdout and stderr to log file
    if not cfg['debug']:
        f = open(os.path.join(tile['directory'], 'stdout.log'), 'w', 0)  # 0 for no buffering
        sys.stdout = f
        sys.stderr = f

    try:
        preprocess.pointing_correction(tile)
    except Exception:
        print("Exception in preprocessing tile:")
        traceback.print_exc()
        raise

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        f.close()


def global_values(tiles):
    """
    Compute the global pointing correction and extrema intensities for the ROI.
    """
    globalvalues.pointing_correction(tiles)


def process_tile_pair(tile, pair_id):
    """
    Process a pair of images on a given tile.

    It includes rectification, disparity estimation and triangulation.

    Args:
        tile: dictionary containing all the information needed to process a
            tile.
        pair_id: index of the pair to process
    """
    # read all the information
    tile_dir = tile['directory']
    col, row, tw, th = tile['coordinates']
    images = cfg['images']
    img1, rpc1 = images[0]['img'], images[0]['rpc']
    img2, rpc2 = images[pair_id]['img'], images[pair_id]['rpc']
    out_dir = os.path.join(tile_dir, 'pair_%d' % pair_id)
    A_global = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_%d.txt' % pair_id)

    print 'processing tile %d %d...' % (col, row)

    # rectification
    if (cfg['skip_existing'] and
            os.path.isfile(os.path.join(out_dir, 'disp_min_max.txt')) and
            os.path.isfile(os.path.join(out_dir, 'rectified_ref.tif')) and
            os.path.isfile(os.path.join(out_dir, 'rectified_sec.tif'))):
        print '\trectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\trectifying tile %d %d (pair %d)...' % (col, row, pair_id)
        process.rectify(out_dir, np.loadtxt(A_global), img1, rpc1,
                        img2, rpc2, col, row, tw, th, cfg['horizontal_margin'])

    # disparity estimation
    if (cfg['skip_existing'] and
            os.path.isfile(os.path.join(out_dir, 'rectified_mask.png')) and
            os.path.isfile(os.path.join(out_dir, 'rectified_disp.tif'))):
        print '\tdisparity estimation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\testimating disparity on tile %d %d (pair %d)...' % (col, row, pair_id)
        process.disparity(out_dir)

    # triangulation
    height_map = os.path.join(out_dir, 'height_map.tif')
    if cfg['skip_existing'] and os.path.isfile(height_map):
        print '\tfile %s already there, skip triangulation' % height_map
    else:
        print '\ttriangulating tile %d %d (pair %d)...' % (col, row, pair_id)
        H_ref = os.path.join(out_dir, 'H_ref.txt')
        H_sec = os.path.join(out_dir, 'H_sec.txt')
        disp = os.path.join(out_dir, 'rectified_disp.tif')
        mask = os.path.join(out_dir, 'rectified_mask.png')
        out_mask = os.path.join(os.path.dirname(out_dir),
                                'cloud_water_image_domain_mask.png')
        rpc_err = os.path.join(out_dir, 'rpc_err.tif')
        triangulation.height_map(height_map, col, row, tw, th,
                                 cfg['subsampling_factor'], rpc1, rpc2, H_ref,
                                 H_sec, disp, mask, rpc_err, out_mask, A_global)


def process_tile(tile):
    """
    Compute a height map on the tile for each image pair.

    Args:
        tile: a dictionary that provides all you need to process a tile
    """
    tile_dir = tile['directory']
    w, h = tile['coordinates'][2:]
    n = len(cfg['images']) - 1

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        l = open(os.path.join(tile_dir, 'stdout.log'), 'a', 0)  # '0' for no buffering
        sys.stdout = l
        sys.stderr = l

    try:  # process each pair to get a height map
        for i in xrange(n):
            process_tile_pair(tile, i + 1)
    except Exception:
        print("Exception in processing tile:")
        traceback.print_exc()
        raise

    maps = np.empty((h, w, n))
    for i in xrange(n):
        f = gdal.Open(os.path.join(tile_dir, 'pair_%d' % (i + 1), 'height_map.tif'))
        maps[:, :, i] = f.GetRasterBand(1).ReadAsArray()
        f = None  # this is the gdal way of closing files

    validity_mask = maps.sum(axis=2)  # sum to propagate nan values
    validity_mask += 1 - validity_mask  # 1 on valid pixels, and nan on invalid
    mean_heights = [np.nanmean(validity_mask * maps[:, :, i]) for i in range(n)]

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        l.close()

    return mean_heights


def tile_fusion_and_ply(tile, mean_heights_global):
    """
    Merge the height maps computed for each image pair and generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
        mean_heights_global: list containing the means of all the global height
            maps
    """
    tile_dir = tile['directory']
    nb_pairs = len(mean_heights_global)

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        f = open(os.path.join(tile_dir, 'stdout.log'), 'a', 0)  # '0' for no buffering
        sys.stdout = f
        sys.stderr = f

    height_maps = [os.path.join(tile_dir, 'pair_%d' % (i + 1), 'height_map.tif')
                   for i in xrange(nb_pairs)]
    try:
        # remove spurious matches
        if cfg['cargarse_basura']:
            for img in height_maps:
                process.cargarse_basura(img, img)

        # merge the height maps (applying mean offset to register)
        fusion.merge_n(os.path.join(tile_dir, 'height_map.tif'), height_maps,
                       mean_heights_global, averaging=cfg['fusion_operator'],
                       threshold=cfg['fusion_thresh'])

        # compute ply: H is the homography transforming the coordinates system of
        # the original full size image into the coordinates system of the crop
        x, y, w, h = tile['coordinates']
        z = cfg['subsampling_factor']
        H = np.dot(np.diag([1 / z, 1 / z, 1]), common.matrix_translation(-x, -y))
        colors = os.path.join(tile_dir, 'crop_ref.tif')
        common.image_crop_tif(cfg['images'][0]['img'], x, y, w, h, colors)
        common.image_qauto(colors, colors)
        triangulation.compute_point_cloud(os.path.join(tile_dir, 'cloud.ply'),
                                          os.path.join(tile_dir, 'height_map.tif'),
                                          cfg['images'][0]['rpc'], H,
                                          colors,
                                          utm_zone=cfg['utm_zone'],
                                          llbbx=tuple(cfg['ll_bbx']))
    except Exception:
        print("Exception in tile_fusion_and_ply:")
        traceback.print_exc()
        raise

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        f.close()


def compute_dsm(tiles):
    """
    """
    out_dsm = os.path.join(cfg['out_dir'], 'dsm.tif')
    clouds = ' '.join(os.path.join(t['directory'], 'cloud.ply') for t in
                      tiles)
    if 'utm_bbx' in cfg:
        bbx = cfg['utm_bbx']
        common.run("ls %s | plyflatten -bb \"%f %f %f %f \" %f %s" % (clouds,
                                                                      bbx[0],
                                                                      bbx[1],
                                                                      bbx[2],
                                                                      bbx[3],
                                                                      cfg['dsm_resolution'],
                                                                      out_dsm))
    else:
        common.run("ls %s | plyflatten %f %s" % (clouds, cfg['dsm_resolution'],
                                                 out_dsm))
        # ls files | ./bin/plyflatten [-c column] [-bb "xmin xmax ymin ymax"] resolution out.tif


def lidar_preprocessor(tiles):
    """
    Produce a single height map, DSM and point cloud for the whole ROI.

    The height maps associated to each pair, as well as the height map obtained
    by merging all the pairs, are stored as VRT files. The final DSM is obtained
    by projecting the 3D points from the ply files obtained on each tile. The
    final point cloud is obtained as the union of all the locally merged point
    clouds, in the LidarViewer format.

    Args:
        tiles: dictionary providing all the information about the
            processed tiles
    """
    # whole point cloud (LidarViewer format)
    if common.which('LidarPreprocessor'):
        out = os.path.join(cfg['out_dir'], 'cloud.lidar_viewer')
        plys = []
        for tile in tiles:
            plys.append(os.path.join(os.path.abspath(tile['directory']),
                                     'cloud.ply'))
        common.lidar_preprocessor(out, plys)

    # copy RPC xml files in the output directory
    for img in cfg['images']:
        shutil.copy2(img['rpc'], cfg['out_dir'])


def launch_parallel_calls(fun, list_of_args, nb_workers, extra_args=None):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel.
        list_of_args: list of (first positional) arguments passed to fun, one
            per call
        nb_workers: number of calls run simultaneously
        extra_args (optional, default is None): tuple containing extra arguments
            to be passed to fun (same value for all calls)
    """
    results = []
    outputs = []
    show_progress.counter = 0
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        args = (x,) + extra_args if extra_args else (x,)
        results.append(pool.apply_async(fun, args=args, callback=show_progress))

    for r in results:
        try:
            outputs.append(r.get(600))  # wait at most 10 min per call
        except multiprocessing.TimeoutError:
            print "Timeout while running %s" % str(r)
            outputs.append(None)
        except common.RunFailure as e:
            print "FAILED call: ", e.args[0]["command"]
            print "\toutput: ", e.args[0]["output"]
            outputs.append(None)
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(1)

    pool.close()
    pool.join()
    return outputs


def main(config_file, steps=range(1, 9)):
    """
    Launch the entire s2p pipeline with the parameters given in a json file.

    It is a succession of six steps:
        initialization
        preprocessing
        global_values
        processing
        compute dsms
        global_finalization

    Args:
        config_file: path to a json configuration file
        steps: list of integers between 1 and 8 specifying which steps to run.
            By default all steps are run.
    """
    print_elapsed_time.t0 = datetime.datetime.now()

    # initialization (has to be done whatever the queried steps)
    initialization.build_cfg(config_file)
    initialization.make_dirs()
    tiles = initialization.tiles_full_info()

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_nb_threads']:
        nb_workers = min(nb_workers, cfg['max_nb_threads'])

    # omp_num_threads: should not exceed nb_workers when multiplied by the
    # number of tiles
    # cfg['omp_num_threads'] = max(1, int(nb_workers / len(tiles)))

    # do the job
    if 2 in steps:
        print '\npreprocessing tiles...'
        show_progress.total = len(tiles)
        launch_parallel_calls(preprocess_tile, tiles, nb_workers)
        print_elapsed_time()

    if 3 in steps:
        print '\ncomputing global values...'
        global_values(tiles)
        print_elapsed_time()

    if 4 in steps:
        print '\nprocessing tiles...'
        show_progress.total = len(tiles)
        mean_heights_local = launch_parallel_calls(process_tile,
                                                   tiles, nb_workers)
        print_elapsed_time()

    if 5 in steps:
        print '\ncompute global pairwise height offsets...'
        mean_heights_global = np.mean(mean_heights_local, axis=0)
        print_elapsed_time()

    if 6 in steps:
        print '\nmerge height maps and compute ply clouds...'
        launch_parallel_calls(tile_fusion_and_ply, tiles, nb_workers,
                              (mean_heights_global,))
        print_elapsed_time()

    if 7 in steps:
        print '\ncompute dsm...'
        compute_dsm(tiles)
        print_elapsed_time()

    if 8 in steps:
        print '\nlidar preprocessor...'
        lidar_preprocessor(tiles)
        print_elapsed_time()

    # cleanup
    print_elapsed_time(since_first_call=True)
    common.garbage_cleanup()


def print_help_and_exit(script_name):
    """
    """
    print """
    Incorrect syntax, use:
      > %s config.json [steps (list of integer between 1 and 8)]
        1: initialization
        2: preprocessing (tilewise sift, local pointing correction)
        3: global-pointing
        4: processing (tilewise rectification, matching and triangulation)
        5: global height maps registration
        6: heights map merging and ply generation
        7: compute dsm from ply files
        8: lidarviewer
        Launches the s2p pipeline.

      All the parameters, paths to input and output files, are defined in
      the json configuration file.

    """ % script_name


if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3 and int(sys.argv[2]) in range(1, 9):
        main(sys.argv[1], int(sys.argv[2]))
    else:
        print_help_and_exit(sys.argv[0])
