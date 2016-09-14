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
from python import pointing_accuracy
from python import rectification
from python import block_matching
from python import masking
from python import triangulation
from python import fusion
from python import visualisation


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


def pointing_correction_pair(tile, i=None):
    """
    Compute the translation that correct the pointing error on a pair of tiles.

    Args:
        tile: dictionary containing the information needed to process the tile
        i: index of the processed pair. If None, there's only one pair.
    """
    x, y, w, h = tile['coordinates']
    out_dir = os.path.join(tile['directory'], 'pair_{}'.format(i)) if i else tile['directory']
    img1 = cfg['images'][0]['img']
    rpc1 = cfg['images'][0]['rpc']
    img2 = cfg['images'][i]['img'] if i else cfg['images'][1]['img']
    rpc2 = cfg['images'][i]['rpc'] if i else cfg['images'][1]['rpc']

    if cfg['skip_existing'] and os.path.isfile(os.path.join(out_dir,
                                                            'pointing.txt')):
        print 'pointing correction done on tile {} {}'.format(x, y),
        print 'pair {}'.format(i) if i else ''
        return

    # correct pointing error
    print 'correcting pointing on tile {} {}'.format(x, y),
    print 'pair {}...'.format(i) if i else '...'
    A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, x, y, w, h)

    if A is not None:  # A is the correction matrix
        np.savetxt(os.path.join(out_dir, 'pointing.txt'), A, fmt='%6.3f')
    if m is not None:  # m is the list of sift matches
        np.savetxt(os.path.join(out_dir, 'sift_matches.txt'), m, fmt='%9.3f')
        np.savetxt(os.path.join(out_dir, 'center_keypts_sec.txt'),
                   np.mean(m[:, 2:], 0), fmt='%9.3f')
        if cfg['debug']:
            visualisation.plot_matches(img1, img2, rpc1, rpc2, m, x, y, w, h,
                                       os.path.join(out_dir,
                                                    'sift_matches_plot.png'))


def global_pointing_correction(tiles):
    """
    Compute the global pointing corrections for each pair of images.

    Args:
        tiles: list of tile dictionaries
    """
    if len(cfg['images']) == 2:
        out = os.path.join(cfg['out_dir'], 'global_pointing.txt')
        if not (os.path.isfile(out) and cfg['skip_existing']):
            np.savetxt(out, pointing_accuracy.global_from_local(t['directory']
                                                                for t in tiles),
                       fmt='%12.6f')
    else:
        for i in xrange(1, len(cfg['images'])):
            out = os.path.join(cfg['out_dir'], 'global_pointing_pair_%d.txt' % i)
            if not (os.path.isfile(out) and cfg['skip_existing']):
                l = [os.path.join(t['directory'], 'pair_%d' % i) for t in tiles]
                np.savetxt(out, pointing_accuracy.global_from_local(l),
                           fmt='%12.6f')


def rectification_pair(tile, i=None):
    """
    Rectify a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair. If None, there's only one pair.
    """
    out_dir = os.path.join(tile['directory'], 'pair_{}'.format(i)) if i else tile['directory']
    x, y, w, h = tile['coordinates']
    img1 = cfg['images'][0]['img']
    rpc1 = cfg['images'][0]['rpc']
    img2 = cfg['images'][i]['img'] if i else cfg['images'][1]['img']
    rpc2 = cfg['images'][i]['rpc'] if i else cfg['images'][1]['rpc']
    pointing = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_{}.txt'.format(i) if i else
                            'global_pointing.txt')

    outputs = ['disp_min_max.txt', 'rectified_ref.tif', 'rectified_sec.tif']
    if cfg['skip_existing'] and all(os.path.isfile(os.path.join(out_dir, f)) for
                                    f in outputs):
        print 'rectification done on tile {} {}'.format(x, y),
        print 'pair {}'.format(i) if i else ''
        return

    print 'rectifying tile {} {}'.format(x, y),
    print 'pair {}...'.format(i) if i else '...'
    try:
        A = np.loadtxt(os.path.join(out_dir, 'pointing.txt'))
    except IOError:
        A = pointing
    try:
        m = np.loadtxt(os.path.join(out_dir, 'sift_matches.txt'))
    except IOError:
        m = None
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m,
                                                            margin=cfg['horizontal_margin'])
    np.savetxt(os.path.join(out_dir, 'H_ref.txt'), H1, fmt='%12.6f')
    np.savetxt(os.path.join(out_dir, 'H_sec.txt'), H2, fmt='%12.6f')
    np.savetxt(os.path.join(out_dir, 'disp_min_max.txt'), [disp_min, disp_max],
                            fmt='%3.1f')


def disparity_pair(tile, i=None):
    """
    Compute the disparity of a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair. If None, there's only one pair.
    """
    out_dir = os.path.join(tile['directory'], 'pair_{}'.format(i)) if i else tile['directory']
    x, y = tile['coordinates'][:2]

    outputs = ['rectified_mask.png', 'rectified_disp.tif']
    if cfg['skip_existing'] and all(os.path.isfile(os.path.join(out_dir, f)) for
                                    f in outputs):
        print 'disparity estimation done on tile {} {}'.format(x, y),
        print 'pair {}'.format(i) if i else ''
        return

    print 'estimating disparity on tile {} {}'.format(x, y),
    print 'pair {}...'.format(i) if i else '...'
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    disp_min, disp_max = np.loadtxt(os.path.join(out_dir, 'disp_min_max.txt'))
    if cfg['disp_min'] is not None: disp_min = cfg['disp_min']
    if cfg['disp_max'] is not None: disp_max = cfg['disp_max']
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # add margin around masked pixels
    masking.erosion(mask, mask, cfg['msk_erosion'])


def triangulation_pair(tile, i=None):
    """
    Triangulate a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair. If None, there's only one pair.
    """
    out_dir = os.path.join(tile['directory'], 'pair_{}'.format(i)) if i else tile['directory']
    height_map = os.path.join(out_dir, 'height_map.tif')
    x, y, w, h = tile['coordinates']

    if cfg['skip_existing'] and os.path.isfile(height_map):
        print 'triangulation done on tile {} {}'.format(x, y),
        print 'pair {}'.format(i) if i else ''
        return

    print 'triangulating tile {} {}'.format(x, y),
    print 'pair {}...'.format(i) if i else '...'
    rpc1 = cfg['images'][0]['rpc']
    rpc2 = cfg['images'][i]['rpc'] if i else cfg['images'][1]['rpc']
    H_ref = os.path.join(out_dir, 'H_ref.txt')
    H_sec = os.path.join(out_dir, 'H_sec.txt')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    rpc_err = os.path.join(out_dir, 'rpc_err.tif')
    out_mask = os.path.join(tile['directory'], 'cloud_water_image_domain_mask.png')
    pointing = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_{}.txt'.format(i) if i else
                            'global_pointing.txt')
    triangulation.height_map(height_map, x, y, w, h, cfg['subsampling_factor'],
                             rpc1, rpc2, H_ref, H_sec, disp, mask, rpc_err,
                             out_mask, pointing)


def compute_mean_heights(tile):
    """
    """
    w, h = tile['coordinates'][2:]
    n = len(cfg['images']) - 1
    maps = np.empty((h, w, n))
    for i in xrange(n):
        f = gdal.Open(os.path.join(tile['directory'], 'pair_%d' % (i + 1), 'height_map.tif'))
        maps[:, :, i] = f.GetRasterBand(1).ReadAsArray()
        f = None  # this is the gdal way of closing files

    validity_mask = maps.sum(axis=2)  # sum to propagate nan values
    validity_mask += 1 - validity_mask  # 1 on valid pixels, and nan on invalid
    return [np.nanmean(validity_mask * maps[:, :, i]) for i in range(n)]


def tile_fusion(tile, mean_heights_global):
    """
    Merge the height maps computed for each image pair and generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
        mean_heights_global: list containing the means of all the global height
            maps
    """
    tile_dir = tile['directory']
    nb_pairs = len(mean_heights_global)
    height_maps = [os.path.join(tile_dir, 'pair_%d' % (i + 1), 'height_map.tif')
                   for i in xrange(nb_pairs)]

    # remove spurious matches
    if cfg['cargarse_basura']:
        for img in height_maps:
            common.cargarse_basura(img, img)

    # merge the height maps (applying mean offset to register)
    fusion.merge_n(os.path.join(tile_dir, 'height_map.tif'), height_maps,
                   mean_heights_global, averaging=cfg['fusion_operator'],
                   threshold=cfg['fusion_thresh'])


def compute_ply(tile):
    """
    Generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
    """
    tile_dir = tile['directory']
    x, y, w, h = tile['coordinates']
    z = cfg['subsampling_factor']

    # H is the homography transforming the coordinates system of the original
    # full size image into the coordinates system of the crop
    H = np.dot(np.diag([1 / z, 1 / z, 1]), common.matrix_translation(-x, -y))
    colors = os.path.join(tile_dir, 'ref.png')
    if cfg['images'][0]['clr']:
        common.image_crop_tif(cfg['images'][0]['clr'], x, y, w, h, colors)
    else:
        common.image_qauto(common.image_crop_tif(cfg['images'][0]['img'], x, y,
                                                 w, h), colors)
    triangulation.compute_point_cloud(os.path.join(tile_dir, 'cloud.ply'),
                                      os.path.join(tile_dir, 'height_map.tif'),
                                      cfg['images'][0]['rpc'], H, colors,
                                      utm_zone=cfg['utm_zone'],
                                      llbbx=tuple(cfg['ll_bbx']))


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


def tilewise_wrapper(fun, *args, **kwargs):
    """
    """
    if not cfg['debug']:  # redirect stdout and stderr to log file
        # the last argument '0' disables buffering
        f = open(kwargs['stdout'], 'a', 0)
        sys.stdout = f
        sys.stderr = f

    try:
        out = fun(*args)
    except Exception:
        print("Exception in %s" % fun.__name__)
        traceback.print_exc()
        raise

    common.garbage_cleanup()
    if not cfg['debug']:  # close logs
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        f.close()

    return out


def launch_parallel_calls(fun, list_of_args, nb_workers, *extra_args):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel.
        list_of_args: list of (first positional) arguments passed to fun, one
            per call
        nb_workers: number of calls run simultaneously
        extra_args (optional): tuple containing extra arguments to be passed to
            fun (same value for all calls)
    """
    results = []
    outputs = []
    show_progress.counter = 0
    show_progress.total = len(list_of_args)
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        if type(x) == tuple:  # we expect x = (tile_dictionary, pair_id)
            args = (fun,) + x + extra_args
            log = os.path.join(x[0]['directory'], 'pair_%d' % x[1], 'stdout.log')
        else:  # we expect x = tile_dictionary
            args = (fun, x) + extra_args
            log = os.path.join(x['directory'], 'stdout.log')
        results.append(pool.apply_async(tilewise_wrapper, args=args,
                                        kwds={'stdout': log},
                                        callback=show_progress))

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


def main(config_file):
    """
    Launch the entire s2p pipeline with the parameters given in a json file.

    Args:
        config_file: path to a json configuration file
    """
    print_elapsed_time.t0 = datetime.datetime.now()

    # initialization
    initialization.build_cfg(config_file)
    initialization.make_dirs()
    tiles = initialization.tiles_full_info()
    if len(cfg['images']) > 2:
        tiles_pairs = [(t, i) for i in xrange(1, len(cfg['images'])) for t in
                       tiles]
    else:
        tiles_pairs = tiles

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_nb_threads']:
        nb_workers = min(nb_workers, cfg['max_nb_threads'])

    # omp_num_threads: should not exceed nb_workers when multiplied by the
    # number of tiles
    # cfg['omp_num_threads'] = max(1, int(nb_workers / len(tiles)))

    # do the job
    print '\ncorrecting pointing locally...'
    launch_parallel_calls(pointing_correction_pair, tiles_pairs, nb_workers)
    print_elapsed_time()

    print '\ncorrecting pointing globally...'
    global_pointing_correction(tiles)
    print_elapsed_time()

    print '\nrectifying tiles...'
    launch_parallel_calls(rectification_pair, tiles_pairs, nb_workers)
    print_elapsed_time()

    print '\nrunning stereo matching...'
    launch_parallel_calls(disparity_pair, tiles_pairs, nb_workers)
    print_elapsed_time()

    print '\ntriangulating tiles...'
    launch_parallel_calls(triangulation_pair, tiles_pairs, nb_workers)
    print_elapsed_time()

    if len(cfg['images']) > 2:
        print '\nregistering height maps...'
        mean_heights_local = launch_parallel_calls(compute_mean_heights, tiles,
                                                   nb_workers)
        print_elapsed_time()

        print '\ncompute global pairwise height offsets...'
        mean_heights_global = np.mean(mean_heights_local, axis=0)
        print_elapsed_time()

        print '\nmerging height maps...'
        launch_parallel_calls(tile_fusion, tiles, nb_workers,
                              mean_heights_global)
        print_elapsed_time()

    print '\ncomputing point clouds...'
    launch_parallel_calls(compute_ply, tiles, nb_workers)
    print_elapsed_time()

    print '\ncompute dsm...'
    compute_dsm(tiles)
    print_elapsed_time()

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
      > %s config.json
        Launches the s2p pipeline.

      All the parameters, paths to input and output files, are defined in
      the json configuration file.

    """ % script_name
    sys.exit()


if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print_help_and_exit(sys.argv[0])
