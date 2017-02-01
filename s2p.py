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

from __future__ import print_function
import sys
import os.path
import datetime
import argparse
import numpy as np
import multiprocessing
from osgeo import gdal

gdal.UseExceptions()

from s2plib.config import cfg
from s2plib import common
from s2plib import parallel
from s2plib import initialization
from s2plib import pointing_accuracy
from s2plib import rectification
from s2plib import block_matching
from s2plib import masking
from s2plib import triangulation
from s2plib import fusion
from s2plib import visualisation


def pointing_correction(tile, i=None):
    """
    Compute the translation that corrects the pointing error on a pair of tiles.

    Args:
        tile: dictionary containing the information needed to process the tile
        i: index of the processed pair. If None, there's only one pair.
    """
    x, y, w, h = tile['coordinates']
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i) if i else '')
    img1 = cfg['images'][0]['img']
    rpc1 = cfg['images'][0]['rpc']
    img2 = cfg['images'][i]['img'] if i else cfg['images'][1]['img']
    rpc2 = cfg['images'][i]['rpc'] if i else cfg['images'][1]['rpc']

    if cfg['skip_existing'] and os.path.isfile(os.path.join(out_dir,
                                                            'pointing.txt')):
        print('pointing correction done on tile {} {}'.format(x, y), end=' ')
        print('pair {}'.format(i) if i else '')
        return

    # correct pointing error
    print('correcting pointing on tile {} {}'.format(x, y), end=' ')
    print('pair {}...'.format(i) if i else '...')
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
                                                    'sift_matches_pointing.png'))


def global_pointing_correction(tiles):
    """
    Compute the global pointing corrections for each pair of images.

    Args:
        tiles: list of tile dictionaries
    """
    if len(cfg['images']) == 2:
        out = os.path.join(cfg['out_dir'], 'global_pointing.txt')
        if not (os.path.isfile(out) and cfg['skip_existing']):
            np.savetxt(out, pointing_accuracy.global_from_local(t['dir']
                                                                for t in tiles),
                       fmt='%12.6f')
    else:
        for i in range(1, len(cfg['images'])):
            out = os.path.join(cfg['out_dir'], 'global_pointing_pair_%d.txt' % i)
            if not (os.path.isfile(out) and cfg['skip_existing']):
                l = [os.path.join(t['dir'], 'pair_%d' % i) for t in tiles]
                np.savetxt(out, pointing_accuracy.global_from_local(l),
                           fmt='%12.6f')


def rectification_pair(tile, i=None):
    """
    Rectify a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair. If None, there's only one pair.
    """
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i) if i else '')
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
        print('rectification done on tile {} {}'.format(x, y), end=' ')
        print('pair {}'.format(i) if i else '')
        return

    print('rectifying tile {} {}'.format(x, y), end=' ')
    print('pair {}...'.format(i) if i else '...')
    try:
        A = np.loadtxt(os.path.join(out_dir, 'pointing.txt'))
    except IOError:
        A = np.loadtxt(pointing)
    try:
        m = np.loadtxt(os.path.join(out_dir, 'sift_matches.txt'))
    except IOError:
        m = None
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m,
                                                            hmargin=cfg['horizontal_margin'],
                                                            vmargin=cfg['vertical_margin'])
    np.savetxt(os.path.join(out_dir, 'H_ref.txt'), H1, fmt='%12.6f')
    np.savetxt(os.path.join(out_dir, 'H_sec.txt'), H2, fmt='%12.6f')
    np.savetxt(os.path.join(out_dir, 'disp_min_max.txt'), [disp_min, disp_max],
                            fmt='%3.1f')


def stereo_matching(tile, i=None):
    """
    Compute the disparity of a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair. If None, there's only one pair.
    """
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i) if i else '')
    x, y = tile['coordinates'][:2]

    outputs = ['rectified_mask.png', 'rectified_disp.tif']
    if cfg['skip_existing'] and all(os.path.isfile(os.path.join(out_dir, f)) for
                                    f in outputs):
        print('disparity estimation done on tile {} {}'.format(x, y), end=' ')
        print('pair {}'.format(i) if i else '')
        return

    print('estimating disparity on tile {} {}'.format(x, y), end=' ')
    print('pair {}...'.format(i) if i else '...')
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


def disparity_to_height(tile, i):
    """
    Compute a height map from the disparity map of a pair of image tiles.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair.
    """
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i))
    height_map = os.path.join(out_dir, 'height_map.tif')
    x, y, w, h = tile['coordinates']

    if cfg['skip_existing'] and os.path.isfile(height_map):
        print('triangulation done on tile {} {} pair {}'.format(x, y, i))
        return

    print('triangulating tile {} {} pair {}...'.format(x, y, i))
    rpc1 = cfg['images'][0]['rpc']
    rpc2 = cfg['images'][i]['rpc']
    H_ref = os.path.join(out_dir, 'H_ref.txt')
    H_sec = os.path.join(out_dir, 'H_sec.txt')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    rpc_err = os.path.join(out_dir, 'rpc_err.tif')
    out_mask = os.path.join(tile['dir'], 'cloud_water_image_domain_mask.png')
    pointing = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_{}.txt'.format(i))
    triangulation.height_map(height_map, x, y, w, h, cfg['subsampling_factor'],
                             rpc1, rpc2, H_ref, H_sec, disp, mask, rpc_err,
                             out_mask, pointing)


def disparity_to_ply(tile):
    """
    Compute a point cloud from the disparity map of a pair of image tiles.

    Args:
        tile: dictionary containing the information needed to process a tile.
    """
    out_dir = tile['dir']
    ply_file = os.path.join(out_dir, 'cloud.ply')
    x, y, w, h = tile['coordinates']

    if cfg['skip_existing'] and os.path.isfile(ply_file):
        print('triangulation done on tile {} {}'.format(x, y))
        return

    print('triangulating tile {} {}...'.format(x, y))
    rpc1 = cfg['images'][0]['rpc']
    rpc2 = cfg['images'][1]['rpc']
    H_ref = os.path.join(out_dir, 'H_ref.txt')
    H_sec = os.path.join(out_dir, 'H_sec.txt')
    pointing = os.path.join(cfg['out_dir'], 'global_pointing.txt')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask_rect = os.path.join(out_dir, 'rectified_mask.png')
    mask_orig = os.path.join(out_dir, 'cloud_water_image_domain_mask.png')

    # prepare the image needed to colorize point cloud
    colors = os.path.join(out_dir, 'rectified_ref.png')
    if cfg['images'][0]['clr']:
        hom = np.loadtxt(H_ref)
        roi = [[x, y], [x+w, y], [x+w, y+h], [x, y+h]]
        ww, hh = common.bounding_box2D(common.points_apply_homography(hom, roi))[2:]
        tmp = common.tmpfile('.tif')
        common.image_apply_homography(tmp, cfg['images'][0]['clr'], hom,
                                      ww + 2*cfg['horizontal_margin'],
                                      hh + 2*cfg['vertical_margin'])
        common.image_qauto(tmp, colors)
    else:
        common.image_qauto(os.path.join(out_dir, 'rectified_ref.tif'), colors)

    # compute the point cloud
    triangulation.disp_map_to_point_cloud(ply_file, disp, mask_rect, rpc1, rpc2,
                                          H_ref, H_sec, pointing, colors,
                                          utm_zone=cfg['utm_zone'],
                                          llbbx=tuple(cfg['ll_bbx']),
                                          xybbx=(x, x+w, y, y+h),
                                          xymsk=mask_orig)


def mean_heights(tile):
    """
    """
    w, h = tile['coordinates'][2:]
    z = cfg['subsampling_factor']
    n = len(cfg['images']) - 1
    maps = np.empty((int(h/z), int(w/z), n))
    for i in range(n):
        f = gdal.Open(os.path.join(tile['dir'], 'pair_%d' % (i + 1),
                                   'height_map.tif'))
        maps[:, :, i] = f.GetRasterBand(1).ReadAsArray()
        f = None  # this is the gdal way of closing files

    validity_mask = maps.sum(axis=2)  # sum to propagate nan values
    validity_mask += 1 - validity_mask  # 1 on valid pixels, and nan on invalid
    return [np.nanmean(validity_mask * maps[:, :, i]) for i in range(n)]


def heights_fusion(tile, mean_heights_global):
    """
    Merge the height maps computed for each image pair and generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
        mean_heights_global: list containing the means of all the global height
            maps
    """
    tile_dir = tile['dir']
    nb_pairs = len(mean_heights_global)
    height_maps = [os.path.join(tile_dir, 'pair_%d' % (i + 1), 'height_map.tif')
                   for i in range(nb_pairs)]

    # remove spurious matches
    if cfg['cargarse_basura']:
        for img in height_maps:
            common.cargarse_basura(img, img)

    # merge the height maps (applying mean offset to register)
    fusion.merge_n(os.path.join(tile_dir, 'height_map.tif'), height_maps,
                   mean_heights_global, averaging=cfg['fusion_operator'],
                   threshold=cfg['fusion_thresh'])


def heights_to_ply(tile):
    """
    Generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
    """
    out_dir = tile['dir']
    x, y, w, h = tile['coordinates']
    z = cfg['subsampling_factor']
    plyfile = os.path.join(out_dir, 'cloud.ply')
    if cfg['skip_existing'] and os.path.isfile(plyfile):
        print('ply file already exists for tile {} {}'.format(x, y))
        return

    # H is the homography transforming the coordinates system of the original
    # full size image into the coordinates system of the crop
    H = np.dot(np.diag([1 / z, 1 / z, 1]), common.matrix_translation(-x, -y))
    colors = os.path.join(out_dir, 'ref.png')
    if cfg['images'][0]['clr']:
        common.image_crop_gdal(cfg['images'][0]['clr'], x, y, w, h, colors)
    else:
        common.image_qauto(common.image_crop_gdal(cfg['images'][0]['img'], x, y,
                                                 w, h), colors)
    common.image_safe_zoom_fft(colors, z, colors)
    triangulation.height_map_to_point_cloud(plyfile, os.path.join(out_dir,
                                                                  'height_map.tif'),
                                            cfg['images'][0]['rpc'], H, colors,
                                            utm_zone=cfg['utm_zone'],
                                            llbbx=tuple(cfg['ll_bbx']))


def plys_to_dsm(tiles):
    """
    """
    out_dsm = os.path.join(cfg['out_dir'], 'dsm.tif')
    clouds = ' '.join(os.path.join(t['dir'], 'cloud.ply') for t in tiles)
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
    Produce a single multiscale point cloud for the whole processed region.

    Args:
        tiles: list of tiles dictionaries
    """
    if common.which('LidarPreprocessor') is None:
        return
    plys = [os.path.join(os.path.abspath(t['dir']), 'cloud.ply') for t in tiles]
    common.lidar_preprocessor(os.path.join(cfg['out_dir'],
                                           'cloud.lidar_viewer'), plys)


ALL_STEPS = ['local-pointing', 'global-pointing', 'rectification', 'matching',
             'triangulation', 'dsm-rasterization', 'lidar-preprocessor']


def main(config_file, steps=ALL_STEPS):
    """
    Launch the s2p pipeline with the parameters given in a json file.

    Args:
        config_file: path to a json configuration file
        steps: either a string (single step) or a list of strings (several
            steps)
    """
    common.print_elapsed_time.t0 = datetime.datetime.now()
    initialization.build_cfg(config_file)
    initialization.make_dirs()

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes']:
        nb_workers = min(nb_workers, cfg['max_processes'])
    cfg['max_processes'] = nb_workers

    print('\ndiscarding masked tiles...')
    tw, th = initialization.adjust_tile_size()
    tiles = initialization.tiles_full_info(tw, th)
    n = len(cfg['images'])
    if n > 2:
        tiles_pairs = [(t, i) for i in range(1, n) for t in tiles]
    else:
        tiles_pairs = tiles

    # omp_num_threads should not exceed nb_workers when multiplied by len(tiles)
    cfg['omp_num_threads'] = max(1, int(nb_workers / len(tiles_pairs)))

    if 'local-pointing' in steps:
        print('correcting pointing locally...')
        parallel.launch_calls(pointing_correction, tiles_pairs, nb_workers)

    if 'global-pointing' in steps:
        print('correcting pointing globally...')
        global_pointing_correction(tiles)
        common.print_elapsed_time()

    if 'rectification' in steps:
        print('rectifying tiles...')
        parallel.launch_calls(rectification_pair, tiles_pairs, nb_workers)

    if 'matching' in steps:
        print('running stereo matching...')
        parallel.launch_calls(stereo_matching, tiles_pairs, nb_workers)

    if 'triangulation' in steps:
        if n > 2:
            print('computing height maps...')
            parallel.launch_calls(disparity_to_height, tiles_pairs, nb_workers)

            print('registering height maps...')
            mean_heights_local = parallel.launch_calls(mean_heights, tiles,
                                                       nb_workers)

            print('computing global pairwise height offsets...')
            mean_heights_global = np.nanmean(mean_heights_local, axis=0)

            print('merging height maps...')
            parallel.launch_calls(heights_fusion, tiles, nb_workers,
                                  mean_heights_global)

            print('computing point clouds...')
            parallel.launch_calls(heights_to_ply, tiles, nb_workers)

        else:
            print('triangulating tiles...')
            parallel.launch_calls(disparity_to_ply, tiles, nb_workers)

    if 'dsm-rasterization' in steps:
        print('computing DSM...')
        plys_to_dsm(tiles)
        common.print_elapsed_time()

    if 'lidar-preprocessor' in steps:
        print('lidar preprocessor...')
        lidar_preprocessor(tiles)
        common.print_elapsed_time()

    # cleanup
    common.garbage_cleanup()
    common.print_elapsed_time(since_first_call=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: Satellite Stereo '
                                                  'Pipeline'))
    parser.add_argument('config', metavar='config.json',
                        help=('path to a json file containing the paths to '
                              'input and output files and the algorithm '
                              'parameters'))
    parser.add_argument('--step', type=str, choices=ALL_STEPS,
                        default=ALL_STEPS)
    args = parser.parse_args()
    main(args.config, args.step)
