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
import json
import datetime
import argparse
import numpy as np
import subprocess
import multiprocessing
from osgeo import gdal
import collections

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


def pointing_correction(tile, i):
    """
    Compute the translation that corrects the pointing error on a pair of tiles.

    Args:
        tile: dictionary containing the information needed to process the tile
        i: index of the processed pair
    """
    x, y, w, h = tile['coordinates']
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i))
    img1 = cfg['images'][0]['img']
    rpc1 = cfg['images'][0]['rpc']
    img2 = cfg['images'][i]['img']
    rpc2 = cfg['images'][i]['rpc']

    if cfg['skip_existing'] and os.path.isfile(os.path.join(out_dir,
                                                            'pointing.txt')):
        print('pointing correction done on tile {} {} pair {}'.format(x, y, i))
        return

    # correct pointing error
    print('correcting pointing on tile {} {} pair {}...'.format(x, y, i))
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
    for i in range(1, len(cfg['images'])):
        out = os.path.join(cfg['out_dir'], 'global_pointing_pair_%d.txt' % i)
        if not (os.path.isfile(out) and cfg['skip_existing']):
            l = [os.path.join(t['dir'], 'pair_%d' % i) for t in tiles]
            np.savetxt(out, pointing_accuracy.global_from_local(l),
                       fmt='%12.6f')
            if cfg['clean_intermediate']:
                for d in l:
                    common.remove(os.path.join(d, 'center_keypts_sec.txt'))


def rectification_pair(tile, i):
    """
    Rectify a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair
    """
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i))
    x, y, w, h = tile['coordinates']
    img1 = cfg['images'][0]['img']
    rpc1 = cfg['images'][0]['rpc']
    img2 = cfg['images'][i]['img']
    rpc2 = cfg['images'][i]['rpc']
    pointing = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_{}.txt'.format(i))

    outputs = ['disp_min_max.txt', 'rectified_ref.tif', 'rectified_sec.tif']
    if cfg['skip_existing'] and all(os.path.isfile(os.path.join(out_dir, f)) for
                                    f in outputs):
        print('rectification done on tile {} {} pair {}'.format(x, y, i))
        return

    print('rectifying tile {} {} pair {}...'.format(x, y, i))
    try:
        A = np.loadtxt(os.path.join(out_dir, 'pointing.txt'))
    except IOError:
        A = np.loadtxt(pointing)
    try:
        m = np.loadtxt(os.path.join(out_dir, 'sift_matches.txt'))
    except IOError:
        m = None

    x, y, w, h = tile['coordinates']

    for n in tile['neighborhood_dirs']:
        if n != tile['dir']:
            nei_dir = os.path.join(n, 'pair_{}'.format(i))
            sift_from_neighborhood = os.path.join(nei_dir, 'sift_matches.txt')
            try:
                m_n = np.loadtxt(sift_from_neighborhood)
                # added sifts in the ellipse of semi axes : (3*w/4, 3*h/4)
                m_n = m_n[np.where(np.linalg.norm([(m_n[:,0]-(x+w/2))/w,
                                                   (m_n[:,1]-(y+h/2))/h],
                                                  axis=0) < 3.0/4)]

                if m is None:
                    m = m_n
                else:
                    m = np.concatenate((m, m_n))
            except IOError:
                print('%s does not exist' % sift_from_neighborhood)
                pass

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

    if cfg['clean_intermediate']:
        common.remove(os.path.join(out_dir,'pointing.txt'))
        common.remove(os.path.join(out_dir,'sift_matches.txt'))


def stereo_matching(tile,i):
    """
    Compute the disparity of a pair of images on a given tile.

    Args:
        tile: dictionary containing the information needed to process a tile.
        i: index of the processed pair
    """
    out_dir = os.path.join(tile['dir'], 'pair_{}'.format(i))
    x, y = tile['coordinates'][:2]

    outputs = ['rectified_mask.png', 'rectified_disp.tif']
    if cfg['skip_existing'] and all(os.path.isfile(os.path.join(out_dir, f)) for
                                    f in outputs):
        print('disparity estimation done on tile {} {} pair {}'.format(x, y, i))
        return

    print('estimating disparity on tile {} {} pair {}...'.format(x, y, i))
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    disp_min, disp_max = np.loadtxt(os.path.join(out_dir, 'disp_min_max.txt'))

    # verifying non-epipolar_rectification is possible
    if cfg['matching_algorithm'] not in ['tvl1']:
        cfg['epipolar_rectification'] = True
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # add margin around masked pixels
    masking.erosion(mask, mask, cfg['msk_erosion'])

    if cfg['clean_intermediate']:
        if len(cfg['images']) > 2:
            common.remove(rect1)
        common.remove(rect2)
        common.remove(os.path.join(out_dir,'disp_min_max.txt'))


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

    if cfg['clean_intermediate']:
        common.remove(H_ref)
        common.remove(H_sec)
        common.remove(disp)
        common.remove(mask)
        common.remove(rpc_err)


def disparity_to_ply(tile):
    """
    Compute a point cloud from the disparity map of a pair of image tiles.

    Args:
        tile: dictionary containing the information needed to process a tile.
    """
    out_dir = os.path.join(tile['dir'])
    ply_file = os.path.join(out_dir, 'cloud.ply')
    plyextrema = os.path.join(out_dir, 'plyextrema.txt')
    x, y, w, h = tile['coordinates']
    rpc1 = cfg['images'][0]['rpc']
    rpc2 = cfg['images'][1]['rpc']

    if cfg['skip_existing'] and os.path.isfile(ply_file):
        print('triangulation done on tile {} {}'.format(x, y))
        return

    print('triangulating tile {} {}...'.format(x, y))
    # This function is only called when there is a single pair (pair_1)
    H_ref = os.path.join(out_dir, 'pair_1', 'H_ref.txt')
    H_sec = os.path.join(out_dir, 'pair_1', 'H_sec.txt')
    pointing = os.path.join(cfg['out_dir'], 'global_pointing_pair_1.txt')
    disp = os.path.join(out_dir, 'pair_1', 'rectified_disp.tif')
    mask_rect = os.path.join(out_dir, 'pair_1', 'rectified_mask.png')
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
        common.image_qauto(os.path.join(out_dir, 'pair_1', 'rectified_ref.tif'), colors)

    # compute the point cloud
    triangulation.disp_map_to_point_cloud(ply_file, disp, mask_rect, rpc1, rpc2,
                                          H_ref, H_sec, pointing, colors,
                                          utm_zone=cfg['utm_zone'],
                                          llbbx=tuple(cfg['ll_bbx']),
                                          xybbx=(x, x+w, y, y+h),
                                          xymsk=mask_orig)

    # compute the point cloud extrema (xmin, xmax, xmin, ymax)
    common.run("plyextrema %s %s" % (ply_file, plyextrema))

    if cfg['clean_intermediate']:
        common.remove(H_ref)
        common.remove(H_sec)
        common.remove(disp)
        common.remove(mask_rect)
        common.remove(mask_orig)
        common.remove(colors)
        common.remove(os.path.join(out_dir, 'pair_1', 'rectified_ref.tif'))


def mean_heights(tile):
    """
    """
    w, h = tile['coordinates'][2:]
    z = cfg['subsampling_factor']
    n = len(cfg['images']) - 1
    maps = np.empty((int(h/z), int(w/z), n))
    for i in range(n):
        try:
            f = gdal.Open(os.path.join(tile['dir'], 'pair_{}'.format(i + 1),
                                       'height_map.tif'))
            maps[:, :, i] = f.GetRasterBand(1).ReadAsArray()
            f = None  # this is the gdal way of closing files
        except RuntimeError:  # the file is not there
            maps[:, :, i] *= np.nan

    validity_mask = maps.sum(axis=2)  # sum to propagate nan values
    validity_mask += 1 - validity_mask  # 1 on valid pixels, and nan on invalid

    # save the n mean height values to a txt file in the tile directory
    np.savetxt(os.path.join(tile['dir'], 'local_mean_heights.txt'),
               [np.nanmean(validity_mask * maps[:, :, i]) for i in range(n)])


def global_mean_heights(tiles):
    """
    """
    local_mean_heights = [np.loadtxt(os.path.join(t['dir'], 'local_mean_heights.txt'))
                          for t in tiles]
    global_mean_heights = np.nanmean(local_mean_heights, axis=0)
    for i in range(len(cfg['images']) - 1):
        np.savetxt(os.path.join(cfg['out_dir'],
                                'global_mean_height_pair_{}.txt'.format(i+1)),
                   [global_mean_heights[i]])


def heights_fusion(tile):
    """
    Merge the height maps computed for each image pair and generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
    """
    tile_dir = tile['dir']
    height_maps = [os.path.join(tile_dir, 'pair_%d' % (i + 1), 'height_map.tif')
                   for i in range(len(cfg['images']) - 1)]

    # remove spurious matches
    if cfg['cargarse_basura']:
        for img in height_maps:
            common.cargarse_basura(img, img)

    # load global mean heights
    global_mean_heights = []
    for i in range(len(cfg['images']) - 1):
        x = np.loadtxt(os.path.join(cfg['out_dir'],
                                    'global_mean_height_pair_{}.txt'.format(i+1)))
        global_mean_heights.append(x)

    # merge the height maps (applying mean offset to register)
    fusion.merge_n(os.path.join(tile_dir, 'height_map.tif'), height_maps,
                   global_mean_heights, averaging=cfg['fusion_operator'],
                   threshold=cfg['fusion_thresh'])

    if cfg['clean_intermediate']:
        for f in height_maps:
            common.remove(f)


def heights_to_ply(tile):
    """
    Generate a ply cloud.

    Args:
        tile: a dictionary that provides all you need to process a tile
    """
    # merge the n-1 height maps of the tile (n = nb of images)
    heights_fusion(tile)

    # compute a ply from the merged height map
    out_dir = tile['dir']
    x, y, w, h = tile['coordinates']
    z = cfg['subsampling_factor']
    plyfile = os.path.join(out_dir, 'cloud.ply')
    plyextrema = os.path.join(out_dir, 'plyextrema.txt')
    height_map = os.path.join(out_dir, 'height_map.tif')
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
    triangulation.height_map_to_point_cloud(plyfile, height_map,
                                            cfg['images'][0]['rpc'], H, colors,
                                            utm_zone=cfg['utm_zone'],
                                            llbbx=tuple(cfg['ll_bbx']))

    # compute the point cloud extrema (xmin, xmax, xmin, ymax)
    common.run("plyextrema %s %s" % (plyfile, plyextrema))

    if cfg['clean_intermediate']:
        common.remove(height_map)
        common.remove(colors)
        common.remove(os.path.join(out_dir,
                                   'cloud_water_image_domain_mask.png'))

def global_srcwin(tiles):
    """
    """
    res = cfg['dsm_resolution']
    if 'utm_bbx' in cfg:
        bbx = cfg['utm_bbx']
        global_xoff = bbx[0]
        global_yoff = bbx[3]
        global_xsize = int(np.ceil((bbx[1]-bbx[0]) / res))
        global_ysize = int(np.ceil((bbx[3]-bbx[2]) / res))
    else:
        extrema = list()
        for t in tiles:
            plyextrema_file = os.path.join(t['dir'], "plyextrema.txt")
            if os.path.exists(plyextrema_file):
                extrema.append(np.loadtxt(plyextrema_file))
            else:
                extrema.append([np.nan]*4)

        xmin = np.nanmin(map(lambda x:x[0], extrema))
        xmax = np.nanmax(map(lambda x:x[1], extrema))
        ymin = np.nanmin(map(lambda x:x[2], extrema))
        ymax = np.nanmax(map(lambda x:x[3], extrema))

        global_xsize = int(1 + np.floor((xmax - xmin) / res))
        global_ysize = int(1 + np.floor((ymax - ymin) / res))
        global_xoff = (xmax + xmin - res * global_xsize) / 2
        global_yoff = (ymax + ymin + res * global_ysize) / 2

    np.savetxt(os.path.join(cfg['out_dir'], "global_srcwin.txt"),
               [global_xoff, global_yoff, global_xsize, global_ysize])

def plys_to_dsm(tile):
    """
    """
    out_dsm = os.path.join(tile['dir'], 'dsm.tif')
    global_srcwin = np.loadtxt(os.path.join(cfg['out_dir'],
                                            "global_srcwin.txt"))

    res = cfg['dsm_resolution']
    global_xoff, global_yoff, global_xsize, global_ysize = global_srcwin

    xmin, xmax, ymin, ymax = np.loadtxt(os.path.join(tile['dir'], "plyextrema.txt"))

    # compute xoff, yoff, xsize, ysize considering final dsm
    local_xoff = max(global_xoff,
                     global_xoff + np.floor((xmin - global_xoff) / res) * res)
    local_xsize = int(1 + np.floor((min(global_xoff + global_xsize * res, xmax) - local_xoff) / res))

    local_yoff = min(global_yoff,
                     global_yoff + np.ceil((ymax - global_yoff) / res) * res)
    local_ysize = int(1 - np.floor((max(global_yoff - global_ysize * res, ymin) - local_yoff) / res))

    clouds = '\n'.join(os.path.join(n_dir, 'cloud.ply') for n_dir in tile['neighborhood_dirs'])

    cmd = ['plyflatten', str(cfg['dsm_resolution']), out_dsm]
    cmd += ['-srcwin', '{} {} {} {}'.format(local_xoff, local_yoff,
                                            local_xsize, local_ysize)]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         stdin=subprocess.PIPE)
    q = p.communicate(input=clouds.encode())

    run_cmd = "ls %s | %s" % (clouds.replace('\n', ' '), " ".join(cmd))
    print ("\nRUN: %s" % run_cmd)

    if p.returncode != 0:
        raise common.RunFailure({"command": run_cmd, "environment": os.environ,
                                 "output": q})

    # ls files | ./bin/plyflatten [-c column] [-srcwin "xoff yoff xsize ysize"] resolution out.tif

def global_dsm(tiles):
    """
    """
    dsm_list = [os.path.join(t['dir'], 'dsm.tif') for t in tiles]
    out_dsm_vrt = os.path.join(cfg['out_dir'], 'dsm.vrt')
    out_dsm_tif = os.path.join(cfg['out_dir'], 'dsm.tif')

    common.run("gdalbuildvrt -vrtnodata nan  %s %s" % (out_dsm_vrt,
                                                       " ".join(dsm_list)))

    global_srcwin = np.loadtxt(os.path.join(cfg['out_dir'],
                                            "global_srcwin.txt"))
    res = cfg['dsm_resolution']
    xoff, yoff, xsize, ysize = global_srcwin

    common.run("gdal_translate -projwin %s %s %s %s %s %s" % (xoff,
                                                              yoff,
                                                              xoff + xsize * res,
                                                              yoff - ysize * res,
                                                              out_dsm_vrt, out_dsm_tif))

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


# ALL_STEPS is a ordonned dictionary : key = 'stepname' : value = is_distributed (True/False)
# initialization : pass in a sequence of tuples
ALL_STEPS = [('initialisation', False),
             ('local-pointing', True),
             ('global-pointing', False),
             ('rectification', True),
             ('matching', True),
             ('triangulation', True),
             ('disparity-to-height', True),
             ('global-mean-heights', False),
             ('heights-to-ply', True),
             ('global-srcwin', False),
             ('local-dsm-rasterization', True),
             ('global-dsm-rasterization', False),
             ('lidar-preprocessor', False)]
ALL_STEPS = collections.OrderedDict(ALL_STEPS)


def main(user_cfg, steps=ALL_STEPS):
    """
    Launch the s2p pipeline with the parameters given in a json file.

    Args:
        user_cfg: user config dictionary
        steps: either a string (single step) or a list of strings (several
            steps)
    """
    common.print_elapsed_time.t0 = datetime.datetime.now()
    initialization.build_cfg(user_cfg)
    if 'initialisation' in steps:
        initialization.make_dirs()

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes']:
        nb_workers = min(nb_workers, cfg['max_processes'])
    cfg['max_processes'] = nb_workers

    tw, th = initialization.adjust_tile_size()
    print('\ndiscarding masked tiles...')
    tiles = initialization.tiles_full_info(tw, th)

    if 'initialisation' in steps:
        # Write the list of json files to outdir/tiles.txt
        with open(os.path.join(cfg['out_dir'],'tiles.txt'),'w') as f:
            for t in tiles:
                f.write(t['json']+os.linesep)

    n = len(cfg['images'])
    tiles_pairs = [(t, i) for i in range(1, n) for t in tiles]

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

    if n > 2:
        if 'disparity-to-height' in steps:
            print('computing height maps...')
            parallel.launch_calls(disparity_to_height, tiles_pairs, nb_workers)

            print('computing local pairwise height offsets...')
            parallel.launch_calls(mean_heights, tiles, nb_workers)

        if 'global-mean-heights' in steps:
            print('computing global pairwise height offsets...')
            global_mean_heights(tiles)

        if 'heights-to-ply' in steps:
            print('merging height maps and computing point clouds...')
            parallel.launch_calls(heights_to_ply, tiles, nb_workers)

    else:
        if 'triangulation' in steps:
            print('triangulating tiles...')
            parallel.launch_calls(disparity_to_ply, tiles, nb_workers)

    if 'global-srcwin' in steps:
        print('computing global source window (xoff, yoff, xsize, ysize)...')
        global_srcwin(tiles)
        common.print_elapsed_time()

    if 'local-dsm-rasterization' in steps:
        print('computing DSM by tile...')
        parallel.launch_calls(plys_to_dsm, tiles, nb_workers)

    if 'global-dsm-rasterization' in steps:
        print('computing global DSM...')
        global_dsm(tiles)
        common.print_elapsed_time()

    if 'lidar-preprocessor' in steps:
        if cfg['run_lidar_preprocessor']:
            print('lidar preprocessor...')
            lidar_preprocessor(tiles)
            common.print_elapsed_time()
        else:
            print("LidarPreprocessor explicitly disabled in config.json")

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

    # read the json configuration file
    with open(args.config, 'r') as f:
        user_cfg = json.load(f)

    main(user_cfg, args.step)
