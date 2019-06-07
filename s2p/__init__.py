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
import collections
import shutil
import rasterio


from s2p.config import cfg
from s2p import common
from s2p import parallel
from s2p import initialization
from s2p import pointing_accuracy
from s2p import rectification
from s2p import block_matching
from s2p import masking
from s2p import triangulation
from s2p import fusion
from s2p import rasterization
from s2p import visualisation


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

    # correct pointing error
    print('correcting pointing on tile {} {} pair {}...'.format(x, y, i))
    try:
        A, m = pointing_accuracy.compute_correction(img1, img2, x, y, w, h)
    except common.RunFailure as e:
        stderr = os.path.join(out_dir, 'stderr.log')
        with open(stderr, 'w') as f:
            f.write('ERROR during pointing correction with cmd: %s\n' % e[0]['command'])
            f.write('Stop processing this pair\n')
        return

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

    if os.path.exists(os.path.join(out_dir, 'stderr.log')):
        print('rectification: stderr.log exists')
        print('pair_{} not processed on tile {} {}'.format(i, x, y))
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

    cur_dir = os.path.join(tile['dir'],'pair_{}'.format(i))
    for n in tile['neighborhood_dirs']:
        nei_dir = os.path.join(tile['dir'], n, 'pair_{}'.format(i))
        if os.path.exists(nei_dir) and not os.path.samefile(cur_dir, nei_dir):
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

    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2,
                                                            x, y, w, h,
                                                            rect1, rect2, A, m,
                                                            method=cfg['rectification_method'],
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

    if os.path.exists(os.path.join(out_dir, 'stderr.log')):
        print('disparity estimation: stderr.log exists')
        print('pair_{} not processed on tile {} {}'.format(i, x, y))
        return

    print('estimating disparity on tile {} {} pair {}...'.format(x, y, i))
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    disp_min, disp_max = np.loadtxt(os.path.join(out_dir, 'disp_min_max.txt'))

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

    if os.path.exists(os.path.join(out_dir, 'stderr.log')):
        print('triangulation: stderr.log exists')
        print('pair_{} not processed on tile {} {}'.format(i, x, y))
        return

    print('triangulating tile {} {} pair {}...'.format(x, y, i))
    rpc1 = cfg['images'][0]['rpc']
    rpc2 = cfg['images'][i]['rpc']
    H_ref = os.path.join(out_dir, 'H_ref.txt')
    H_sec = os.path.join(out_dir, 'H_sec.txt')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    rpc_err = os.path.join(out_dir, 'rpc_err.tif')
    out_mask = os.path.join(tile['dir'], 'mask.png')
    pointing = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_{}.txt'.format(i))
    triangulation.height_map(height_map, x, y, w, h,
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

    if os.path.exists(os.path.join(out_dir, 'stderr.log')):
        print('triangulation: stderr.log exists')
        print('pair_1 not processed on tile {} {}'.format(x, y))
        return

    print('triangulating tile {} {}...'.format(x, y))
    # This function is only called when there is a single pair (pair_1)
    H_ref = os.path.join(out_dir, 'pair_1', 'H_ref.txt')
    H_sec = os.path.join(out_dir, 'pair_1', 'H_sec.txt')
    pointing = os.path.join(cfg['out_dir'], 'global_pointing_pair_1.txt')
    disp  = os.path.join(out_dir, 'pair_1', 'rectified_disp.tif')
    extra = os.path.join(out_dir, 'pair_1', 'rectified_disp_confidence.tif')
    if not os.path.exists(extra):
        extra = ''
    mask_rect = os.path.join(out_dir, 'pair_1', 'rectified_mask.png')
    mask_orig = os.path.join(out_dir, 'mask.png')

    # prepare the image needed to colorize point cloud
    colors = os.path.join(out_dir, 'rectified_ref.png')
    if cfg['images'][0]['clr']:
        hom = np.loadtxt(H_ref)
        # We want rectified_ref.png and rectified_ref.tif to have the same size
        with rasterio.open(os.path.join(out_dir, 'pair_1', 'rectified_ref.tif')) as f:
            ww, hh = f.width, f.height
        common.image_apply_homography(colors, cfg['images'][0]['clr'], hom, ww, hh)
    else:
        common.image_qauto(os.path.join(out_dir, 'pair_1', 'rectified_ref.tif'), colors)

    # compute the point cloud
    triangulation.disp_map_to_point_cloud(ply_file, disp, mask_rect, rpc1, rpc2,
                                          H_ref, H_sec, pointing, colors, extra,
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
    n = len(cfg['images']) - 1
    maps = np.empty((h, w, n))
    for i in range(n):
        try:
            with rasterio.open(os.path.join(tile['dir'], 'pair_{}'.format(i + 1),
                                            'height_map.tif'), 'r') as f:
                maps[:, :, i] = f.read(1)
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
    plyfile = os.path.join(out_dir, 'cloud.ply')
    plyextrema = os.path.join(out_dir, 'plyextrema.txt')
    height_map = os.path.join(out_dir, 'height_map.tif')

    # H is the homography transforming the coordinates system of the original
    # full size image into the coordinates system of the crop
    H = np.dot(np.diag([1, 1, 1]), common.matrix_translation(-x, -y))
    colors = os.path.join(out_dir, 'ref.png')
    if cfg['images'][0]['clr']:
        common.image_crop_gdal(cfg['images'][0]['clr'], x, y, w, h, colors)
    else:
        common.image_qauto(common.image_crop_gdal(cfg['images'][0]['img'], x, y,
                                                 w, h), colors)

    triangulation.height_map_to_point_cloud(plyfile, height_map,
                                            cfg['images'][0]['rpc'], H, colors,
                                            utm_zone=cfg['utm_zone'],
                                            llbbx=tuple(cfg['ll_bbx']))

    # compute the point cloud extrema (xmin, xmax, xmin, ymax)
    common.run("plyextrema %s %s" % (plyfile, plyextrema))

    if cfg['clean_intermediate']:
        common.remove(height_map)
        common.remove(colors)
        common.remove(os.path.join(out_dir, 'mask.png'))

def plys_to_dsm(tile):
    """
    Generates DSM from plyfiles (cloud.ply)

    Args:
        tile: a dictionary that provides all you need to process a tile
    """
    out_dsm  = os.path.join(tile['dir'], 'dsm.tif')
    out_conf = os.path.join(tile['dir'], 'confidence.tif')

    res = cfg['dsm_resolution']
    if 'utm_bbx' in cfg:
        bbx = cfg['utm_bbx']
        global_xoff = bbx[0]
        global_yoff = bbx[3]
    else:
        global_xoff = 0 # arbitrary reference
        global_yoff = 0 # arbitrary reference

    res = cfg['dsm_resolution']

    xmin, xmax, ymin, ymax = np.loadtxt(os.path.join(tile['dir'], "plyextrema.txt"))

    # compute xoff, yoff, xsize, ysize considering final dsm
    local_xoff = global_xoff + np.floor((xmin - global_xoff) / res) * res
    local_xsize = int(1 + np.floor((xmax - local_xoff) / res))

    local_yoff = global_yoff + np.ceil((ymax - global_yoff) / res) * res
    local_ysize = int(1 - np.floor((ymin - local_yoff) / res))

    clouds = [os.path.join(tile['dir'],n_dir, 'cloud.ply') for n_dir in tile['neighborhood_dirs']]
    roi = (local_xoff, local_yoff, local_xsize, local_ysize)

    raster, profile = rasterization.plyflatten_from_plyfiles_list(clouds,
                                                                  resolution=res,
                                                                  roi=roi,
                                                                  radius=cfg['dsm_radius'],
                                                                  sigma=cfg['dsm_sigma'])

    # save_output_image_with_utm_georeferencing
    common.rasterio_write(out_dsm, raster[:, :, 0], profile=profile)

    # export confidence (optional)
    if raster.shape[-1] == 5:
        common.rasterio_write(out_conf, raster[:, :, 4], profile=profile)


def global_dsm(tiles):
    """
    """
    out_dsm_vrt = os.path.join(cfg['out_dir'], 'dsm.vrt')
    out_dsm_tif = os.path.join(cfg['out_dir'], 'dsm.tif')

    dsms_list = [os.path.join(t['dir'], 'dsm.tif') for t in tiles]
    dsms = '\n'.join(d for d in dsms_list if os.path.exists(d))

    input_file_list = os.path.join(cfg['out_dir'], 'gdalbuildvrt_input_file_list.txt')

    with open(input_file_list, 'w') as f:
        f.write(dsms)

    common.run("gdalbuildvrt -vrtnodata nan -input_file_list %s %s" % (input_file_list,
                                                                       out_dsm_vrt))

    res = cfg['dsm_resolution']

    if 'utm_bbx' in cfg:
        bbx = cfg['utm_bbx']
        xoff = bbx[0]
        yoff = bbx[3]
        xsize = int(np.ceil((bbx[1]-bbx[0]) / res))
        ysize = int(np.ceil((bbx[3]-bbx[2]) / res))
        projwin = "-projwin %s %s %s %s" % (xoff, yoff,
                                            xoff + xsize * res,
                                            yoff - ysize * res)
    else:
        projwin = ""

    common.run(" ".join(["gdal_translate",
                         "-co TILED=YES -co BIGTIFF=IF_SAFER",
                         "%s %s %s" % (projwin, out_dsm_vrt, out_dsm_tif)]))

    # EXPORT CONFIDENCE
    out_conf_vrt = os.path.join(cfg['out_dir'], 'confidence.vrt')
    out_conf_tif = os.path.join(cfg['out_dir'], 'confidence.tif')

    dsms_list = [os.path.join(t['dir'], 'confidence.tif') for t in tiles]
    dems_list_ok = [d for d in dsms_list if os.path.exists(d)]
    dsms = '\n'.join(dems_list_ok)

    input_file_list = os.path.join(cfg['out_dir'], 'gdalbuildvrt_input_file_list2.txt')

    if len(dems_list_ok) > 0:

        with open(input_file_list, 'w') as f:
            f.write(dsms)

        common.run("gdalbuildvrt -vrtnodata nan -input_file_list %s %s" % (input_file_list,
                                                                           out_conf_vrt))

        common.run(" ".join(["gdal_translate",
                             "-co TILED=YES -co BIGTIFF=IF_SAFER",
                             "%s %s %s" % (projwin, out_conf_vrt, out_conf_tif)]))


def main(user_cfg):
    """
    Launch the s2p pipeline with the parameters given in a json file.

    Args:
        user_cfg: user config dictionary
    """
    common.print_elapsed_time.t0 = datetime.datetime.now()
    initialization.build_cfg(user_cfg)
    initialization.make_dirs()

    # multiprocessing setup
    nb_workers = multiprocessing.cpu_count()  # nb of available cores
    if cfg['max_processes'] is not None:
        nb_workers = cfg['max_processes']

    tw, th = initialization.adjust_tile_size()
    tiles_txt = os.path.join(cfg['out_dir'],'tiles.txt')
    tiles = initialization.tiles_full_info(tw, th, tiles_txt, create_masks=True)

    # initialisation step:
    # Write the list of json files to outdir/tiles.txt
    with open(tiles_txt,'w') as f:
        for t in tiles:
            f.write(t['json']+os.linesep)

    n = len(cfg['images'])
    tiles_pairs = [(t, i) for i in range(1, n) for t in tiles]

    # local-pointing step:
    print('correcting pointing locally...')
    parallel.launch_calls(pointing_correction, tiles_pairs, nb_workers)

    # global-pointing step:
    print('correcting pointing globally...')
    global_pointing_correction(tiles)
    common.print_elapsed_time()

    # rectification step:
    print('rectifying tiles...')
    parallel.launch_calls(rectification_pair, tiles_pairs, nb_workers)

    # matching step:
    print('running stereo matching...')
    parallel.launch_calls(stereo_matching, tiles_pairs, nb_workers)

    if n > 2:
        # disparity-to-height step:
        print('computing height maps...')
        parallel.launch_calls(disparity_to_height, tiles_pairs, nb_workers)

        print('computing local pairwise height offsets...')
        parallel.launch_calls(mean_heights, tiles, nb_workers)

        # global-mean-heights step:
        print('computing global pairwise height offsets...')
        global_mean_heights(tiles)

        # heights-to-ply step:
        print('merging height maps and computing point clouds...')
        parallel.launch_calls(heights_to_ply, tiles, nb_workers)
    else:
        # triangulation step:
        print('triangulating tiles...')
        parallel.launch_calls(disparity_to_ply, tiles, nb_workers)

    # local-dsm-rasterization step:
    print('computing DSM by tile...')
    parallel.launch_calls(plys_to_dsm, tiles, nb_workers)

    # global-dsm-rasterization step:
    print('computing global DSM...')
    global_dsm(tiles)
    common.print_elapsed_time()

    # cleanup
    common.garbage_cleanup()
    common.print_elapsed_time(since_first_call=True)


def make_path_relative_to_file(path, f):
    return os.path.join(os.path.abspath(os.path.dirname(f)), path)


def read_tiles(tiles_file):
    tiles = []
    outdir = os.path.dirname(tiles_file)

    with open(tiles_file) as f:
        tiles = f.readlines()

    # Strip trailing \n
    tiles = list(map(str.strip,tiles))
    tiles = [os.path.join(outdir, t) for t in tiles]

    return tiles


def read_config_file(config_file):
    """
    Read a json configuration file and interpret relative paths.

    If any input or output path is a relative path, it is interpreted as
    relative to the config_file location (and not relative to the current
    working directory). Absolute paths are left unchanged.
    """
    with open(config_file, 'r') as f:
        user_cfg = json.load(f)

    # output paths
    if not os.path.isabs(user_cfg['out_dir']):
        print('WARNING: out_dir is a relative path. It is interpreted with '
              'respect to {} location (not cwd)'.format(config_file))
        user_cfg['out_dir'] = make_path_relative_to_file(user_cfg['out_dir'],
                                                         config_file)
        print('out_dir is: {}'.format(user_cfg['out_dir']))

    # ROI paths
    for k in ["roi_kml", "roi_geojson"]:
        if k in user_cfg and isinstance(user_cfg[k], str) and not os.path.isabs(user_cfg[k]):
            user_cfg[k] = make_path_relative_to_file(user_cfg[k], config_file)

    # input paths
    for img in user_cfg['images']:
        for d in ['img', 'rpc', 'clr', 'cld', 'roi', 'wat']:
            if d in img and img[d] is not None and not os.path.isabs(img[d]):
                img[d] = make_path_relative_to_file(img[d], config_file)

    return user_cfg
