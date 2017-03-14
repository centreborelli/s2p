# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from __future__ import print_function
import os
import sys
import utm
import json
import copy
import shutil
import numpy as np

from s2plib import piio
from s2plib import common
from s2plib import srtm
from s2plib import rpc_utils
from s2plib import masking
from s2plib import parallel
from s2plib.config import cfg


def dict_has_keys(d, l):
    """
    Return True if the dict d contains all the keys of the input list l.
    """
    return all(k in d for k in l)


def check_parameters(d):
    """
    Check that the provided dictionary defines all mandatory s2p arguments.

    Args:
        d: python dictionary
    """
    # verify that input files paths are defined
    if 'images' not in d or len(d['images']) < 2:
        print('ERROR: missing paths to input images')
        sys.exit(1)
    for img in d['images']:
        if not dict_has_keys(img, ['img', 'rpc']):
            print('ERROR: missing img or rpc paths for image', img)
            sys.exit(1)

    # verify that roi or path to preview file are defined
    if 'full_img' in d and d['full_img']:
        sz = common.image_size_gdal(d['images'][0]['img'])
        d['roi'] = {'x': 0, 'y': 0, 'w': sz[0], 'h': sz[1]}
    elif 'roi' in d and dict_has_keys(d['roi'], ['x', 'y', 'w', 'h']):
        pass
    elif 'roi_utm' in d and dict_has_keys(d['roi_utm'], ['utm_band',
                                                         'hemisphere',
                                                         'x', 'y', 'w', 'h']):
        d['roi'] = rpc_utils.utm_roi_to_img_roi(d['images'][0]['rpc'],
                                                d['roi_utm'])
    elif 'roi_kml' in d:
        # this call defines cfg['ll_bbx'] and cfg['utm_bbx'] as side effects
        d['roi'] = rpc_utils.kml_roi_process(d['images'][0]['rpc'],
                                             d['roi_kml'])
    elif 'prv' in d['images'][0]:
        x, y, w, h = common.get_roi_coordinates(d['images'][0]['img'],
                                                d['images'][0]['prv'])
        d['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
    else:
        print('ERROR: missing or incomplete roi definition')
        sys.exit(1)

    # warn about unknown parameters. The known parameters are those defined in
    # the global config.cfg dictionary, plus the mandatory 'images' and 'roi' or
    # 'roi_utm'
    for k in d.keys():
        if k not in ['images', 'roi', 'roi_kml', 'roi_utm', 'utm_zone']:
            if k not in cfg:
                print('WARNING: ignoring unknown parameter {}.'.format(k))


def build_cfg(user_cfg):
    """
    Populate a dictionary containing the s2p parameters from a user config file.

    This dictionary is contained in the global variable 'cfg' of the config
    module.

    Args:
        user_cfg: user config dictionary
    """
    # check that all the mandatory arguments are defined
    check_parameters(user_cfg)

    # fill the config module: updates the content of the config.cfg dictionary
    # with the content of the user_cfg dictionary
    cfg.update(user_cfg)

    # sets keys 'clr', 'cld' and 'roi' of the reference image to None if they
    # are not already defined. The default values of these optional arguments
    # can not be defined directly in the config.py module. They would be
    # overwritten by the previous update, because they are in a nested dict.
    cfg['images'][0].setdefault('clr')
    cfg['images'][0].setdefault('cld')
    cfg['images'][0].setdefault('roi')
    cfg['images'][0].setdefault('wat')

    # check the zoom factor
    z = cfg['subsampling_factor']
    assert(z > 0 and z == np.floor(z))

    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    x = cfg['roi']['x']
    y = cfg['roi']['y']
    w = cfg['roi']['w']
    h = cfg['roi']['h']
    x, y, w, h = common.round_roi_to_nearest_multiple(z, x,y, w, h)
    cfg['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}

    # if srtm is disabled set disparity range method to sift
    if 'disable_srtm' in cfg and cfg['disable_srtm']:
        cfg['disp_range_method'] = 'sift'

    # get utm zone
    if 'utm_zone' not in cfg or cfg['utm_zone'] is None:
        cfg['utm_zone'] = rpc_utils.utm_zone(cfg['images'][0]['rpc'], x, y, w, h)


def make_dirs():
    """
    Create directories needed to run s2p.
    """
    common.mkdir_p(cfg['out_dir'])
    common.mkdir_p(os.path.expandvars(cfg['temporary_dir']))

    # store a json dump of the config.cfg dictionary
    with open(os.path.join(cfg['out_dir'], 'config.json'), 'w') as f:
        json.dump(cfg, f, indent=2)

    # copy RPC xml files in the output directory
    for img in cfg['images']:
        shutil.copy2(img['rpc'], cfg['out_dir'])

    # download needed srtm tiles
    if not cfg['disable_srtm']:
        x = cfg['roi']['x']
        y = cfg['roi']['y']
        w = cfg['roi']['w']
        h = cfg['roi']['h']
        for s in srtm.list_srtm_tiles(cfg['images'][0]['rpc'], x, y, w, h):
            srtm.get_srtm_tile(s, cfg['srtm_dir'])


def adjust_tile_size():
    """
    Adjust the size of the tiles.
    """
    zoom = cfg['subsampling_factor']
    tile_w = min(cfg['roi']['w'], zoom * cfg['tile_size'])  # tile width
    ntx = int(np.round(float(cfg['roi']['w']) / tile_w))
    # ceil so that, if needed, the last tile is slightly smaller
    tile_w = int(np.ceil(float(cfg['roi']['w']) / ntx))

    tile_h = min(cfg['roi']['h'], zoom * cfg['tile_size'])  # tile height
    nty = int(np.round(float(cfg['roi']['h']) / tile_h))
    tile_h = int(np.ceil(float(cfg['roi']['h']) / nty))

    print('tile size: {} {}'.format(tile_w, tile_h))
    n = len(cfg['images'])
    if n == 2:
        print('total number of tiles: {} ({} x {})'.format(ntx * nty, ntx, nty))
    else:
        print('total number of tiles: {} ({} x {}) x {} pairs'.format(ntx*nty*(n-1),
                                                                      ntx, nty, n-1))
    return tile_w, tile_h


def compute_tiles_coordinates(rx, ry, rw, rh, tw, th, z=1):
    """
    """
    out = []
    neighborhood_dict = dict()

    for y in np.arange(ry, ry + rh, th):
        h = min(th, ry + rh - y)
        for x in np.arange(rx, rx + rw, tw):
            w = min(tw, rx + rw - x)

            # ensure that tile coordinates are multiples of the zoom factor
            x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)

            out.append((x, y, w, h))

            # get coordinates of tiles from neighborhood
            out2 = []
            for y2 in [y - th, y, y + th]:
                h2 = min(th, ry + rh - y2)
                for x2 in [x - tw, x, x + tw]:
                    w2 = min(tw, rx + rw - x2)
                    if rx + rw > x2 >= rx:
                        if ry + rh > y2 >= ry:
                            x2, y2, w2, h2 = common.round_roi_to_nearest_multiple(z, x2, y2, w2, h2)
                            out2.append((x2, y2, w2, h2))

            neighborhood_dict[str((x, y, w, h))] = out2

    return out, neighborhood_dict


def get_tile_dir(x, y, w, h):
    """
    Get the name of a tile directory
    """
    return os.path.join(cfg['out_dir'],
                        'tiles_row_{}_height_{}'.format(y, h),
                        'col_{}_width_{}'.format(x, w))


def tiles_full_info(tw, th):
    """
    List the tiles to process and prepare their output directories structures.

    Most of the time is spent discarding tiles that are masked by water
    (according to SRTM data).

    Returns:
        a list of dictionaries. Each dictionary contains the image coordinates
        and the output directory path of a tile.
    """
    rpc = cfg['images'][0]['rpc']
    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    z =  cfg['subsampling_factor']
    rx = cfg['roi']['x']
    ry = cfg['roi']['y']
    rw = cfg['roi']['w']
    rh = cfg['roi']['h']

    # list tiles coordinates
    tiles_coords, neighborhood_coords_dict = compute_tiles_coordinates(rx, ry, rw, rh, tw, th, z)

    # compute all masks in parallel as numpy arrays
    tiles_masks = parallel.launch_calls_simple(masking.cloud_water_image_domain,
                                               tiles_coords,
                                               cfg['max_processes'], rpc,
                                               roi_msk, cld_msk, wat_msk,
                                               cfg['use_srtm_for_water'])

    # build a tile dictionary for all non-masked tiles and store them in a list
    tiles = []
    for coords, mask in zip(tiles_coords,
                            tiles_masks):
        if mask.any():  # there's at least one non-masked pixel in the tile
            tile = {}
            x, y, w, h = coords
            tile['dir'] = get_tile_dir(x, y, w, h)
            tile['coordinates'] = coords
            tile['mask'] = mask
            tile['neighborhood_dirs'] = list()
            key = str((x, y, w, h))

            if 'neighborhood_dirs' in cfg:
                tile['neighborhood_dirs'] = cfg['neighborhood_dirs']
            elif key in neighborhood_coords_dict:
                for coords2 in neighborhood_coords_dict[key]:
                    x2, y2, w2, h2 = coords2
                    tile['neighborhood_dirs'].append(get_tile_dir(x2,
                                                                  y2,
                                                                  w2,
                                                                  h2))
            tiles.append(tile)

    # make tiles directories and store json configuration dumps
    for tile in tiles:
        common.mkdir_p(tile['dir'])
        for i in range(1, len(cfg['images'])):
            common.mkdir_p(os.path.join(tile['dir'], 'pair_{}'.format(i)))

        # save a json dump of the tile configuration
        tile_cfg = copy.deepcopy(cfg)
        x, y, w, h = tile['coordinates']
        tile_cfg['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
        tile_cfg['full_img'] = False
        tile_cfg['max_processes'] = 1
        tile_cfg['omp_num_threads'] = 1
        tile_cfg['neighborhood_dirs'] = tile['neighborhood_dirs']

        tile_json = os.path.join(tile['dir'], 'config.json')
        tile['json'] = tile_json

        with open(tile_json, 'w') as f:
            json.dump(tile_cfg, f, indent=2)

        # save the mask
        piio.write(os.path.join(tile['dir'],
                                'cloud_water_image_domain_mask.png'),
                   tile['mask'].astype(np.uint8))

    return tiles
