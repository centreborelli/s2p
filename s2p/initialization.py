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
import warnings
import numpy as np

from s2p import common
from s2p import rpc_utils
from s2p import masking
from s2p import parallel
from s2p.config import cfg

# This function is here as a workaround to python bug #24313 When
# using python3, json does not know how to serialize numpy.int64 on
# some platform numpy also decides to go for int64 when numpy.arange
# is called. This results in our json not being serializable anymore
# Calling json.dump(..,default=workaround_json_int64) fixes this
# https://bugs.python.org/issue24313
def workaround_json_int64(o):
    if isinstance(o, np.integer) : return int(o)
    raise TypeError

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
    for i in range(len(d['images'])):
        img = d['images'][i]
        if not dict_has_keys(img, ['img']):
            print('ERROR: missing img paths for image', img)
            sys.exit(1)
        if not dict_has_keys(img, ['rpc']) or img['rpc'] == '':
            img['rpc'] = rpc_utils.rpc_from_geotiff(img['img'])

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
    elif 'roi_geojson' in d:
        # this call defines cfg['ll_bbx'] and cfg['utm_bbx'] as side effects
        d['roi'] = rpc_utils.geojson_roi_process(d['images'][0]['rpc'],
                                                 d['roi_geojson'])
    elif 'prv' in d['images'][0]:
        x, y, w, h = common.get_roi_coordinates(d['images'][0]['img'],
                                                d['images'][0]['prv'])
        d['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
    else:
        print('ERROR: missing or incomplete roi definition')
        sys.exit(1)

    # d['roi'] : all the values must be integers
    d['roi']['x'] = int(np.floor(d['roi']['x']))
    d['roi']['y'] = int(np.floor(d['roi']['y']))
    d['roi']['w'] = int(np.ceil(d['roi']['w']))
    d['roi']['h'] = int(np.ceil(d['roi']['h']))

    # warn about unknown parameters. The known parameters are those defined in
    # the global config.cfg dictionary, plus the mandatory 'images' and 'roi' or
    # 'roi_utm'
    for k in d.keys():
        if k not in ['images', 'roi', 'roi_kml', 'roi_geojson', 'roi_utm', 'utm_zone']:
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

    # Make sure that input data have absolute paths
    for i in range(len(cfg['images'])):
        for d in ['clr', 'cld', 'roi', 'wat', 'img']:
            if d in cfg['images'][i] and cfg['images'][i][d] is not None and not os.path.isabs(cfg['images'][i][d]):
                cfg['images'][i][d] = os.path.abspath(cfg['images'][i][d])

    x = cfg['roi']['x']
    y = cfg['roi']['y']
    w = cfg['roi']['w']
    h = cfg['roi']['h']

    cfg['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}

    # get utm zone
    if 'utm_zone' not in cfg or cfg['utm_zone'] is None:
        cfg['utm_zone'] = rpc_utils.utm_zone(cfg['images'][0]['img'], x, y, w, h)


def make_dirs():
    """
    Create directories needed to run s2p.
    """
    common.mkdir_p(cfg['out_dir'])
    common.mkdir_p(os.path.expandvars(cfg['temporary_dir']))

    # store a json dump of the config.cfg dictionary
    with open(os.path.join(cfg['out_dir'], 'config.json'), 'w') as f:
        cfg_copy = copy.deepcopy(cfg)
        cfg_copy['out_dir'] = '.'
        for img in cfg_copy['images']:
            img.pop('rpc', None)
        json.dump(cfg_copy, f, indent=2, default=workaround_json_int64)


def adjust_tile_size():
    """
    Adjust the size of the tiles.
    """

    tile_w = min(cfg['roi']['w'], cfg['tile_size'])  # tile width
    ntx = int(np.round(float(cfg['roi']['w']) / tile_w))
    # ceil so that, if needed, the last tile is slightly smaller
    tile_w = int(np.ceil(float(cfg['roi']['w']) / ntx))

    tile_h = min(cfg['roi']['h'], cfg['tile_size'])  # tile height
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


def compute_tiles_coordinates(rx, ry, rw, rh, tw, th):
    """
    """
    out = []
    neighborhood_dict = dict()

    for y in np.arange(ry, ry + rh, th):
        h = min(th, ry + rh - y)
        for x in np.arange(rx, rx + rw, tw):
            w = min(tw, rx + rw - x)

            out.append((x, y, w, h))

            # get coordinates of tiles from neighborhood
            out2 = []
            for y2 in [y - th, y, y + th]:
                h2 = min(th, ry + rh - y2)
                for x2 in [x - tw, x, x + tw]:
                    w2 = min(tw, rx + rw - x2)
                    if rx + rw > x2 >= rx:
                        if ry + rh > y2 >= ry:
                            out2.append((x2, y2, w2, h2))

            neighborhood_dict[str((x, y, w, h))] = out2

    return out, neighborhood_dict


def get_tile_dir(x, y, w, h):
    """
    Get the name of a tile directory
    """
    return os.path.join('tiles','row_{:07d}_height_{}'.format(y, h),
                        'col_{:07d}_width_{}'.format(x, w))


def create_tile(coords, neighborhood_coords_dict):
    tile = {}
    x, y, w, h = coords
    tile['dir'] = os.path.join(cfg['out_dir'], get_tile_dir(x, y, w, h))
    tile['coordinates'] = coords
    tile['neighborhood_dirs'] = list()
    key = str((x, y, w, h))

    if 'neighborhood_dirs' in cfg:
        tile['neighborhood_dirs'] = cfg['neighborhood_dirs']
    elif key in neighborhood_coords_dict:
        for coords2 in neighborhood_coords_dict[key]:
            x2, y2, w2, h2 = coords2
            tile['neighborhood_dirs'].append(os.path.join('../../..', get_tile_dir(x2,
                                                                                   y2,
                                                                                   w2,
                                                                                   h2)))

    tile_json = os.path.join(get_tile_dir(x,y,w,h),'config.json')
    tile['json'] = tile_json

    return tile

def tiles_full_info(tw, th, tiles_txt, create_masks=False):
    """
    List the tiles to process and prepare their output directories structures.

    Most of the time is spent discarding tiles that are masked by water
    (according to exogenous dem).

    Returns:
        a list of dictionaries. Each dictionary contains the image coordinates
        and the output directory path of a tile.
    """
    rpc = cfg['images'][0]['rpc']
    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    rx = cfg['roi']['x']
    ry = cfg['roi']['y']
    rw = cfg['roi']['w']
    rh = cfg['roi']['h']

    # build a tile dictionary for all non-masked tiles and store them in a list
    tiles = []
    # list tiles coordinates
    tiles_coords, neighborhood_coords_dict = compute_tiles_coordinates(rx, ry, rw, rh, tw, th)

    if os.path.exists(tiles_txt) is False or create_masks is True:
        print('\ndiscarding masked tiles...')
        # compute all masks in parallel as numpy arrays
        tiles_masks = parallel.launch_calls_simple(masking.cloud_water_image_domain,
                                                   tiles_coords,
                                                   cfg['max_processes'], rpc,
                                                   roi_msk, cld_msk, wat_msk)

        for coords, mask in zip(tiles_coords,
                                tiles_masks):
            if mask.any():  # there's at least one non-masked pixel in the tile
                tile = create_tile(coords, neighborhood_coords_dict)
                tiles.append(tile)

                # make tiles directories and store json configuration dumps
                common.mkdir_p(tile['dir'])
                for i in range(1, len(cfg['images'])):
                    common.mkdir_p(os.path.join(tile['dir'], 'pair_{}'.format(i)))

                # save a json dump of the tile configuration
                tile_cfg = copy.deepcopy(cfg)
                x, y, w, h = tile['coordinates']
                for img in tile_cfg['images']:
                    img.pop('rpc', None)
                tile_cfg['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
                tile_cfg['full_img'] = False
                tile_cfg['max_processes'] = 1
                tile_cfg['neighborhood_dirs'] = tile['neighborhood_dirs']
                tile_cfg['out_dir'] = '../../..'

                with open(os.path.join(cfg['out_dir'], tile['json']), 'w') as f:
                    json.dump(tile_cfg, f, indent=2,default=workaround_json_int64)

                # save the mask
                common.rasterio_write(os.path.join(tile['dir'],
                                                   'cloud_water_image_domain_mask.png'),
                                      mask.astype(np.uint8))
    else:
        if len(tiles_coords) == 1:
            tiles.append(create_tile(tiles_coords[0], neighborhood_coords_dict))
        else:
            with open(tiles_txt, 'r') as f_tiles:
                for config_json in f_tiles:
                    tile = {}
                    with open(os.path.join(cfg['out_dir'],
                                           config_json.rstrip(os.linesep)), 'r') as f_config:
                        tile_cfg = json.load(f_config)
                        roi = tile_cfg['roi']
                        coords = roi['x'], roi['y'], roi['w'], roi['h']
                        tiles.append(create_tile(coords, neighborhood_coords_dict))

    return tiles
