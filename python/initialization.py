# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import json
import numpy as np

from python import common
from python import srtm
from python import tee
from python import rpc_utils
from config import cfg


def dict_has_keys(d, l):
    """
    Return True if the dict d contains all the keys of the input list l.
    """
    return all(k in d for k in l)


def check_parameters(d):
    """
    Checks that the provided dictionary defines all the mandatory arguments.

    Args:
        d: python dictionary
    """
    # verify that input files paths are defined
    if 'images' not in d or len(d['images']) < 2:
        print 'ERROR: missing paths to input images'
        sys.exit(1)
    for img in d['images']:
        if not dict_has_keys(img, ['img', 'rpc']):
            print 'ERROR: missing img or rpc paths for image', img
            sys.exit(1)

    # verify that roi or path to preview file are defined
    if 'full_img' in d and d['full_img']:
        sz = common.image_size_tiffinfo(d['images'][0]['img'])
        d['roi'] = {'x': 0, 'y': 0, 'w': sz[0], 'h': sz[1]}
    elif 'roi' in d and dict_has_keys(d['roi'], ['x', 'y', 'w', 'h']):
        pass
    elif 'roi_utm' in d and dict_has_keys(d['roi_utm'], ['utm_band',
                                                         'hemisphere',
                                                         'x', 'y', 'w', 'h']):
        d['roi'] = rpc_utils.utm_roi_to_img_roi(d['images'][0]['rpc'],
                                                d['roi_utm'])
    elif 'prv' in d['images'][0]:
        x, y, w, h = common.get_roi_coordinates(d['images'][0]['img'],
                                                d['images'][0]['prv'])
        d['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
    else:
        print 'ERROR: missing or incomplete roi definition'
        sys.exit(1)

    # warn about unknown parameters. The known parameters are those defined in
    # the global config.cfg dictionary, plus the mandatory 'images' and 'roi' or
    # 'roi_utm'
    for k in d.keys():
        if k not in ['images', 'roi', 'roi_utm']:
            if k not in cfg:
                print 'WARNING: ignoring unknown parameter {}.'.format(k)


def build_cfg(config_file):
    """
    Update the dictionary containing all the parameters controlling s2p.

    This dictionary is contained in the global variable 'cfg' of the config
    module.

    Args:
        config_file: path to a json configuration file
    """
    # read the json configuration file
    f = open(config_file)
    user_cfg = json.load(f)
    f.close()

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
    x, y, w, h = common.round_roi_to_nearest_multiple(z, *cfg['roi'].values())
    cfg['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}

    # get utm zone
    cfg['utm_zone'] = rpc_utils.utm_zone(cfg['images'][0]['rpc'], x, y, w, h)


def make_dirs():
    """
    Create directories needed to run s2p.
    """
    if not os.path.exists(cfg['out_dir']):
        os.makedirs(cfg['out_dir'])

    if not os.path.exists(os.path.join(cfg['out_dir'], 'dsm')):
        os.makedirs(os.path.join(cfg['out_dir'], 'dsm'))

    if not os.path.exists(cfg['temporary_dir']):
        os.makedirs(cfg['temporary_dir'])

    if not os.path.exists(os.path.join(cfg['temporary_dir'], 'meta')):
        os.makedirs(os.path.join(cfg['temporary_dir'], 'meta'))

    # store a json dump of the config.cfg dictionary
    f = open(os.path.join(cfg['out_dir'], 'config.json'), 'w')
    json.dump(cfg, f, indent=2)
    f.close()

    # duplicate stdout and stderr to log file
    tee.Tee(os.path.join(cfg['out_dir'], 'stdout.log'), 'w')

    # download needed srtm tiles
    for s in srtm.list_srtm_tiles(cfg['images'][0]['rpc'],
                                  *cfg['roi'].values()):
        srtm.get_srtm_tile(s, cfg['srtm_dir'])


def adjust_tile_size():
    """
    Adjust the size of the tiles.
    """
    zoom = cfg['subsampling_factor']
    overlap = 100 * zoom  # overlap between tiles
    tile_w = min(cfg['roi']['w'], zoom * cfg['tile_size'])  # tile width
    tile_h = min(cfg['roi']['h'], zoom * cfg['tile_size'])  # tile height
    print 'tile size: {} {}'.format(tile_w, tile_h)

    ntx = int(np.ceil(float(cfg['roi']['w'] - overlap) / (tile_w - overlap)))
    nty = int(np.ceil(float(cfg['roi']['h'] - overlap) / (tile_h - overlap)))
    print 'total number of tiles: {} ({} x {})'.format(ntx * nty, ntx, nty)
    
    return tile_w, tile_h, overlap


def tiles_full_info():
    """
    Prepare the entire process.

    Build tiles_full_info: a list of dictionaries, one per tile, providing all you need to process a tile
       * col/row : position of the tile (upper left corner)
       * tw/th : size of the tile
       * ov : size of the overlapping
       * i/j : relative position of the tile
       * pos : position inside the ROI : UL for a tile place at th Upper Left corner, M for the ones placed in the middle, and so forth.
       * x/y/w/h : information about the ROI
       * images : a dictionary directly given by the json config file, that store the information about all the involved images, their rpc, and so forth.
       * nb_pairs : number of pairs
       * cld_msk/roi_msk : path to a gml file containing a cloud mask/ defining the area contained in the full image

    Returns:
        tiles_full_info: list containing dictionaries
    """
    tw, th, ov = adjust_tile_size()
    x, y, w, h = cfg['roi'].values()
    z =  cfg['subsampling_factor']

    # build tile_info dictionaries and store them in a list
    tiles_full_info = list()
    range_y = np.arange(y, y + h - ov, th - ov)
    range_x = np.arange(x, x + w - ov, tw - ov)
    rowmin, rowmax = range_y[0], range_y[-1]
    colmin, colmax = range_x[0], range_x[-1]
       
    for i, row in enumerate(range_y):
        for j, col in enumerate(range_x):
            # ensure that tile coordinates are multiples of the zoom factor
            col, row, tw, th = common.round_roi_to_nearest_multiple(z, col, row,
                                                                    tw, th)
            tile_dir = os.path.join(cfg['out_dir'], 'tile_%d_%d_row_%d' % (tw,
                                                                           th,
                                                                           row),
                                    'col_%d' % col)

            if row == rowmin and col == colmin:
                pos = 'UL'
            elif row == rowmin and col == colmax:
                pos = 'UR'
            elif row == rowmax and col == colmax:
                pos = 'BR'
            elif row == rowmax and col == colmin:
                pos = 'BL'
            elif row == rowmin and col > colmin:
                pos = 'U'
            elif col == colmin and row > rowmin:
                pos = 'L'
            elif row == rowmax and col > colmin:
                pos = 'B'
            elif col == colmax and row > rowmin:
                pos = 'R'
            else:
                pos = 'M'

            tile_info = {}
            tile_info['directory'] = tile_dir
            tile_info['coordinates'] = (col, row, tw, th)
            tile_info['index_in_roi'] = (i, j)
            tile_info['position_type'] = pos
            tile_info['roi_coordinates'] = (x, y, w, h)
            tile_info['overlap'] = ov
            tile_info['number_of_pairs'] = len(cfg['images']) - 1
            tile_info['images'] = cfg['images']
            tiles_full_info.append(tile_info)

    if len(tiles_full_info) == 1:
        tiles_full_info[0]['position_type'] = 'Single'

    return tiles_full_info
