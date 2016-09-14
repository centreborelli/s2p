# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import utm
import json
import shutil
import numpy as np

from python import common
from python import srtm
from python import tee
from python import rpc_utils
from python import masking
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
    elif 'roi_kml' in d:
        # this call defines cfg['ll_bbx'] and cfg['utm_bbx'] as side effects
        d['roi'] = rpc_utils.kml_roi_process(d['images'][0]['rpc'],
                                             d['roi_kml'])
    elif 'prv' in d['images'][0]:
        x, y, w, h = common.get_roi_coordinates(d['images'][0]['img'],
                                                d['images'][0]['prv'])
        d['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
    else:
        print 'ERROR: missing or incomplete roi definition'
        sys.exit(1)

    # if srtm is disable set disparity range method to sift
    if 'disable_srtm' in d:
        d['disp_range_method'] = 'sift'

    # warn about unknown parameters. The known parameters are those defined in
    # the global config.cfg dictionary, plus the mandatory 'images' and 'roi' or
    # 'roi_utm'
    for k in d.keys():
        if k not in ['images', 'roi', 'roi_kml', 'roi_utm']:
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
    if cfg.has_key('roi_utm'):
        cfg['utm_zone'] = cfg['roi_utm']['utm_band']
    elif cfg.has_key('roi_kml'):
        lon = (cfg['ll_bbx'][0] + cfg['ll_bbx'][1]) * .5
        lat = (cfg['ll_bbx'][2] + cfg['ll_bbx'][3]) * .5
        z = utm.conversion.latlon_to_zone_number(lat, lon)
        cfg['utm_zone'] = '%dS' % z if lat < 0 else '%dN' % z
    else:
        cfg['utm_zone'] = rpc_utils.utm_zone(cfg['images'][0]['rpc'], x, y, w, h)


def make_dirs():
    """
    Create directories needed to run s2p.
    """
    common.mkdir_p(cfg['out_dir'])
    common.mkdir_p(cfg['temporary_dir'])
    common.mkdir_p(os.path.join(cfg['temporary_dir'], 'meta'))

    # store a json dump of the config.cfg dictionary
    f = open(os.path.join(cfg['out_dir'], 'config.json'), 'w')
    json.dump(cfg, f, indent=2)
    f.close()

    # duplicate stdout and stderr to log file
    tee.Tee(os.path.join(cfg['out_dir'], 'stdout.log'), 'w')

    # download needed srtm tiles
    if not cfg['disable_srtm']:
        for s in srtm.list_srtm_tiles(cfg['images'][0]['rpc'],
                                      *cfg['roi'].values()):
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

    print 'tile size: {} {}'.format(tile_w, tile_h)
    n = len(cfg['images'])
    if n == 2:
        print 'total number of tiles: {} ({} x {})'.format(ntx * nty, ntx, nty)
    else:
        print 'total number of tiles: {} ({} x {} x {})'.format(ntx*nty*(n-1),
                                                                ntx, nty, n-1)
    return tile_w, tile_h


def tiles_full_info():
    """
    Prepare the entire process.

    Build tiles_full_info: a list of dictionaries, one per tile, providing all you need to process a tile
       * col/row : position of the tile (upper left corner)
       * tw/th : size of the tile

    Returns:
        tiles_full_info: list containing dictionaries
    """
    rpc = cfg['images'][0]['rpc']
    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    z =  cfg['subsampling_factor']
    rx, ry, rw, rh = cfg['roi'].values()  # roi coordinates
    tw, th = adjust_tile_size()  # default tile size

    # build tile_info dictionaries and store them in a list
    tiles_full_info = list()
    for y in np.arange(ry, ry + rh, th):
        h = min(th, ry + rh - y)
        for x in np.arange(rx, rx + rw, tw):
            w = min(tw, rx + rw - x)

            # ensure that tile coordinates are multiples of the zoom factor
            x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)

            # check if the tile is entirely masked by water, clouds, or img mask
            H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
            msk = common.tmpfile('.png')
            if masking.cloud_water_image_domain(msk, w, h, H, rpc, roi_msk, cld_msk, wat_msk):
                print 'tile {} {} {} {} masked, we skip it'.format(x, y, w, h)
                continue

            # add the tile to the list
            tile_info = {}
            tile_info['directory'] = os.path.join(cfg['out_dir'],
                                                  'tiles_row_%d_height_%d' % (y,
                                                                              h),
                                                  'col_%d_width_%d' % (x, w))
            tile_info['coordinates'] = (x, y, w, h)
            tiles_full_info.append(tile_info)

            # make the directories
            common.mkdir_p(tile_info['directory'])
            if len(cfg['images']) > 2:
                for i in xrange(1, len(cfg['images'])):
                    common.mkdir_p(os.path.join(tile_info['directory'],
                                                'pair_{}'.format(i)))

            # keep the mask
            shutil.copy(msk, os.path.join(tile_info['directory'],
                                          'cloud_water_image_domain_mask.png'))

    return tiles_full_info
