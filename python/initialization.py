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
from config import cfg


def check_parameters(usr_cfg):
    """
    Checks that the provided dictionary defines all the mandatory
    arguments, and warns about unknown optional arguments.

    Args:
        usr_cfg: python dict read from the json input file
    """

    # verify that i/o files are defined
    if 'out_dir' not in usr_cfg:
        print "missing output dir: abort"
        sys.exit(1)
    if 'images' not in usr_cfg or len(usr_cfg['images']) < 2:
        print "missing input data paths: abort"
        sys.exit(1)
    for i in range(0, len(usr_cfg['images'])):
        if 'img' not in usr_cfg['images'][i] or 'rpc' not in usr_cfg['images'][i]:
            errMsg = 'missing input data paths for image ' + str(i) + ': abort'
            print errMsg
            sys.exit(1)

    # verify that roi or path to preview file are defined
    if ('full_img' not in usr_cfg) or (not usr_cfg['full_img']):
        if 'roi' not in usr_cfg or any(p not in usr_cfg['roi'] for p in ['x',
                                                                         'y',
                                                                         'w',
                                                                         'h']):
            if 'prv' not in usr_cfg['images'][0]:
                print """missing or incomplete roi definition, and no preview
                file is specified: abort"""
                sys.exit(1)

    # warn about unknown optional parameters: these parameters have no default
    # value in the global config.cfg dictionary, and thus they are not used
    # anywhere.  They may appear in the usr_cfg because of a typo.
    l = usr_cfg.keys()

    # remove mandatory parameters (they are not in config.cfg)
    l.remove('out_dir')
    l.remove('images')
    if 'roi' in l:
        l.remove('roi')

    # check
    for k in l:
        if k not in cfg:
            print """parameter %s unknown: you should remove it from the input
            json file. It will be ignored.""" % k


def init_roi(config_file):
    """
    1) Loads configuration file
    2) Checks parameters
    3) Selects the ROI
    4) Checks the zoom factor
    """
	
    # read the json configuration file
    f = open(config_file)
    user_cfg = json.load(f)
    f.close()

    # Check that all the mandatory arguments are defined, and warn about
    # 'unknown' params
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

    # update roi definition if the full_img flag is set to true
    if ('full_img' in cfg) and cfg['full_img']:
        sz = common.image_size_tiffinfo(cfg['images'][0]['img'])
        cfg['roi'] = {}
        cfg['roi']['x'] = 0
        cfg['roi']['y'] = 0
        cfg['roi']['w'] = sz[0]
        cfg['roi']['h'] = sz[1]

    # check that the roi is well defined
    if 'roi' not in cfg or any(p not in cfg['roi'] for p in ['x', 'y', 'w',
                                                             'h']):
        print "missing or incomplete ROI definition"
        print "ROI will be redefined by interactive selection"
        x, y, w, h = common.get_roi_coordinates(cfg['images'][0]['img'],
                                                cfg['images'][0]['prv'])
        cfg['roi'] = {}
        cfg['roi']['x'] = x
        cfg['roi']['y'] = y
        cfg['roi']['w'] = w
        cfg['roi']['h'] = h
    else :
        x = cfg['roi']['x']    
        y = cfg['roi']['y']
        w = cfg['roi']['w']
        h = cfg['roi']['h']

    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        print 'Neither a ROI nor a preview file are defined. Aborting.'
        return

    # check the zoom factor
    z = cfg['subsampling_factor']
    assert(z > 0 and z == np.floor(z))
       
    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.          
    x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)
    cfg['roi']['x'] = x    
    cfg['roi']['y'] = y
    cfg['roi']['w'] = w
    cfg['roi']['h'] = h


def init_dirs_srtm(config_file):
    """
    1) Creates different directories : output, temp...
    2) Downloads srtm files
    
    Args:
        - config_file : a json configuratio file
    """

    init_roi(config_file)


    # create tmp dir and output directory for the experiment, and store a json
    # dump of the config.cfg dictionary there, download srtm files...
    if not os.path.exists(cfg['out_dir']):
        os.makedirs(cfg['out_dir'])

        if not os.path.exists(cfg['temporary_dir']):
            os.makedirs(cfg['temporary_dir'])

        if not os.path.exists(os.path.join(cfg['temporary_dir'], 'meta')):
            os.makedirs(os.path.join(cfg['temporary_dir'], 'meta'))
        f = open('%s/config.json' % cfg['out_dir'], 'w')
        json.dump(cfg, f, indent=2)
        f.close()

        # duplicate stdout and stderr to log file
        tee.Tee('%s/stdout.log' % cfg['out_dir'], 'w')

        # needed srtm tiles
        srtm_tiles = srtm.list_srtm_tiles(cfg['images'][0]['rpc'],
                                          *cfg['roi'].values())
        for s in srtm_tiles:
            srtm.get_srtm_tile(s, cfg['srtm_dir'])


def init_tiles_full_info(config_file):
    """
    Prepare the entire process.

    1) Make sure coordinates of the ROI are multiples of the zoom factor
    2) Compute optimal size for tiles, get the number of pairs
    3) Build tiles_full_info: a list of dictionaries, one per tile, providing all you need to process a tile
       * col/row : position of the tile (upper left corner)
       * tw/th : size of the tile
       * ov : size of the overlapping
       * i/j : relative position of the tile
       * pos : position inside the ROI : UL for a tile place at th Upper Left corner, M for the ones placed in the middle, and so forth.
       * x/y/w/h : information about the ROI
       * images : a dictionary directly given by the json config file, that store the information about all the involved images, their rpc, and so forth.
       * nb_pairs : number of pairs
       * cld_msk/roi_msk : path to a gml file containing a cloud mask/ defining the area contained in the full image

    Args:
         config_file: path to a json configuration file

    Returns:
        tiles_full_info: list containing dictionaries
    """

    init_roi(config_file)

    #Get ROI
    x = cfg['roi']['x']    
    y = cfg['roi']['y']
    w = cfg['roi']['w']
    h = cfg['roi']['h']    
    z = cfg['subsampling_factor']

    # Automatically compute optimal size for tiles
    # tw, th : dimensions of the tiles
    # ov : width of overlapping bands between tiles
    ov = z * 100
    if w <= z * cfg['tile_size']:
        tw = w
    else:
        tw = z * cfg['tile_size']
    if h <= z * cfg['tile_size']:
        th = h
    else:
        th = z * cfg['tile_size']

    ntx = np.ceil(float(w - ov) / (tw - ov))
    nty = np.ceil(float(h - ov) / (th - ov))
    nt = ntx * nty

    print 'tiles size: (%d, %d)' % (tw, th)
    print 'total number of tiles: %d (%d x %d)' % (nt, ntx, nty)
    nb_pairs = len(cfg['images']) - 1
    print 'total number of pairs: %d' % nb_pairs

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
            tile_info['number_of_pairs'] = nb_pairs
            tile_info['images'] = cfg['images']
            tiles_full_info.append(tile_info)

    if len(tiles_full_info) == 1:
        tiles_full_info[0]['position_type'] = 'Single'

    cutting_info=open(os.path.join(cfg['out_dir'],'list_of_tiles.txt'),'w')
    for tile_info in tiles_full_info:
        cutting_info.write( '%s\n' % (tile_info['directory']))
    cutting_info.close()
    
    return tiles_full_info
