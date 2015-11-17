# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import numpy as np
from config import cfg
import sys
import os
import json


from python import common
from python import srtm
from python import tee



def check_parameters(usr_cfg):
    """
    Checks that the provided dictionary defines all the mandatory
    arguments, and warns about unknown optional arguments.

    Args:
        usr_cfg: python dict read from the json input file

    Returns:
        nothing
    """

    # verify that i/o files are defined
    if 'out_dir' not in usr_cfg:
        print "missing output dir: abort"
        sys.exit(1)
    if 'images' not in usr_cfg or len(usr_cfg['images']) < 2:
        print "missing input data paths: abort"
        sys.exit(1)
    for i in range(0,len(usr_cfg['images'])):
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
        
        
def prepare_config(config_file):
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

    prepare_config(config_file)

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
		


def init_tilesFullInfo(config_file):
    """
    Prepares the entire process : 
    1) Makes sure coordinates of the ROI are multiples of the zoom factor
    2) Computes optimal size for tiles, get the number of pairs
    3) Builds tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]. USED EVERYWHERE IN THIS CODE.
       * col/row : position of the tile (upper left corner)
       * tw/th : size of the tile
       * ov : size of the overlapping
       * i/j : relative position of the tile
       * pos : position inside the ROI : UL for a tile place at th Upper Left corner, M for the ones placed in the middle, and so forth.
       * images : a dictionary directly given by the json config file, that store the information about all the involved images, their rpc, and so forth.
       
    Args:
         - config_file : a json configuratio file
         - write : if False, initialize will just create the output dir and nothing else. 
         
    Returns:
         - tilesFullInfo
    """          

    prepare_config(config_file)

    #Get ROI
    x = cfg['roi']['x']    
    y = cfg['roi']['y']
    w = cfg['roi']['w']
    h = cfg['roi']['h']    
    z = cfg['subsampling_factor']


    # Automatically compute optimal size for tiles
    #tw, th : dimensions of the tiles
    #ov : width of overlapping bands between tiles
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
    NbPairs = len(cfg['images'])-1
    print 'total number of pairs: %d ' % NbPairs


    # Build tiles dicos
    tilesFullInfo={}
    
    rangey = np.arange(y, y + h - ov, th - ov)
    rangex = np.arange(x, x + w - ov, tw - ov)
    rowmin,rowmax = rangey[0],rangey[-1]
    colmin,colmax = rangex[0],rangex[-1]
    

    for pair_id in range(1,len(cfg['images'])) :
        for i, row in enumerate(rangey):
            for j, col in enumerate(rangex):
                # ensure that the coordinates of the ROI are multiples of the zoom factor
                col, row, tw, th = common.round_roi_to_nearest_multiple(z, col, row, tw, th)
                tile_dir = '%s/tile_%d_%d_row_%d/col_%d/' % (cfg['out_dir'], tw, th, row, col)
                paired_tile_dir = '%s/tile_%d_%d_row_%d/col_%d/pair_%d' % (cfg['out_dir'], tw, th, row, col, pair_id)
                
                if row==rowmin and col==colmin:
                    pos='UL'
                elif row==rowmin and col==colmax:
                    pos='UR'
                elif row==rowmax and col==colmax:
                    pos='BR'
                elif row==rowmax and col==colmin:
                    pos='BL'
                elif row==rowmin and col>colmin:
                    pos='U'
                elif col==colmin and row>rowmin:
                    pos='L'
                elif row==rowmax and col>colmin:
                    pos='B'
                elif col==colmax and row>rowmin:
                    pos='R'
                else:
                    pos='M'
                
                tilesFullInfo[tile_dir]=[col,row,tw,th,ov,i,j,pos,x,y,w,h,cfg['images'],NbPairs,cfg['images'][0]['cld'],cfg['images'][0]['roi']]


    if len(tilesFullInfo) == 1 : #Only one tile : put position 'Single'
        tile_dir,info=tilesFullInfo.items()[0]
        col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=info
        tilesFullInfo[tile_dir]=[col,row,tw,th,ov,i,j,'Single',x,y,w,h,images,NbPairs,cld_msk,roi_msk]

    return tilesFullInfo
