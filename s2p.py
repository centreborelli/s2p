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
# F
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import multiprocessing
import sys
import json
import numpy as np
import os.path
from os import mkdir
from shutil import rmtree
import copy
import glob
import time

from python import tee
from python import srtm
from python import common
from python import masking
from python import rpc_model
from python import rpc_utils
from python import geographiclib
from python import pointing_accuracy
from python import visualisation
from python import rectification
from python import block_matching
from python import triangulation
from python import tile_composer
from python import fusion
from python.config import cfg



def show_progress(a):
    show_progress.counter += 1
    if show_progress.counter > 1:
        print "Processed %d tiles" % show_progress.counter
    else:
        print "Processed 1 tile"


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


# ----------------------------------------------------------------------------------------------------
# ---------------------------------------  initialize ---------------------------------------
# ----------------------------------------------------------------------------------------------------

def initialize(config_file):
    """
    Prepares the entire process : 
    1) Select a ROI, check it, make sure its coordinates are multiples of the zoom factor, and so forth.
    2) Compute optimal size for tiles, get the number of pairs
    3) Build the tileComposerInfo : a simple list that gives the size of the tiles (after having removed the overlapping areas), regarding their positions inside the ROI.
       It will be useful for the contruction of the vrt file that merge all the local height maps and other data.
    4) Build tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]. USED EVERYWHERE IN THIS CODE.
       * col/row : position of the tile (upper left corner)
       * tw/th : size of the tile
       * ov : size of the overlapping
       * i/j : relative position of the tile
       * pos : position inside the ROI : UL for a tile place at th Upper Left corner, M for the ones placed in the middle, and so forth.
       * images : a dictionary directly given by the json config file, that store the information about all the involved images, their rpc, and so forth.
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
    


    # create tmp dir and output directory for the experiment, and store a json
    # dump of the config.cfg dictionary there
    if not os.path.exists(cfg['temporary_dir']):
        os.makedirs(cfg['temporary_dir'])
    if not os.path.exists(os.path.join(cfg['temporary_dir'], 'meta')):
        os.makedirs(os.path.join(cfg['temporary_dir'], 'meta'))
    if not os.path.exists(cfg['out_dir']):
        os.makedirs(cfg['out_dir'])
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


    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    z = cfg['subsampling_factor']
    x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)
    

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

# ----------------------------------------------------------------------------------------------------------
# ---------------------------------------  initialize (end) ---------------------------------------
# ----------------------------------------------------------------------------------------------------------
    

# ----------------------------------------------------------------------------------------------------
# ---------------------------------------  preprocess_tiles ------------------------------------
# ----------------------------------------------------------------------------------------------------


def preprocess_tile(tile_dir, tilesFullInfo):
    """
    Compute pointing corrections, crop ref image into tile, and get min/max intensities values 
    """
    
    pointing_correction(tile_dir, tilesFullInfo)
    common.getMinMaxFromExtract(tile_dir,tilesFullInfo)


def pointing_correction(tile_dir, tilesFullInfo):
    """                        
    Computes pointing corrections
    
    Args: 
         - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]         
    """     
                                                          
    
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=tilesFullInfo[tile_dir]
            

    for i in range(0,NbPairs):
        pair_id = i+1
        paired_tile_dir=tile_dir + 'pair_%d' % pair_id

        img1,rpc1,img2,rpc2 = images[0]['img'],images[0]['rpc'],images[pair_id]['img'],images[pair_id]['rpc']
        
        
        # create a directory for the experiment
        if not os.path.exists(paired_tile_dir):
            os.makedirs(paired_tile_dir)
        
        # redirect stdout and stderr to log file
        if not cfg['debug']:
            fout = open('%s/stdout.log' % paired_tile_dir, 'w', 0)  # '0' for no buffering
            sys.stdout = fout
            sys.stderr = fout
    
        # debug print
        print 'tile %d %d running on process %s' % (x, y,
                                                    multiprocessing.current_process())
    
        # output files
        cwid_msk = '%s/cloud_water_image_domain_mask.png' % (paired_tile_dir)
        pointing = '%s/pointing.txt' % paired_tile_dir
        center = '%s/center_keypts_sec.txt' % paired_tile_dir
        sift_matches = '%s/sift_matches.txt' % paired_tile_dir
        sift_matches_plot = '%s/sift_matches_plot.png' % paired_tile_dir


        # check if the tile is already done
        if os.path.isfile('%s/pointing.txt' % paired_tile_dir) and cfg['skip_existing']:
            print "pointing correction on tile %d %d (pair %d) already done, skip" % (col,row,pair_id)
        else:
        
            # check if the ROI is completely masked (water, or outside the image domain)
            H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
            if masking.cloud_water_image_domain(cwid_msk, w, h, H, rpc1, roi_msk,cld_msk):
                print "Tile masked by water or outside definition domain, skip"
                open("%s/this_tile_is_masked.txt" % paired_tile_dir, 'a').close()
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                if not cfg['debug']:
                    fout.close()
            else:
            
                # correct pointing error
                # A is the correction matrix and m is the list of sift matches
                A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, x,
                                                                y, w, h)
                if A is not None:
                    np.savetxt(pointing, A)
                if m is not None:
                    np.savetxt(sift_matches, m)
                    np.savetxt(center, np.mean(m[:, 2:4], 0))
                    visualisation.plot_matches_pleiades(img1, img2, rpc1, rpc2, m, x, y,
                                                            w, h, sift_matches_plot)
       
        # close logs
        common.garbage_cleanup()
        if not cfg['debug']:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            fout.close()

    return
        
# ----------------------------------------------------------------------------------------------------
# ------------------------------------------  preprocess_tiles (end) ---------------------------------
# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# --------------------------------------------  global_values  ---------------------------------------
# ----------------------------------------------------------------------------------------------------


def global_values(tilesFullInfo):
    """
    Compute the global pointing correction + Compute the min and max intensities from the tiles that will be processed.
    """
    global_pointing_correction(tilesFullInfo)
    global_minmax_intensities(tilesFullInfo)
    


def global_pointing_correction(tilesFullInfo):
    """
    Compute the global pointing correction. For that purpose, global_pointing_correction needs all the tile related to one pair (given by pairedTilesPerPairId),
    and the total number of pairs as well.
    """

    # Build pairedTilesPerPairId : a dictionary that provides the list of all the tile for a given pair_id (key)
    pairedTilesPerPairId = {}
    for tile_dir in tilesFullInfo:
        col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=tilesFullInfo[tile_dir]

        for i in range(0,NbPairs):
            pair_id = i+1
            pair_dir = tile_dir + 'pair_%d' % pair_id
            if pairedTilesPerPairId.has_key(pair_dir):
                pairedTilesPerPairId[pair_id].append(pair_dir)
            else:
                pairedTilesPerPairId[pair_id]=[]
     # Build pairedTilesPerPairId (end)

    for pair_id in pairedTilesPerPairId:
        tiles = pairedTilesPerPairId[pair_id]
        A_globalMat = pointing_accuracy.global_from_local(tiles)
        np.savetxt('%s/global_pointing_pair_%d.txt' % (cfg['out_dir'],pair_id), A_globalMat)
		

def global_minmax_intensities(tilesFullInfo):
    """
    Compute the min and max intensities from the tiles that will be processed.
    This will allow to re-code colors by using 8-bits instead of 12-bits or more, and to better vizualise the ply files.

    Args:
        - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]

    """

    minlist=[]
    maxlist=[]
    for tile_dir in tilesFullInfo:
        minmax = np.loadtxt(tile_dir + '/local_minmax.txt')
        minlist.append(minmax[0])
        maxlist.append(minmax[1])
		
    global_minmax=[min(minlist),max(maxlist)]
	
    np.savetxt(cfg['out_dir']+'/global_minmax.txt',global_minmax)
        
# ----------------------------------------------------------------------------------------------------
# ---------------------------------------------  global_values --------------------------------------
# ----------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------
# ---------------------------------------------- process_tiles ---------------------------------------
# ----------------------------------------------------------------------------------------------------

def color_crop_ref(tile_dir,tilesFullInfo,clr=None):

    """
    Colorization of a crop_ref (for a given tile)
    
    Args:
        - tile_dir : a key for the dictionary tilesFullInfo; refers to a particular tile
        - tilesFullInfo : a dictionary that provides all you need to process a tile -> col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk
        
        - clr (optional) : if crop_ref is a pan image, will perform the pansharpening with the color image clr
            
            If clr is None then:
                case 1 : tile is an RGBI image, so removes I channel, and perform rescaling of the remaining channels
                case 2 : tile is already an RGB image, so just perform rescaling

                Note that if rescaling is already performed, then the file applied_minmax.txt exists :
                - if applied_minmax.txt exists and cfg['skip_existing'] is True, rescaling won't be performed again
                - if applied_minmax.txt exists and is different from global_minmax, rescaling will be compulsorily performed (can occur if a new tile is added)
    """

    # Get info
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    z = cfg['subsampling_factor']

    # Paths
    #crop_ref = tile_dir + '/roi_ref.tif'
    crop_ref = tile_dir + '/roi_ref_crop.tif'
    global_minmax = cfg['out_dir']+'/global_minmax.txt'
    applied_minmax = tile_dir + '/applied_minmax.txt'
    
    global_minmax_arr = np.loadtxt(global_minmax)
            
    if cfg['color_ply']:
        
        doProcess=False
        if not os.path.exists(applied_minmax):
            doProcess=True
            applied_minmax_arr = global_minmax_arr
        else:     
            applied_minmax_arr = np.loadtxt(applied_minmax)
        
            if (applied_minmax_arr[0]!=global_minmax_arr[0]) or (applied_minmax_arr[1]!=global_minmax_arr[1]):
                doProcess=True
                applied_minmax_arr = global_minmax_arr
        

        if not doProcess and cfg['skip_existing']:
                print 'Rescaling of tile %s already done, skip' % tile_dir
        else:

            crop_color = tile_dir + '/roi_color_ref.tif'
            if clr is not None:
                triangulation.colorize(crop_ref, clr, col,row, z, crop_color, applied_minmax_arr[0],applied_minmax_arr[1])
            else: # use of image_rescaleintensities
            
                np.savetxt(applied_minmax,applied_minmax_arr)
    
                if common.image_pix_dim_tiffinfo(crop_ref) == 4:
                    print 'the image is pansharpened fusioned'
                    tmp = common.rgbi_to_rgb(crop_ref, out=None, tilewise=True)
                    #common.image_qauto(tmp, crop_color, tilewise=False)
                    common.image_rescaleintensities(tmp,crop_color,applied_minmax_arr[0],applied_minmax_arr[1])
                else:
                    print 'no color data'
                    #common.image_qauto(crop_ref, crop_color, tilewise=False)
                    common.image_rescaleintensities(crop_ref,crop_color,applied_minmax_arr[0],applied_minmax_arr[1])
        
        

def generate_cloud(tile_dir,tilesFullInfo, do_offset=False):
    """
    Args:
        - tile_dir : a key for the dictionary tilesFullInfo; refers to a particular tile
        - tilesFullInfo : a dictionary that provides all you need to process a tile -> col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk
        - do_offset (optional, default: False): boolean flag to decide wether the
            x, y coordinates of points in the ply file will be translated or
            not (translated to be close to 0, to avoid precision loss due to
            huge numbers)
    """
    print "\nComputing point cloud..."
    
    # Get info
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    img1,rpc1 = images[0]['img'],images[0]['rpc']
    
    #height_map = tile_dir + '/local_merged_height_map.tif'
    height_map = tile_dir +'/local_merged_height_map_crop.tif'
    crop_color = tile_dir + '/roi_color_ref.tif'
    if not os.path.exists(crop_color):
        crop_color=''

    z = cfg['subsampling_factor']

    #Compute the homography transforming the coordinates system of the
    #original full size image into the coordinates system 
    #of the crop we are dealing with
    A = common.matrix_translation(-col*z, -row*z) # col and row have been divided by z inside 'finalize_tile' for convinience; col*z and row*z allow to get the initial values back.
    z = cfg['subsampling_factor']
    f = 1.0/z
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = tile_dir + '/trans.txt'
    np.savetxt(trans, A)

    # compute coordinates (offsets) of the point we want to use as origin in the local coordinate system of the computed cloud
    if do_offset:
        r = rpc_model.RPCModel(rpc1)
        lat = r.latOff
        lon = r.lonOff
        off_x, off_y = geographiclib.geodetic_to_utm(lat, lon)[0:2]
    else:
        off_x, off_y = 0, 0



    # output 
    cloud = tile_dir + '/cloud.ply'
    
    triangulation.compute_point_cloud(cloud, height_map, rpc1, trans, crop_color,
                                      off_x, off_y)
                                                                            
    common.garbage_cleanup()
    
    
def merge_height_maps(height_maps,tile_dir,thresh,conservative,k=1,garbage=[]):
    """
    Merge a list of height maps recursively, computed for one tile from N image pairs.
  
    Args :
         - height_maps : list of height map directories
         - tile_dir : directory of the tile from which to get a merged height map
         - thresh : threshold used for the fusion algorithm, in meters.
         - conservative (optional, default is False): if True, keep only the
            pixels where the two height map agree (fusion algorithm)
         - k : used to identify the current call of merge_height_maps (default = 1, first call)
         - garbage : a list used to remove temp data (default = [], first call)

    """

    #output file
    local_merged_height_map = tile_dir +'/local_merged_height_map.tif'

    if os.path.isfile(local_merged_height_map) and cfg['skip_existing']:
        print 'final height map %s already done, skip' % local_merged_height_map
    else:
        list_height_maps=[]
        for i in range(len(height_maps)-1):
            height_map = tile_dir +'/height_map_'+str(i)+'_'+str(i+1)+'_'+str(k)+'.tif'
            fusion.merge(height_maps[i], height_maps[i+1], thresh, height_map,
                     conservative)
            list_height_maps.append(height_map)
            garbage.append(height_map)
        
        if len(list_height_maps) > 1:
            merge_height_maps(list_height_maps,tile_dir,thresh,conservative,k+1,garbage)
        else:
            common.run('cp %s %s' % (list_height_maps[0],local_merged_height_map))
            for imtemp in garbage:
                common.run('rm -f %s' % imtemp )    


def finalize_tile(tile_dir, height_maps, tilesFullInfo):
    """
    Finalize the processing of a tile : merge the height maps from the N pairs, remove overlapping areas, 
    get the colors from a XS image and use it to color and generate a ply file (colorization is not mandatory)

    Args:
        - tile_dir : directory of the tile to be processed
        - height_maps : list of the height maps generated from N pairs
        - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    """

    # Get all info
    fullInfo=tilesFullInfo[tile_dir]
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=fullInfo
    img1,rpc1 = images[0]['img'],images[0]['rpc']


    # merge the n height maps
    local_merged_height_map = tile_dir +'/local_merged_height_map.tif'
    if len(height_maps)>1:
        merge_height_maps(height_maps,tile_dir,cfg['fusion_thresh'],cfg['fusion_conservative'],1,[])
    else:    
        common.run('cp %s %s' % (height_maps[0],local_merged_height_map))

    # Remove overlapping areas
    #By tile
    local_merged_height_map = tile_dir +'/local_merged_height_map.tif'
    local_merged_height_map_crop = tile_dir +'/local_merged_height_map_crop.tif'
    crop_ref = tile_dir + '/roi_ref.tif'
    crop_ref_crop = tile_dir + '/roi_ref_crop.tif'


    dicoPos={}
    dicoPos['M']  = [ov/2,ov/2,-ov,-ov]
    dicoPos['L']  = [0,ov/2,-ov/2,-ov]
    dicoPos['R']  = [ov/2,ov/2,-ov/2,-ov]
    dicoPos['U']  = [ov/2,0,-ov,-ov/2]
    dicoPos['B']  = [ov/2,ov/2,-ov,-ov/2]
    dicoPos['UL'] = [0,0,-ov/2,-ov/2]
    dicoPos['UR'] = [ov/2,0,-ov/2,-ov/2]
    dicoPos['BR'] = [ov/2,ov/2,-ov/2,-ov/2]
    dicoPos['BL'] = [0,ov/2,-ov/2,-ov/2]
    dicoPos['Single'] = [0,0,0,0]  

    z = cfg['subsampling_factor']
    #tilesFullInfo[tile_dir] = dicoPos[pos][0]
    newcol,newrow,difftw,diffth = np.array(dicoPos[pos])/z
    info=tilesFullInfo[tile_dir]
    info[0] = info[0]/z + newcol
    info[1] = info[1]/z + newrow
    info[2] = info[2]/z + difftw
    info[3] = info[3]/z + diffth
    tilesFullInfo[tile_dir]=info   
   
    # z=1 beacause local_merged_height_map, crop_ref (and so forth) have already been zoomed. So don't zoom again to crop these images.
    common.cropImage(local_merged_height_map,local_merged_height_map_crop,newcol,newrow,info[2],info[3],1)
    common.cropImage(crop_ref,crop_ref_crop,newcol,newrow,info[2],info[3],1)   

        
    #By pair    
    for i in range(0,NbPairs):
        pair_id = i+1
            
        single_height_map = tile_dir + '/pair_%d/height_map.tif' % pair_id
        single_height_map_crop =tile_dir + '/pair_%d/height_map_crop.tif' % pair_id
          
        single_rpc_err = tile_dir + '/pair_%d/rpc_err.tif' % pair_id
        single_rpc_err_crop =tile_dir + '/pair_%d/rpc_err_crop.tif' % pair_id
            
        common.cropImage(single_height_map,single_height_map_crop,newcol,newrow,info[2],info[3],1)
        common.cropImage(single_rpc_err,single_rpc_err_crop,newcol,newrow,info[2],info[3],1)

    
    
    #if not cfg['debug']:
     #    common.run('rm -f %s' % (local_merged_height_map)) 
         #don't remove crop_ref !!
                   
    
    # Colors 
    color_crop_ref(tile_dir,tilesFullInfo,cfg['images'][0]['clr'])
    
    # Generate cloud 
    generate_cloud(tile_dir,tilesFullInfo,cfg['offset_ply'])
    
    

def process_pair(tile_dir,pair_id,tilesFullInfo):
    """
    For a given tile, process pair #pair_id : rectification, disparity map, triangulation, 

    Args:
         - tile_dir : directory of the tile 
         - pair_id : id of the pair to be processed
         - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    """

    #Get all info
    fullInfo=tilesFullInfo[tile_dir]
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=fullInfo
    img1,rpc1 = images[0]['img'],images[0]['rpc']

    paired_tile_dir=tile_dir + 'pair_%d' % pair_id

    img2,rpc2 = images[pair_id]['img'],images[pair_id]['rpc']
    A_global = '%s/global_pointing_pair_%d.txt' % (cfg['out_dir'],pair_id)

    if os.path.isfile('%s/this_tile_is_masked.txt' % tile_dir):
        print "tile %s already masked, skip" % paired_tile_dir
        return None

    else:
        print 'process tile %d %d...' % (col,row)
       
        #height_maps.append('%s/height_map.tif' % paired_tile_dir)

        # Rectification
        if os.path.isfile('%s/rectified_ref.tif' % paired_tile_dir) and os.path.isfile('%s/rectified_sec.tif' % paired_tile_dir) and cfg['skip_existing']:
            print 'rectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print 'rectification : process tile %d %d (pair %d)...' % (col,row,pair_id)
            rectify(paired_tile_dir, A_global, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk)
        
        # Disparity
        if os.path.isfile('%s/rectified_disp.tif' % paired_tile_dir) and cfg['skip_existing']:
            print 'disparity on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print 'disparity : process tile %d %d (pair %d)...' % (col,row,pair_id)
            disparity(paired_tile_dir, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk)
            
        # Triangulate
        if os.path.isfile('%s/height_map.tif' % paired_tile_dir) and cfg['skip_existing']:
            print 'Triangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print 'Triangulation : process tile %d %d (pair %d)...' % (col,row,pair_id)    
            triangulate(paired_tile_dir, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk, np.loadtxt(A_global))

        return '%s/height_map.tif' % paired_tile_dir
        
        

def process_tile(tile_dir,tilesFullInfo):
    """
    Process a tile : compute the height maps from the N pairs, and finalize the processing, ie.:
    1) Produce a merged height map, without overlapping areas, 
    2) And produce a ply file.

    Args:
        - tile_dir : directory of the tile to be processed
        - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    """
    
    # Process each pair : get a height map
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=tilesFullInfo[tile_dir]
    height_maps = []
    for i in range(0,NbPairs):
        pair_id = i+1
        height_map = process_pair(tile_dir,pair_id,tilesFullInfo)
        if height_map:
            height_maps.append(height_map)


    ## Finalization ##
    finalize_tile(tile_dir, height_maps, tilesFullInfo)


    
# ----------------------------------------------------------------------------------------------------
# ------------------------------------------- process_tiles ------------------------------------------
# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# ----------------------------------- rectify disparity triangulate ----------------------------------
# ----------------------------------------------------------------------------------------------------

def rectify(out_dir, A_global, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None):
    """
    Computes rectifications, without tiling

    Args:
        out_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.

    Returns:
        nothing
    """
    
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir


    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())
                                                
    
    A,m=None,None

    if os.path.isfile('%s/pointing.txt' % out_dir): 
        A=np.loadtxt('%s/pointing.txt' % out_dir)
    else:
        A=A_global
    if os.path.isfile('%s/sift_matches.txt' % out_dir):
        m=np.loadtxt('%s/sift_matches.txt' % out_dir)

    # rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m)
    
    z = cfg['subsampling_factor']                           
    np.savetxt(subsampling, np.array([z]))                                                        
    np.savetxt(H_ref, H1)
    np.savetxt(H_sec, H2)
    np.savetxt(disp_min_max, np.array([disp_min, disp_max]))



    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return


def disparity(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None):
    """
    Computes a disparity map from a Pair of Pleiades images, without tiling

    Args:
        out_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.

    Returns:
        nothing
    """
    
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir


    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())
                                                

    # disparity (block-matching)
    
    disp_min,disp_max=np.loadtxt(disp_min_max)
    
    if cfg['disp_min'] is not None:
        disp_min = cfg['disp_min']
    if cfg['disp_max'] is not None:
        disp_max = cfg['disp_max']
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # intersect mask with the cloud_water_image_domain mask (recomputed here to
    # get to be sampled on the epipolar grid)
    ww, hh = common.image_size(rect1)
    H1 = np.loadtxt(H_ref) 
    masking.cloud_water_image_domain(cwid_msk, ww, hh, H1, rpc1, roi_msk,
                                     cld_msk)
    try:
        masking.intersection(mask, mask, cwid_msk)
        masking.erosion(mask, mask, cfg['msk_erosion'])
    except OSError:
        print "file %s not produced" % mask


    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return




def triangulate(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None, A=None):
    """
    Computes triangulations, without tiling

    Args:
        out_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.
        A (optional, default None): pointing correction matrix. 

    Returns:
        nothing
    """
    
    # output files
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    rpc_err = '%s/rpc_err.tif' % out_dir
    height_map = '%s/height_map.tif' % out_dir

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())

    
    # triangulation
    z = cfg['subsampling_factor']
    triangulation.compute_dem(height_map, x, y, w, h, z,
                                              rpc1, rpc2, H_ref, H_sec, disp, mask,
                                              rpc_err, A)
                                                            
    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return

# ----------------------------------------------------------------------------------------------------
# ------------------------------- rectify disparity triangulate (end) --------------------------------
# ----------------------------------------------------------------------------------------------------					 
        

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------- global finalization --------------------------------------
# ----------------------------------------------------------------------------------------------------   

def write_vrt_files(tilesFullInfo):
    """
    Merge pieces of data into single VRT files : height map comprising the N pairs, height map for each signle pair, and err_rpc 

    Args:
         - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    """
    #VRT file : height map (N pairs)    

    #-------------- tileComposerInfo -------------- 
    # tileComposerInfo : a simple list that gives the size of the tiles (after having removed the overlapping areas), 
    # regarding their positions inside the ROI --> ULw , ULh, Lw , Lh, Uw , Uh, Mw , Mh
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=tilesFullInfo.values()[0]
    ULw , ULh = tw-ov/2, th-ov/2 #Size of Corner tiles (UL UR BL BR)
    Lw , Lh   = tw-ov/2,th-ov    #Size of Left/Right tile (L R)
    Uw , Uh   = tw-ov,th-ov/2    #Size of Upper/Bottom tile (U B)
    Mw , Mh   = tw-ov,th-ov      #Size of Middle tile

    rangey = np.arange(y, y + h - ov, th - ov)
    rangex = np.arange(x, x + w - ov, tw - ov)
    rowmin,rowmax = rangey[0],rangey[-1]
    colmin,colmax = rangex[0],rangex[-1]
    
    #Tile composer info (bis)
    imax=len(rangey)-1
    jmax=len(rangex)-1
    fh = 2*ULh+(imax-1)*Mh #height of the final height map after removing margins
    fw = 2*ULw+(jmax-1)*Mh #width of the final height map after removing margins
    #-------------- tileComposerInfo (end) --------------
    
    
    tileSizesAndPositions={}
    for tile_dir in tilesFullInfo:

        col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=tilesFullInfo[tile_dir]
        
        dicoPos={}
        dicoPos['M']  = [ULw+(j-1)*Mw, ULh+(i-1)*Mh , Mw, Mh]
        dicoPos['UL'] = [0,0, ULw, ULh]
        dicoPos['UR'] = [fw-ULw,0, ULw, ULh]
        dicoPos['BR'] = [fw-ULw,fh-ULh, ULw, ULh]
        dicoPos['BL'] = [0,fh-ULh, ULw, ULh]
        dicoPos['L']  = [0,ULh+(i-1)*Lh, Lw, Lh]
        dicoPos['R']  = [fw-ULw,ULh+(i-1)*Lh, Lw, Lh]
        dicoPos['U']  = [ULw+(j-1)*Uw, 0, Uw, Uh]
        dicoPos['B']  = [ULw+(j-1)*Uw, fh-ULh, Uw, Uh]
        dicoPos['Single'] = [0,0, tw, th]
        
        tile_reldir = 'tile_%d_%d_row_%d/col_%d/' % (tw, th, row, col)

        tileSizesAndPositions[tile_reldir] = dicoPos[pos]  
   
    
    z = cfg['subsampling_factor']
    tile_composer.mosaic_gdal2(cfg['out_dir'] + '/heightMap_N_pairs.vrt', tileSizesAndPositions, 'local_merged_height_map_crop.tif', fw, fh,z)
    

    # VRT file : height map (for each single pair)
    # VRT file : rpc_err (for each single pair)
    for i in range(0,NbPairs):
        
        pair_id = i+1
        
        pairSizesAndPositions={}
        for tile_reldir in tileSizesAndPositions:
            pair_reldir = tile_reldir + '/pair_%d' % (pair_id)
            pairSizesAndPositions[pair_reldir] = tileSizesAndPositions[tile_reldir]  
        
        tile_composer.mosaic_gdal2(cfg['out_dir'] + '/heightMap_pair_%d.vrt' % (pair_id), pairSizesAndPositions, 'height_map_crop.tif', fw, fh,z)
        tile_composer.mosaic_gdal2(cfg['out_dir'] + '/rpc_err_pair_%d.vrt' % (pair_id), pairSizesAndPositions, 'rpc_err_crop.tif', fw, fh,z)    



def write_dsm(tilesFullInfo):
    """
    Write the DSM, from the ply files given by each tile.

    Args :
         - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    """
    clouds_dir = cfg['out_dir']+'/clouds'
    if (os.path.exists(clouds_dir)):
        rmtree(clouds_dir)
    mkdir(clouds_dir)    
    
    for tile_dir in tilesFullInfo:
        
        col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=tilesFullInfo[tile_dir]
        cloud = tile_dir + '/cloud.ply'
        cloud_link_name = clouds_dir + '/cloud_%d_%d_row_%d_col_%d.ply' % (tw, th, row, col)
        common.run('ln %s %s' % (cloud,cloud_link_name) )
        
    out_dsm = '%s/dsm.tif' % (cfg['out_dir'])
    common.run("ls %s | plyflatten %f %s" % (clouds_dir+'/cloud*',cfg['dsm_resolution'], out_dsm))



def global_finalization(tilesFullInfo):
    """
    Merge pieces of data into single VRT files : height map comprising the N pairs, height map for each signle pair, and err_rpc.
    Write the DSM, from the ply files given by each tile.
    Crop corresponding areas in the secondary images.
    Copy the RPC'. 

    Args:
         - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]

    """

    write_vrt_files(tilesFullInfo)
    write_dsm(tilesFullInfo)

    if not cfg['full_img']:
        common.crop_corresponding_areas(cfg['out_dir'], cfg['images'], cfg['roi'])

    for i in range(len(cfg['images'])):
        from shutil import copy2
        copy2(cfg['images'][i]['rpc'], cfg['out_dir'])

# ----------------------------------------------------------------------------------------------------
# --------------------------------------- global finalization (end) ----------------------------------
# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# ------------------------------------------- map_processing -----------------------------------------
# ----------------------------------------------------------------------------------------------------
def map_processing(config_file):
    """
    Init +  reprocessing + global_values + processing + global finalization

    Args: 
         - json config file
    
    """
    
    try:
    
        tilesFullInfo = initialize(config_file)
    
        if cfg['debug']: #monoprocessing
    
            print 'preprocess_tile...'
            for tile_dir in tilesFullInfo:
                preprocess_tile(tile_dir, tilesFullInfo)
            print 'global values...'
            
            global_values(tilesFullInfo)
            
            print 'process_tile...'
            for tile_dir in tilesFullInfo:
                process_tile(tile_dir, tilesFullInfo)
                
            print 'global finalization...'     
            global_finalization(tilesFullInfo)
            
        else: # multiprocessing
        
            # create pool with less workers than available cores
            nb_workers = multiprocessing.cpu_count()
            if cfg['max_nb_threads']:
                nb_workers = min(nb_workers, cfg['max_nb_threads'])
            
            print 'preprocess_tile...'
            results = []
            show_progress.counter = 0
            pool = multiprocessing.Pool(nb_workers)
            for tile_dir in tilesFullInfo:
    
                p = pool.apply_async(preprocess_tile,
                                             args=(tile_dir, tilesFullInfo), callback=show_progress)
                results.append(p)
                
            for r in results:
                try:
                    r.get(3600)  # wait at most one hour per tile
                except multiprocessing.TimeoutError:
                    print "Timeout while computing tile "+str(r)       
                    
            print 'global values...'        
            global_values(tilesFullInfo)        
                    
            print 'process_tile...'
            results = []
            show_progress.counter = 0   
            pool = multiprocessing.Pool(nb_workers)     
            for tile_dir in tilesFullInfo:
                
                p = pool.apply_async(process_tile,
                                             args=(tile_dir, tilesFullInfo), callback=show_progress)
                results.append(p)
                
            for r in results:
                try:
                    r.get(3600)  # wait at most one hour per tile
                except multiprocessing.TimeoutError:
                    print "Timeout while computing tile "+str(r)
            
            print 'global finalization...'        
            global_finalization(tilesFullInfo)
    
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
# ----------------------------------------------------------------------------------------------------
# ------------------------------------- map_processing (end) -----------------------------------------
# ----------------------------------------------------------------------------------------------------


def main(config_file):
    """
    Launches s2p with the parameters given by a json file.

    Args:
        config_file: path to the config json file
    """


    # measure total runtime
    t0 = time.time()


    ###########################
    map_processing(config_file)
    ###########################                
        

    # runtime
    t = int(time.time() - t0)
    h = t/3600
    m = (t/60) % 60
    s = t % 60
    print "Total runtime: %dh:%dm:%ds" % (h, m, s)
    common.garbage_cleanup()


if __name__ == '__main__':

    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print """
        Incorrect syntax, use:
          > %s config.json

          Launches the s2p pipeline. All the parameters, paths to input and
          output files, are defined in the json configuration file.
        """ % sys.argv[0]
        sys.exit(1)
