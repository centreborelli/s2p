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
import numpy as np
import os.path
import time

from python import common
from python.config import cfg


from python import initialization
from python import preprocess
from python import globalvalues
from python import process
from python import globalfinalization


def show_progress(a):
    show_progress.counter += 1
    if show_progress.counter > 1:
        print "Processed %d tiles" % show_progress.counter
    else:
        print "Processed 1 tile"





# ----------------------------------------------------------------------------------------------------
# ------------------------------------------  initialize ---------------------------------------------
# ----------------------------------------------------------------------------------------------------

def init_dirs(config_file):
	
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
		

def initialize(config_file):
    """
    1) Loads configuration file
    2) Checks parameters
    3) Selects a ROI, checks the zoom factor; make sure coordinates of the ROI are multiples of the zoom factor
    4) Creates different directories : output, temp...
    5) Builds tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]. USED EVERYWHERE IN THIS CODE.
       * col/row : position of the tile (upper left corner)
       * tw/th : size of the tile
       * ov : size of the overlapping
       * i/j : relative position of the tile
       * pos : position inside the ROI : UL for a tile place at th Upper Left corner, M for the ones placed in the middle, and so forth.
       * images : a dictionary directly given by the json config file, that store the information about all the involved images, their rpc, and so forth.
    
    Args :
         - config_file : a json configuratio file
    Returns :
         - tilesFullInfo 
    """


    initialization.init_dirs_srtm_roi(config_file)
    tilesFullInfo = initialization.init_tilesFullInfo(config_file)
   
    return tilesFullInfo

# -----------------------------------------------------------------------------------------------------
# ---------------------------------------  initialize (end) -------------------------------------------
# -----------------------------------------------------------------------------------------------------
    

# ----------------------------------------------------------------------------------------------------
# ---------------------------------------  preprocess_tiles ------------------------------------------
# ----------------------------------------------------------------------------------------------------

def preprocess_tile(tile_dir, tilesFullInfo):
    """
    1) Computes pointing corrections, 
    2) Crops ref image into tiles and get min/max intensities values for each one 
    """
    
    preprocess.pointing_correction(tile_dir, tilesFullInfo)
    preprocess.getMinMaxFromExtract(tile_dir,tilesFullInfo)
        
# ----------------------------------------------------------------------------------------------------
# ------------------------------------------  preprocess_tiles (end) ---------------------------------
# ----------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------
# --------------------------------------------  global_values  ---------------------------------------
# ----------------------------------------------------------------------------------------------------

def global_values(tilesFullInfo):
    """
    1) Computes the global pointing correction 
    2) Computes the min and max intensities from the tiles that will be processed (to later rescale colors to 8-bits ones, useful for utils such as cloudcompare)
    """
    globalvalues.global_pointing_correction(tilesFullInfo)
    globalvalues.global_minmax_intensities(tilesFullInfo)
    
# ----------------------------------------------------------------------------------------------------
# ---------------------------------------------  global_values --------------------------------------
# ----------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------
# ---------------------------------------------- process_tiles ---------------------------------------
# ----------------------------------------------------------------------------------------------------

def process_pair(tile_dir,pair_id,tilesFullInfo):
    """
    For a given tile, processes pair #pair_id : rectification, disparity map, triangulation.

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
            process.rectify(paired_tile_dir, A_global, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk)
        
        # process.disparity
        if os.path.isfile('%s/rectified_disp.tif' % paired_tile_dir) and cfg['skip_existing']:
            print 'process.disparity on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print 'process.disparity : process tile %d %d (pair %d)...' % (col,row,pair_id)
            process.disparity(paired_tile_dir, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk)
            
        # process.triangulate
        if os.path.isfile('%s/height_map.tif' % paired_tile_dir) and cfg['skip_existing']:
            print 'Triangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
        else:
            print 'Triangulation : process tile %d %d (pair %d)...' % (col,row,pair_id)    
            process.triangulate(paired_tile_dir, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk, np.loadtxt(A_global))

        return '%s/height_map.tif' % paired_tile_dir
        
        

def process_tile(tile_dir,tilesFullInfo):
    """
    Processes a tile : compute the height maps from the N pairs, and finalize the processing, ie.:
    1) Produce a merged height map, without overlapping areas, 
    2) And produces a ply file.

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
    process.finalize_tile(tile_dir, height_maps, tilesFullInfo)

# ----------------------------------------------------------------------------------------------------
# ------------------------------------------- process_tiles ------------------------------------------
# ----------------------------------------------------------------------------------------------------
			 
        

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------- global finalization --------------------------------------
# ----------------------------------------------------------------------------------------------------   

def global_finalization(tilesFullInfo):
    """
    Merges pieces of data into single VRT files : height map comprising the N pairs, height map for each signle pair, and err_rpc.
    Writes the DSM, from the ply files given by each tile.
    Crops corresponding areas in the secondary images.
    Copies the RPC'. 

    Args:
         - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]

    """

    globalfinalization.write_vrt_files(tilesFullInfo)
    globalfinalization.write_dsm(tilesFullInfo)

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
def map_processing(config_file,steps):
    """
    Initialization + preprocessing + global_values + processing + global finalization

    Args: 
         - json config file
    """
    
    try:
    
        if "init_dirs" in steps:
			initialization.init_dirs_srtm_roi(config_file)
    
        tilesFullInfo = initialization.init_tilesFullInfo(config_file)
    
        if cfg['debug']: #monoprocessing
    
            if "preprocess_tiles" in steps:
                print 'preprocess_tiles...'
                for tile_dir in tilesFullInfo:
                    preprocess_tile(tile_dir, tilesFullInfo)
            
            if "global_values" in steps:
                print 'global values...'
                global_values(tilesFullInfo)
            
            if "process_tiles" in steps:
                print 'process_tiles...'
                for tile_dir in tilesFullInfo:
                    process_tile(tile_dir, tilesFullInfo)
                
            if "global_finalization" in steps:    
                print 'global finalization...'     
                global_finalization(tilesFullInfo)
            
        else: # multiprocessing
        
            # create pool with less workers than available cores
            nb_workers = multiprocessing.cpu_count()
            if cfg['max_nb_threads']:
                nb_workers = min(nb_workers, cfg['max_nb_threads'])
            
            
            results = []
            show_progress.counter = 0
            pool = multiprocessing.Pool(nb_workers)
            
            if "preprocess_tiles" in steps:
                print 'preprocess_tile...'
                for tile_dir in tilesFullInfo:
                  
                    p = pool.apply_async(preprocess_tile,
                                                 args=(tile_dir, tilesFullInfo), callback=show_progress)
                    results.append(p)
                    
                for r in results:
                    try:
                        r.get(3600)  # wait at most one hour per tile
                    except multiprocessing.TimeoutError:
                        print "Timeout while computing tile "+str(r)       
                    
            if "global_values" in steps:        
                print 'global values...'        
                global_values(tilesFullInfo)        
                    
            if "process_tiles" in steps:        
                print 'process_tiles...'
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
            
            if "global_finalization" in steps:
                print 'global finalization...'        
                global_finalization(tilesFullInfo)
    
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
        
def execute_job(config_file,tile_dir,command):
    """
    Execute a job

    Args: 
         - json config file
         - job <==> [tile_dir,step]
    
    """
    #tile_dir = job.split(' ')[0]
    #command = job.split(' ')[1]
    
    try:
    
        tilesFullInfo = initialize(config_file)
    
        if command == "preprocess_tiles":
            print 'preprocess_tiles on %s ...' %tile_dir
            preprocess_tile(tile_dir, tilesFullInfo)
        
        if command == "global_values":
            print 'global values...'
            global_values(tilesFullInfo)
        
        if command == "process_tiles" :
            print 'process_tiles on %s ...' %tile_dir
            process_tile(tile_dir, tilesFullInfo)
            
        if command == "global_finalization":    
            print 'global finalization...'     
            global_finalization(tilesFullInfo)  
                
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
            
        
def list_jobs(config_file,step):

    tilesFullInfo = initialize(config_file)
    filename = step[0] + ".jobs"
    
    if step[0] in ["preprocess_tiles","process_tiles"]:
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        for tile_dir in tilesFullInfo:
            f.write(tile_dir + ' ' + step[0] + '\n')
        f.close()
    elif step[0] in ["global_values","global_finalization"]:
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        f.write('all_tiles ' + step[0] + '\n')
        f.close()
    else:
        print "Unkown step required: %s" %step[0]
# ----------------------------------------------------------------------------------------------------
# ------------------------------------- map_processing (end) -----------------------------------------
# ----------------------------------------------------------------------------------------------------


def main(config_file,steps=None,running_mode=None,job=None):
    """
    Launches s2p with the parameters given by a json file.

    Args:
        config_file: path to the config json file
    """


    # measure total runtime
    t0 = time.time()


    ###########################
    if running_mode in ['all','run']:
        map_processing(config_file,steps)
    elif running_mode == 'job':
        execute_job(config_file,job[0],job[1])    
    elif running_mode == 'list_jobs':
        list_jobs(config_file,steps) 
    ###########################                
        

    # runtime
    t = int(time.time() - t0)
    h = t/3600
    m = (t/60) % 60
    s = t % 60
    print "Total runtime: %dh:%dm:%ds" % (h, m, s)
    common.garbage_cleanup()


if __name__ == '__main__':

    possible_steps = ["init_dirs", "preprocess_tiles","global_values","process_tiles","global_finalization"]

    error = False
    if len(sys.argv) == 2 and sys.argv[1].endswith(".json"):
        main(sys.argv[1],possible_steps,'all')
        
    elif len(sys.argv) > 3 and sys.argv[1] != "job":
    
        if sys.argv[1] not in ['run','list_jobs','job']:
            error = True
        for step in sys.argv[3:]:
                if step not in possible_steps:
                    print "Unkown step required: %s" %step
                    error = True    
        
        if not error:
                main(sys.argv[2],sys.argv[3:],sys.argv[1])
                
    elif len(sys.argv) == 5 and sys.argv[1] == "job":
        main(sys.argv[2],None,'job',sys.argv[3:]) 
                
    else:
        error=True
         
         
    if error:
        print """
        Incorrect syntax, use:
          > %s config.json

          Launches the s2p pipeline. All the parameters, paths to input and
          output files, are defined in the json configuration file.

         > %s run config.json [init_dir preprocess_tiles global_values process_tiles global_finalization]
          
          Run specific steps of the s2p pipeline, while skipping the remaining ones.

        > %s job config.json tile_dir command

          Run a specific job defined by a json string. This mode allows to run jobs returned
          by the list_jobs running mode in configuration file.
          
        > %s list_jobs config.json [preprocess_tiles global_values process_tiles global_finalization]
        
          Return the list of jobs for a specific step.
          
        """ % (sys.argv[0],sys.argv[0],sys.argv[0],sys.argv[0])
        sys.exit(1)
