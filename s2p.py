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

import sys
import shutil
import os.path
import datetime
import traceback
import numpy as np
import multiprocessing

from python.config import cfg
from python import common
from python import rpc_utils
from python import initialization
from python import preprocess
from python import globalvalues
from python import process
from python import globalfinalization


def show_progress(a):
    """
    Print the number of tiles that have been processed.

    Args:
        a: useless argument, but since this function is used as a callback by
            apply_async, it has to take one argument.
    """
    show_progress.counter += 1
    print 'done %d / %d tiles' % (show_progress.counter, show_progress.total)


def print_elapsed_time():
    """
    Print the elapsed time since last call.
    """
    t1 = datetime.datetime.now()
    print "Elapsed time:", t1 - print_elapsed_time.t0
    print_elapsed_time.t0 = t1


def preprocess_tile(tile_info):
    """
    Compute pointing corrections and extrema intensities for a single tile.

    Args:
        tile_info: dictionary containing all the information needed to process a
            tile.
    """
    # create output directory for the tile
    tile_dir = tile_info['directory']
    if not os.path.exists(tile_dir):
        os.makedirs(tile_dir)

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open(os.path.join(tile_dir, 'stdout.log'), 'w', 0)
        # the last arg '0' is for no buffering
        sys.stdout = fout
        sys.stderr = fout

    try:
        preprocess.pointing_correction(tile_info)
        preprocess.minmax_color_on_tile(tile_info)
    except Exception:
        print("Exception in preprocessing tile:")
        traceback.print_exc()
        raise

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()


def global_values(tiles_full_info):
    """
    Compute the global pointing correction and extrema intensities for the ROI.
    """
    globalvalues.pointing_correction(tiles_full_info)
    globalvalues.minmax_intensities(tiles_full_info)


def process_tile_pair(tile_info, pair_id):
    """
    Process a pair of images on a given tile.

    It includes rectification, disparity estimation and triangulation.

    Args:
        tile_info: dictionary containing all the information needed to process a
            tile.
        pair_id: index of the pair to process
    """
    # read all the information
    tile_dir = tile_info['directory']
    col, row, tw, th = tile_info['coordinates']
    images = cfg['images']
    cld_msk = cfg['images'][0]['cld']
    roi_msk = cfg['images'][0]['roi']

    img1, rpc1 = images[0]['img'], images[0]['rpc']
    img2, rpc2 = images[pair_id]['img'], images[pair_id]['rpc']

    out_dir = os.path.join(tile_dir, 'pair_%d' % pair_id)

    A_global = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_%d.txt' % pair_id)

    print 'processing tile %d %d...' % (col, row)
    # rectification
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'rectified_ref.tif')) and
        os.path.isfile(os.path.join(out_dir, 'rectified_sec.tif'))):
        print '\trectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\trectifying tile %d %d (pair %d)...' % (col, row, pair_id)
        process.rectify(out_dir, np.loadtxt(A_global), img1, rpc1,
                        img2, rpc2, col, row, tw, th, None, cld_msk,
                        roi_msk)

    # disparity estimation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'rectified_disp.tif'))):
        print '\tdisparity estimation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\testimating disparity on tile %d %d (pair %d)...' % (col, row, pair_id)
        process.disparity(out_dir, img1, rpc1, img2, rpc2, col, row,
                          tw, th, None, cld_msk, roi_msk)

    # triangulation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'height_map.tif'))):
        print '\ttriangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\ttriangulating tile %d %d (pair %d)...' % (col, row, pair_id)
        process.triangulate(out_dir, img1, rpc1, img2, rpc2, col,
                            row, tw, th, None, cld_msk, roi_msk,
                            np.loadtxt(A_global))


def process_tile(tile_info):
    """
    Process a tile by merging the height maps computed for each image pair.

    Args:
        tile_info: a dictionary that provides all you need to process a tile
    """
    tile_dir = tile_info['directory']

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % tile_dir, 'a', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    try:
        # check that the tile is not masked
        if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
            print 'tile %s already masked, skip' % tile_dir
            return

        # process each pair to get a height map
        nb_pairs = tile_info['number_of_pairs']
        for pair_id in range(1, nb_pairs + 1):
            process_tile_pair(tile_info, pair_id)

        # finalization
        height_maps = [os.path.join(tile_dir, 'pair_%d' % i, 'height_map.tif') for i in range(1, nb_pairs + 1)]
        process.finalize_tile(tile_info, height_maps)
        
        # ply extrema
        out_extrema_dir = os.path.join(tile_dir,'plyextrema.txt')
        common.run("plyextrema %s %s" % (  tile_dir,
                                           out_extrema_dir))

    except Exception:
        print("Exception in processing tile:")
        traceback.print_exc()
        raise

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()


def global_extent(tiles_full_info):
    """
    Compute the global extent from the extrema of each ply file
    """
    xmin,xmax,ymin,ymax = float('inf'),-float('inf'),float('inf'),-float('inf')
    
    for tile in tiles_full_info:
        extremaxy = np.loadtxt(os.path.join(tile['directory'],
                                         'plyextrema.txt'))
            
        xmin=min(xmin,extremaxy[0])
        xmax=max(xmax,extremaxy[1])
        ymin=min(ymin,extremaxy[2])
        ymax=max(ymax,extremaxy[3])
        
    global_extent = [xmin,xmax,ymin,ymax]
    np.savetxt(os.path.join(cfg['out_dir'], 'global_extent.txt'), global_extent,
               fmt='%6.3f') 


def compute_dsm(args):
    """
    Compute the DSMs

    Args: 
         - args  ( <==> [number_of_tiles,current_tile])
    """
    list_of_tiles_dir = os.path.join(cfg['out_dir'],'list_of_tiles.txt')
   
    number_of_tiles,current_tile = args
    
    dsm_dir = os.path.join(cfg['out_dir'],'dsm')
    out_dsm = os.path.join(dsm_dir,'dsm_%d.tif' % (current_tile) )
    
    extremaxy = np.loadtxt(os.path.join(cfg['out_dir'], 'global_extent.txt'))
    
    global_xmin,global_xmax,global_ymin,global_ymax = extremaxy
    
    global_y_diff = global_ymax-global_ymin
    tile_y_size = (global_y_diff)/(number_of_tiles)
    
    # horizontal cuts
    ymin = global_ymin + current_tile*tile_y_size
    ymax = ymin + tile_y_size
    
    if (ymax <= global_ymax):
        common.run("plytodsm %f %s %s %f %f %f %f" % ( 
                                                 cfg['dsm_resolution'], 
                                                 out_dsm, 
                                                 list_of_tiles_dir,
                                                 global_xmin,
                                                 global_xmax,
                                                 ymin,
                                                 ymax))
                                                 
                                             
def global_finalization(tiles_full_info):
    """
    Produce a single height map, DSM and point cloud for the whole ROI.

    The height maps associated to each pair, as well as the height map obtained
    by merging all the pairs, are stored as VRT files. The final DSM is obtained
    by projecting the 3D points from the ply files obtained on each tile. The
    final point cloud is obtained as the union of all the locally merged point
    clouds, in the LidarViewer format.

    Args:
        tiles_full_info: dictionary providing all the information about the
            processed tiles
    """
    globalfinalization.write_vrt_files(tiles_full_info)
    globalfinalization.write_dsm(tiles_full_info)

    # whole point cloud (LidarViewer format)
    if common.which('LidarPreprocessor'):
        out = os.path.join(cfg['out_dir'], 'cloud.lidar_viewer')
        plys = []
        for tile_info in tiles_full_info:
            plys.append(os.path.join(os.path.abspath(tile_info['directory']),
                                     'cloud.ply'))
        globalfinalization.lidar_preprocessor(out, plys)

    # crop area corresponding to the ROI in the secondary images
    if not cfg['full_img']:
        rpc_utils.crop_corresponding_areas(cfg['out_dir'], cfg['images'],
                                           cfg['roi'])

    # copy RPC xml files in the output directory
    for img in cfg['images']:
        shutil.copy2(img['rpc'], cfg['out_dir'])       

def launch_parallel_calls(fun, list_of_args, nb_workers):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel. It must have a
            unique input argument.
        list_of_args: list of arguments passed to fun, one per call
        nb_workers: number of calls run simultaneously
    """
    results = []
    show_progress.counter = 0
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        p = pool.apply_async(fun, args=(x,), callback=show_progress)
        results.append(p)

    for r in results:
        try:
            r.get(3600)  # wait at most one hour per call
        except multiprocessing.TimeoutError:
            print "Timeout while running %s" % str(r)
        except common.RunFailure as e:
            print "FAILED call: ", e.args[0]["command"]
            print "\toutput: ", e.args[0]["output"]
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(1)



def execute_job(config_file,params):
    """
    Execute a job.

    Args:
         - json config file
         - params  ( <==> [tile_dir,step,...])
    """
    tile_dir = params[0]
    step = int(params[1])
    
    tiles_full_info = initialization.init_tiles_full_info(config_file)
    
    if not (tile_dir == 'all_tiles' or 'dsm' in tile_dir ):
        for tile in tiles_full_info:
            if tile_dir == tile['directory']:
                tile_to_process = tile
                break

    try:
        
        if step == 2:#"preprocess_tiles":
            print 'preprocess_tiles on %s ...' % tile_to_process
            preprocess_tile(tile_to_process)
        
        if step == 3:#"global_values":
            print 'global values ...'
            global_values(tiles_full_info)

        if step == 4:#"process_tiles" :
            print 'process_tiles on %s ...' % tile_to_process
            process_tile(tile_to_process)
        
        if step == 5:#"global extent" :
            print 'global extent ...' 
            global_extent(tiles_full_info)
            
        if step == 6:#"compute_dsm" :
            print 'compute_dsm ...'
            current_tile=int(tile_dir.split('_')[1]) # for instance, dsm_2 becomes 2
            compute_dsm([cfg['dsm_nb_tiles'],current_tile])
            
        if step == 7:#"global_finalization":    
            print 'global finalization...'     
            global_finalization(tiles_full_info)  
                
    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]


def list_jobs(config_file, step):

    tiles_full_info = initialization.init_tiles_full_info(config_file)
    filename = str(step) + ".jobs"

    if not (os.path.exists(cfg['out_dir'])):
        os.mkdir(cfg['out_dir'])
    
    if step in [2,4]:           #preprocessing, processing
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        for tile in tilesFullInfo:
            tile_dir = tile['directory']
            f.write(tile_dir + ' ' + str(step) + '\n')
        f.close()
    elif step in [3,5,7]:       # global values, global extent, finalization
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        f.write('all_tiles ' + str(step) + '\n')
        f.close()
    elif step ==6 :             # compute dsm
        f = open(os.path.join(cfg['out_dir'],filename),'w')
        for i in range(cfg['dsm_nb_tiles']):
            f.write('dsm_'+ str(i) + ' ' + str(step) + '\n')
        f.close()
    else:
        print "Unkown step required: %s" % str(step)


def main(config_file, step=None, clusterMode=None, misc=None):
    """
    Launch the entire s2p pipeline with the parameters given in a json file.

    It is a succession of six steps:
        initialization
        preprocessing
        global_values
        processing
        compute dsms
        global_finalization

    Args:
        config_file: path to a json configuration file
        step: integer between 1 and 5 specifying which step to run. Default
        value is None. In that case all the steps are run.
    """    

    if clusterMode == 'list_jobs':
        list_jobs(config_file, step)
    elif clusterMode == 'job':
        cfg['omp_num_threads'] = 1
        execute_job(config_file,misc)
    else:
        # determine which steps to run
        steps = [step] if step else [1, 2, 3, 4, 5, 6, 7]

        # initialization (has to be done whatever the queried steps)
        initialization.init_dirs_srtm(config_file)
        tiles_full_info = initialization.init_tiles_full_info(config_file)
        print_elapsed_time.t0 = datetime.datetime.now()

        # multiprocessing setup
        nb_workers = multiprocessing.cpu_count()  # nb of available cores
        if cfg['max_nb_threads']:
            nb_workers = min(nb_workers, cfg['max_nb_threads'])

        # omp_num_threads: should not exceed nb_workers when multiplied by the
        # number of tiles
        cfg['omp_num_threads'] = max(1, int(nb_workers / len(tiles_full_info)))

        # do the job
        if 2 in steps:
            print '\npreprocessing tiles...'
            show_progress.total = len(tiles_full_info)
            launch_parallel_calls(preprocess_tile, tiles_full_info, nb_workers)
            print_elapsed_time()

        if 3 in steps:
            print '\ncomputing global values...'
            global_values(tiles_full_info)
            print_elapsed_time()

        if 4 in steps:
            print '\nprocessing tiles...'
            show_progress.total = len(tiles_full_info)
            launch_parallel_calls(process_tile, tiles_full_info, nb_workers)
            print_elapsed_time()
           
        if 5 in steps:
            print '\ncomputing global extent...'
            global_extent(tiles_full_info)
            print_elapsed_time()
           
        if 6 in steps:
            print '\ncompute dsm...'
            args=[]
            for i in range(cfg['dsm_nb_tiles']):
                args.append([cfg['dsm_nb_tiles'],i])
            show_progress.total = cfg['dsm_nb_tiles']
            launch_parallel_calls(compute_dsm,args,nb_workers)
            print_elapsed_time()

        if 7 in steps:
            print '\nglobal finalization...'
            global_finalization(tiles_full_info)
            print_elapsed_time()
        
    # cleanup
    common.garbage_cleanup()


if __name__ == '__main__':

    error = False
    steps=[1,2,3,4,5,6,7]

    if len(sys.argv) < 2:
        error = True
        
    elif sys.argv[1].endswith(".json"):
        if len(sys.argv) == 2:
            main(sys.argv[1])
        elif len(sys.argv) == 3 and int(sys.argv[2]) in steps:
            main(sys.argv[1], int(sys.argv[2]))
        else:
            error = True
    else:  # cluster modes
        if sys.argv[1] not in ['list_jobs', 'job']:
            error = True
        else:
            if sys.argv[1] == 'list_jobs':
                if len(sys.argv) == 4 and int(sys.argv[3]) in steps:
                    main(sys.argv[2], int(sys.argv[3]), 'list_jobs')
                else:
                    error = True

            if sys.argv[1] == 'job': 
                if len(sys.argv) >= 5 and int(sys.argv[4]) in steps:
                    main(sys.argv[2], None, 'job', sys.argv[3:])
                else:
                    error = True
    if error:
        print """
        Incorrect syntax, use:
          > %s config.json [step (integer between 1 and 7)]
            1: initialization
            2: preprocessing (tilewise sift, local pointing correction)
            3: global-pointing
            4: processing (tilewise rectification, matching and triangulation)
            5: global-extent
            6: compute dsm from ply files (one per tile)
            7: finalization
            Launches the s2p pipeline.

          > %s list_jobs config.json step (integer between 2 and 7)
            Return the list of jobs for a specific step.

          > %s job config.json tile_dir step (integer between 2 and 7)
            Run a specific job defined by a json string. This mode allows to run jobs returned
            by the list_jobs running mode.


          All the parameters, paths to input and output files, are defined in
          the json configuration file.
          
        """ % (sys.argv[0], sys.argv[0], sys.argv[0])
        sys.exit(1)
