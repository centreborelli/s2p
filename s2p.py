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
from python import initialization
from python import preprocess
from python import globalvalues
from python import process
from python import globalfinalization

global vabcd
vabcd=[]

def show_progress(a):
    """
    Print the number of tiles that have been processed.

    Args:
        a: useless argument, but since this function is used as a callback by
            apply_async, it has to take one argument.
    """
    show_progress.counter += 1
    status = "done {:{fill}{width}} / {} tiles".format(show_progress.counter,
                                                       show_progress.total,
                                                       fill='',
                                                       width=len(str(show_progress.total)))
    if show_progress.counter < show_progress.total:
        status += chr(8) * len(status)
    else:
        status += '\n'
    sys.stdout.write(status)
    sys.stdout.flush()


def print_elapsed_time(since_first_call=False):
    """
    Print the elapsed time since the last call or since the first call.

    Args:
        since_first_call:
    """
    t2 = datetime.datetime.now()
    if since_first_call:
        print "Total elapsed time:", t2 - print_elapsed_time.t0
    else:
        try:
            print "Elapsed time:", t2 - print_elapsed_time.t1
        except AttributeError:
            print t2 - print_elapsed_time.t0
    print_elapsed_time.t1 = t2


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

    img1, rpc1 = images[0]['img'], images[0]['rpc']
    img2, rpc2 = images[pair_id]['img'], images[pair_id]['rpc']

    out_dir = os.path.join(tile_dir, 'pair_%d' % pair_id)



    A_global = os.path.join(cfg['out_dir'],
                            'global_pointing_pair_%d.txt' % pair_id)

    print 'processing tile %d %d...' % (col, row)

    # check that the tile is not masked
    if os.path.isfile(os.path.join(out_dir, 'this_tile_is_masked.txt')):
        print 'tile %s already masked, skip' % out_dir
        return

    # rectification
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'disp_min_max.txt')) and
        os.path.isfile(os.path.join(out_dir, 'rectified_ref.tif')) and
        os.path.isfile(os.path.join(out_dir, 'rectified_sec.tif'))):
        print '\trectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\trectifying tile %d %d (pair %d)...' % (col, row, pair_id)
        process.rectify(out_dir, np.loadtxt(A_global), img1, rpc1,
                        img2, rpc2, col, row, tw, th, None)

    # disparity estimation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'rectified_mask.png')) and
        os.path.isfile(os.path.join(out_dir, 'rectified_disp.tif'))):
        print '\tdisparity estimation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\testimating disparity on tile %d %d (pair %d)...' % (col, row, pair_id)
        process.disparity(out_dir, img1, rpc1, img2, rpc2, col, row,
                          tw, th, None)

    # triangulation
    if (cfg['skip_existing'] and
        os.path.isfile(os.path.join(out_dir, 'height_map.tif'))):
        print '\ttriangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
    else:
        print '\ttriangulating tile %d %d (pair %d)...' % (col, row, pair_id)
        process.triangulate(out_dir, img1, rpc1, img2, rpc2, col,
                            row, tw, th, None, np.loadtxt(A_global))


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

	##### SOMEHOW THIS IS ABSOLUTELY NECESSARY! IT CREATES SOMETHING MAGIC THAT ALLOWS TO RUN global_align
        # finalization
        height_maps = []
        for i in xrange(nb_pairs):
            if not os.path.isfile(os.path.join(tile_dir, 'pair_%d' % (i+1), 'this_tile_is_masked.txt')):
                height_maps.append(os.path.join(tile_dir, 'pair_%d' % (i+1), 'height_map.tif'))
        process.finalize_tile(tile_info, height_maps, cfg['utm_zone'], cfg['ll_bbx'])

        # ply extrema
        common.run("plyextrema {} {}".format(tile_dir, os.path.join(tile_dir, 'plyextrema.txt')))

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

def process_tile_fusion(tile_info):
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
        nb_pairs = tile_info['number_of_pairs']

        # check that the tile is not masked
        if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
            print 'tile %s already masked, skip' % tile_dir
            return

        # finalization
        height_maps = []
        for i in xrange(nb_pairs):
            if not os.path.isfile(os.path.join(tile_dir, 'pair_%d' % (i+1), 'this_tile_is_masked.txt')):
                height_maps.append(os.path.join(tile_dir, 'pair_%d' % (i+1), 'height_map.tif'))
        process.finalize_tile(tile_info, height_maps, cfg['utm_zone'], cfg['ll_bbx'])

        # ply extrema
        common.run("plyextrema {} {}".format(tile_dir, os.path.join(tile_dir, 'plyextrema.txt')))

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


def apply_global_alignment(tile_info):
    """
    Apply global alignment

    Args:
        tile_info: a dictionary that provides all you need to process a tile
    """
    tile_dir = tile_info['directory']

    # check that the tile is not masked
    if os.path.isfile(os.path.join(tile_dir, 'this_tile_is_masked.txt')):
        print 'tile %s already masked, skip' % tile_dir
        return

    # process each pair to get a height map
    nb_pairs = int(tile_info['number_of_pairs'])
    for pair_id in range(nb_pairs):
    	global vabcd
        abcd = vabcd[pair_id]# = tile_info['alignment_correction_parameters'][pair_id]
        tile_dir = tile_info['directory']
        out_dir = os.path.join(tile_dir, 'pair_%d' %(pair_id + 1))
        fname_h = os.path.join(out_dir, 'height_map.tif')
        cmd = 'plambda %s \"%f * %f +\" -o %s' % (fname_h, abcd[0], abcd[3], fname_h)
        common.run(cmd)


def global_extent(tiles_full_info):
    """
    Compute the global extent from the extrema of each ply file
    """
    xmin, xmax, ymin, ymax = float('inf'), -float('inf'), float('inf'), -float('inf')

    for tile in tiles_full_info:
        plyextrema_file = os.path.join(tile['directory'], 'plyextrema.txt')

        if (os.path.exists(plyextrema_file)):
            extremaxy = np.loadtxt(plyextrema_file)
            xmin = min(xmin, extremaxy[0])
            xmax = max(xmax, extremaxy[1])
            ymin = min(ymin, extremaxy[2])
            ymax = max(ymax, extremaxy[3])

    global_extent = [xmin, xmax, ymin, ymax]
    np.savetxt(os.path.join(cfg['out_dir'], 'global_extent.txt'), global_extent,
               fmt='%6.3f')


def global_align(tiles_full_info):
    """
    Align the N-plets vertically to remove the affine bias due to the pointing error, 
    computed for one tile from N image pairs.
    Uses the first height map as a reference. Stores the old (and biased) height maps 
    as height_map_bias_ and overwrites the height_map_

    Args :
         - height_maps : list of height map directories
         - tile_dir : directory of the tile from which to get a merged height map
         - garbage : a list used to remove temp data (default = [], first call)

    Return value :
        - list of N 4-tuples describing the (a,b,c,d) correction for each pair
    """
    from osgeo import gdal,ogr
    ret = [np.array([1,0,0,0])]

    # first assemble the full height maps associated to each pair
    globalfinalization.write_vrt_files(tiles_full_info)

    nb_pairs = tiles_full_info[0]['number_of_pairs']

    if nb_pairs == 1:
       return ret


    # 1. read the reference height map and create the ij grid
    reference_height_map = cfg['out_dir'] + '/heightMap_pair_%d.vrt'%1

    if not os.path.isfile(reference_height_map ):
        print reference_height_map
        print("the VRT file is supposed to be there after calling globalfinalization.write_vrt_files but it's not!")
        exit()

    hhd = gdal.Open(reference_height_map)
    hhr = hhd.GetRasterBand(1)
    hh  = hhr.ReadAsArray().copy()
    hhd = None
    XX, YY = np.meshgrid(range(hh.shape[1]),range(hh.shape[0]))

    # 2. if the reference height map is almost entirely NAN skip the process? TODO: or choose a new reference? 
    #if np.sum(np.isfinite(hh) < ): 
    #    pass

    #from python import piio
    #piio.write(reference_height_map+'.tif', hh)

    # 3. for each remaining pair of height maps
    for i in range(2, nb_pairs + 1):
        height_map = cfg['out_dir'] + '/heightMap_pair_%d.vrt'%i

        # 3.1 read the secondary height map  
        hhd2 = gdal.Open(height_map)
        hhr2 = hhd2.GetRasterBand(1)
        hh2 = hhr2.ReadAsArray().copy()
	hhd2 = None

        # 3.2 use only the non-nan points in both maps
        mask = np.isfinite(hh) & np.isfinite(hh2) 
        HH  = hh[mask]
        HH2 = hh2[mask] 
        XX2 = XX[mask]
        YY2 = YY[mask]

        print HH.mean(), HH2.mean()


        def solve_irls(X,Y,iter=100):
           ''' 
           solves argmin_A ||A X - Y|| with irls
           for more information about IRLS and sparsity:
           http://www.ricam.oeaw.ac.at/people/page/fornasier/DDFG14.pdf
           
           initialize W as identity matrix
           iterate: 
              argmin_A ||A X W - Y W||^2
              A = Y W X^T inv(X W X^T)
              W = diag( ||AX - Y|| )^(p-2) 
           '''

           # initialize with least squares
           pesos = np.ones(X.shape[1])

           p=0.8
           for t in range(iter):
              W = pesos[:,np.newaxis]
              WXt = W*(X.transpose())

              A = Y.dot(WXt).dot( np.linalg.inv(X.dot(WXt)) )
              r = Y - A.dot(X)
              pesos = np.sqrt(np.sum(r*r,axis=0))**(p-2)
           return A, pesos

        # 4. compute the transformation
        # [ h2 X Y 1 ] alpha = h 
        ##### also implement data centering to improve numerical stability
        #Yo = np.hstack([xy1, h1m[:, np.newaxis]]).transpose()
        #Xo = np.hstack([xy2, h2m[:, np.newaxis], np.ones(( xy1.shape[0] ,1))]).transpose()
        ## center the data for numerical stability
        #center = Xo.mean(axis=1); center[3] = 0;
        #X = Xo - center[:,np.newaxis]
        #Y = Yo - center[0:3,np.newaxis]

        #Xo = np.vstack( [ HH2, XX2*0, YY2*0,  np.ones(( XX2.shape[0] ,1)).squeeze() ])
        #Yo = HH
        ### center the data for numerical stability
        ##center = Xo.mean(axis=1); center[3] = 0;
        ##X = Xo - center[:,np.newaxis]
        ##Y = Yo - center[0,np.newaxis]

        ##   solves argmin_A ||A X - Y||^2
        #alpha = np.linalg.lstsq(Xo.transpose(), Yo)[0]
        #alpha,W = solve_irls(Xo,Yo)
        #assert(alpha[1] == 0)
        #assert(alpha[2] == 0)
        alpha = np.array([1,0,0, HH.mean() - HH2.mean() ])

        ret.append(alpha)

        #hh2 = hh2 * alpha[0] + XX*alpha[1] + YY*alpha[2] + alpha[3]
        #from python import piio
        #piio.write(height_map+'.tif', hh2)

    print ret
    return ret


def compute_dsm(args):
    """
    Compute the DSMs from ply files

    Args:
         - args  ( <==> [config_file,number_of_tiles,current_tile])

    Files with input data (assumed to exist):
        ${out_dir}/tile_*_row_*/col_*/plyextrema.ply
        ${out_dir}/tile_*_row_*/col_*/cloud.ply

    Files created by this code:
        ${out_dir}/tile_*_row_*/col_*/dsm.tif

    """
    list_of_tiles_dir = os.path.join(cfg['out_dir'],'list_of_tiles.txt')

    config_file, number_of_tiles, current_tile = args

    dsm_dir = os.path.join(cfg['out_dir'],'dsm')
    out_dsm = os.path.join(dsm_dir,'dsm_%d.tif' % (current_tile) )

    extremaxy = np.loadtxt(os.path.join(cfg['out_dir'], 'global_extent.txt'))

    global_xmin,global_xmax,global_ymin,global_ymax = extremaxy

    global_y_diff = global_ymax-global_ymin
    tile_y_size = (global_y_diff)/(number_of_tiles)

    # horizontal cuts
    ymin = global_ymin + current_tile*tile_y_size
    ymax = ymin + tile_y_size

    # cutting info
    x, y, w, h, z, ov, tw, th, nb_pairs = initialization.cutting(config_file)
    range_y = np.arange(y, y + h - ov, th - ov)
    range_x = np.arange(x, x + w - ov, tw - ov)
    colmin, rowmin, tw, th = common.round_roi_to_nearest_multiple(z, range_x[0], range_y[0], tw, th)
    colmax, rowmax, tw, th = common.round_roi_to_nearest_multiple(z, range_x[-1], range_y[-1], tw, th)
    cutsinf = '%d %d %d %d %d %d %d %d' % (rowmin, th - ov, rowmax, colmin, tw - ov, colmax, tw, th)

    flags = {}
    flags['average-orig'] = 0
    flags['average'] = 1
    flags['variance'] = 2
    flags['min'] = 3
    flags['max'] = 4
    flags['median'] = 5
    flag = "-flag %d" % (flags.get(cfg['dsm_option'], 0))

    if (ymax <= global_ymax):
        common.run("plytodsm %s %f %s %f %f %f %f %s %s" % (flag,
                                                            cfg['dsm_resolution'],
                                                            out_dsm,
                                                            global_xmin,
                                                            global_xmax, ymin,
                                                            ymax, cutsinf,
                                                            cfg['out_dir']))


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
    globalfinalization.write_dsm()

    # whole point cloud (LidarViewer format)
    if common.which('LidarPreprocessor'):
        out = os.path.join(cfg['out_dir'], 'cloud.lidar_viewer')
        plys = []
        for tile_info in tiles_full_info:
            plys.append(os.path.join(os.path.abspath(tile_info['directory']),
                                     'cloud.ply'))
        globalfinalization.lidar_preprocessor(out, plys)

    # copy RPC xml files in the output directory
    for img in cfg['images']:
        shutil.copy2(img['rpc'], cfg['out_dir'])


def launch_parallel_calls(fun, list_of_args, nb_workers, extra_args=None):
    """
    Run a function several times in parallel with different given inputs.

    Args:
        fun: function to be called several times in parallel.
        list_of_args: list of (first positional) arguments passed to fun, one
            per call
        nb_workers: number of calls run simultaneously
        extra_args (optional, default is None): tuple containing extra arguments
            to be passed to fun (same value for all calls)
    """
    results = []
    show_progress.counter = 0
    pool = multiprocessing.Pool(nb_workers)
    for x in list_of_args:
        args = (x,) + extra_args if extra_args else (x,)
        results.append(pool.apply_async(fun, args=args, callback=show_progress))

    for r in results:
        try:
            r.get(3600)  # wait at most one hour per call
        except multiprocessing.TimeoutError:
            print "Timeout while running %s" % str(r)
        except common.RunFailure as e:
            print "FAILED call: ", e.args[0]["command"]
            print "\toutput: ", e.args[0]["output"]
        except ValueError as e:
            print traceback.format_exc()
            print str(r)
            pass
        except KeyboardInterrupt:
            pool.terminate()
            sys.exit(1)


    pool.close()
    pool.join()


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
            compute_dsm([config_file,cfg['dsm_nb_tiles'],current_tile])

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
        for tile in tiles_full_info:
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
    print_elapsed_time.t0 = datetime.datetime.now()

    if clusterMode == 'list_jobs':
        list_jobs(config_file, step)
    elif clusterMode == 'job':
        cfg['omp_num_threads'] = 1
        execute_job(config_file,misc)
    else:
        # determine which steps to run
        steps = [step] if step else [1, 2, 3, 4, 4.5, 5, 6, 7]

        # initialization (has to be done whatever the queried steps)
        initialization.init_dirs_srtm(config_file)
        tiles_full_info = initialization.init_tiles_full_info(config_file)

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

        if 4.5 in steps:
            print '\nsplit global alignment...'
            abcd = global_align(tiles_full_info)
            global vabcd
            vabcd=abcd
            print_elapsed_time()

            print '\napply global alignment...'
            launch_parallel_calls(apply_global_alignment, tiles_full_info, nb_workers)
            print_elapsed_time()

            print '\ncreate ply clouds...'
            launch_parallel_calls(process_tile_fusion,tiles_full_info,nb_workers)
            print_elapsed_time()

        if 5 in steps:
            print '\ncomputing global extent...'
            global_extent(tiles_full_info)
            print_elapsed_time()


        if 6 in steps:
            print '\ncompute dsm...'
            args = []
            for i in range(cfg['dsm_nb_tiles']):
                args.append([config_file, cfg['dsm_nb_tiles'], i])
            show_progress.total = cfg['dsm_nb_tiles']
            launch_parallel_calls(compute_dsm, args, nb_workers)
            print_elapsed_time()

        if 7 in steps:
            print '\nglobal finalization...'
            global_finalization(tiles_full_info)
            print_elapsed_time()

    # cleanup
    print_elapsed_time(since_first_call=True)
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
