#!/usr/bin/env python

# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>
# Copyright (C) 2013, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>

import multiprocessing
import sys
import json
import numpy as np
import os.path

from python import common
from python import rpc_model
from python import rpc_utils
from python import geographiclib
from python import pointing_accuracy
from python import rectification
from python import block_matching
from python import triangulation
from python import tile_composer
from python import fusion
from python import global_params



def process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
        w=None, h=None, A_global=None, prv1=None):
    """
    Computes a height map from a Pair of Pleiades images, without tiling

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
        A_global (optional): global pointing correction matrix, used for
            triangulation (but not for stereo-rectification)
        prv1 (optional): path to a preview of the reference image

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # redirect stdout and stderr to log file
    if not global_params.debug:
        fout = open('%s/stdout.log' % out_dir, 'w', 0) # '0' is for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d, running on process ' % (x, y), multiprocessing.current_process()

    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp    = '%s/rectified_disp.pgm'   % (out_dir)
    mask    = '%s/rectified_mask.png'   % (out_dir)
    height  = '%s/rectified_height.tif' % (out_dir)
    rpc_err = '%s/rpc_err.tif'% (out_dir)
    dem     = '%s/dem.tif' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir

    ## select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        x, y, w, h = common.get_roi_coordinates(rpc1, prv1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor
    z = global_params.subsampling_factor
    assert(z > 0 and z == np.floor(z))
    if (z != 1):
        x = z * np.floor(x / z)
        y = z * np.floor(y / z)
        w = z * np.ceil(w / z)
        h = z * np.ceil(h / z)

    ## correct pointing error
    # A is the correction matrix and m is the list of sift matches
    A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, x, y,
        w, h, first_guess=A_global)

    ## save the subsampling factor, the sift matches and the pointing
    # correction matrix
    np.savetxt(subsampling, np.array([z]))
    np.savetxt(pointing, A)
    np.savetxt(sift_matches, m)

    # ATTENTION if subsampling_factor is set the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact

    ## rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
        rpc2, x, y, w, h, rect1, rect2, A, m)


    if global_params.disp_range_method == "auto_srtm" or global_params.disp_range_method == "wider_sift_srtm":
        # Read models
        rpci1 = rpc_model.RPCModel(rpc1)
        rpci2 = rpc_model.RPCModel(rpc2)
        srtm_disp_min, srtm_disp_max = rpc_utils.rough_disparity_range_estimation(rpci1,rpci2,x,y,w,h,H1,H2,A,global_params.disp_range_srtm_low_margin, global_params.disp_range_srtm_high_margin)
        # srtm_disp_min, srtm_disp_max = estimate_disp_range_from_srtm(img1, img2, rpc1, rpc2, H1, H2, x, y, w, h)

        if global_params.disp_range_method == "auto_srtm":
            disp_min = srtm_disp_min
            disp_max = srtm_disp_max
            print "Auto srtm disp range: ["+str(disp_min)+", "+str(disp_max)+"]"
        elif global_params.disp_range_method == "wider_sift_srtm":
            disp_min = min(srtm_disp_min, disp_min)
            disp_max = max(srtm_disp_max, disp_max)
            print "Wider sift srtm disp range: ["+str(disp_min)+", "+str(disp_max)+"]"
    else:
         print "Auto sift disp range: ["+str(disp_min)+", "+str(disp_max)+"]"
    ## block-matching
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        global_params.matching_algorithm, disp_min, disp_max)

    ## triangulation
    if A_global is not None:
        A = A_global
    triangulation.compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height,
        rpc_err, A)
    triangulation.transfer_map(height, H1, x, y, w, h, z, dem)

    # close logs
    if not global_params.debug:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return dem


def safe_process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
        w=None, h=None, A_global=None, prv1=None):
    """
    Safe call to process_pair_single_tile (all exceptions will be
    caught). Arguments are the same. This safe version is used when
    processing with multiprocessing module, which will silent
    exceptions and break downstream tasks.

    Returns:
    path to the height map, resampled on the grid of the reference image
    """
    dem = ""
    try:
        dem = process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x, y,
        w, h, A_global, prv1)
    # Catch all possible exceptions here
    except:
        e = sys.exc_info()[0]
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        print "Failed to generate tile %i %i %i %i: %s)" %(x,y,w,h,str(e))
        # Append to the list of failed tiles
        return dem

    print "Tile %i %i %i %i generated." %(x,y,w,h)
    return dem

def process_pair(out_dir, img1, rpc1, img2, rpc2, x=None, y=None, w=None,
        h=None, tw=None, th=None, ov=None):
    """
    Computes a height map from a Pair of Pleiades images, using tiles.

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
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        tw, th: dimensions of the tiles
        ov: width of overlapping bands between tiles

    Returns:
        path to the digital elevation model (dem), resampled on the grid of the
        reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    z = global_params.subsampling_factor
    assert(z > 0 and z == np.floor(z))
    if (z != 1):
        x = z * np.floor(float(x) / z)
        y = z * np.floor(float(y) / z)
        w = z * np.ceil(float(w) / z)
        h = z * np.ceil(float(h) / z)

    # TODO: automatically compute optimal size for tiles
    # TODO: impose the constraint that ntx*nty is inferior to or equal to a
    # multiple of the number of cores
    if tw is None and th is None and ov is None:
        ov = z * 100
        if w <= z * global_params.tile_size:
            tw = w
        else:
            tw = z * global_params.tile_size
            #TODO: modify tiles size to be close do a divisor of w
            #while (np.ceil((w - ov) / (tw - ov)) - .2 > (w - ov) / (tw - ov)):
            #    tw += 1
        if h <= z * global_params.tile_size:
            th = h
        else:
            th = z * global_params.tile_size
            #TODO: modify tiles size to be close do a divisor of h
            #hhile (np.ceil((h - ov) / (th - ov)) - .2 > (h - ov) / (th - ov)):
            #    th += 1
    ntx = np.ceil(float(w - ov) / (tw - ov))
    nty = np.ceil(float(h - ov) / (th - ov))

    print 'tiles size is tw, th = (%d, %d)' % (tw, th)
    print 'number of tiles in each dimension is %d, %d' % (ntx, nty)
    print 'total number of tiles is %d' % (ntx * nty)

    # if several tiles, compute global pointing correction (on the whole ROI)
    A = None
    if ntx * nty > 1:
        # the global pointing correction is run in a subprocess. This is a
        # workaround to a nasty bug affecting the Multiprocessing package when used
        # with Numpy on osx, causing Python to 'quit unexpectedly':
        # http://stackoverflow.com/questions/19705200/multiprocessing-with-numpy-makes-python-quit-unexpectedly-on-osx
        manager = multiprocessing.Manager()
        out_dict = manager.dict()
        p = multiprocessing.Process(target=pointing_accuracy.compute_correction,
            args=(img1, rpc1, img2, rpc2, x, y, w, h, out_dict))
        p.start()
        p.join()
        if 'correction_matrix' in out_dict:
            A = out_dict['correction_matrix']
            np.savetxt('%s/pointing_global.txt' % out_dir, A)
        else:
            print """WARNING: global correction matrix not found. The
            estimation process seems to have failed. No global correction will
            be applied."""

    # create pool with less workers than available cores
    PROCESSES = int(0.75 * multiprocessing.cpu_count())
    if global_params.max_nb_threads > 0:
        PROCESSES = min(PROCESSES, global_params.max_nb_threads)

    print 'Creating pool with %d processes\n' % PROCESSES
    pool = multiprocessing.Pool(PROCESSES)

    # process the tiles
    # don't parallellize if in debug mode
    tiles = []
    for j in np.arange(y, y + h - ov, th - ov):
        for i in np.arange(x, x + w - ov, tw - ov):
            tile_dir = '%s/tile_%d_%d_%d_%d' % (out_dir, i, j, tw, th)
            tiles.append('%s/dem.tif' % tile_dir)
            if global_params.debug:
                process_pair_single_tile(tile_dir, img1, rpc1, img2, rpc2, i,
                    j, tw, th, A)
            else:
                pool.apply_async(safe_process_pair_single_tile, args=(tile_dir,
                    img1, rpc1, img2, rpc2, i, j, tw, th, A))

    # wait for all the processes to terminate
    pool.close()
    pool.join()


    # Check if all tiles were computed
    for j in np.arange(y, y + h - ov, th - ov):
        for i in np.arange(x, x + w - ov, tw - ov):
            tile_dir = '%s/tile_%d_%d_%d_%d' % (out_dir, i, j, tw, th)
            dem = '%s/dem.tif' % tile_dir

            if not os.path.exists(dem):
                print "WARNING: Tile %d %d %d %d failed. Retrying..." % (i, j,
                        tw, th)
                process_pair_single_tile(tile_dir, img1, rpc1, img2, rpc2, i,
                        j, tw, th, A)

    # tiles composition
    out = '%s/dem.tif' % out_dir
    tile_composer.mosaic(out, w/z, h/z, tiles, tw/z, th/z, ov/z)

    # cleanup
    if global_params.clean_tmp:
        while common.garbage:
            common.run('rm ' + common.garbage.pop())

    return out


def process_triplet(out_dir, img1, rpc1, img2, rpc2, img3, rpc3, x=None,
        y=None, w=None, h=None, thresh=3, tile_w=None, tile_h=None,
        overlap=None, prv1=None):
    """
    Computes a height map from three Pleiades images.

    Args:
        out_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image of the first pair
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image of the first pair
        img3: path to the secondary image of the second pair
        rpc3: paths to the xml file containing the rpc coefficients of the
            secondary image of the second pair
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        thresh: threshold used for the fusion algorithm, in meters.
        tile_w, tile_h: dimensions of the tiles
        overlap: width of overlapping bands between tiles
        prv1 (optional): path to a preview of the reference image

    Returns:
        path to the digital elevaton model (dem), resampled on the grid of the
        reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        x, y, w, h = common.get_roi_coordinates(rpc1, prv1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # process the two pairs
    out_dir_left = '%s/left' % out_dir
    dem_left = process_pair(out_dir_left, img1, rpc1, img2, rpc2, x, y, w, h,
            tile_w, tile_h, overlap)

    out_dir_right = '%s/right' % out_dir
    dem_right = process_pair(out_dir_right, img1, rpc1, img3, rpc3, x, y, w, h,
            tile_w, tile_h, overlap)

    # merge the two digital elevation models
    dem = '%s/dem.tif' % out_dir
    fusion.merge(dem_left, dem_right, thresh, dem)

    # cleanup
    if global_params.clean_tmp:
        while common.garbage:
            common.run('rm ' + common.garbage.pop())

    return dem


def generate_cloud(out_dir, img, rpc, clr, x, y, w, h, dem, do_offset=False):
    """
    Args:
        out_dir: output directory. The file cloud.ply will be written there
        img: path to the panchro image
        rpc: path to the xml file containing rpc coefficients
        clr: path to the xs (multispectral, ie color) image
        x, y, w, h: four integers defining the rectangular ROI in the original
            panchro image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.
        dem: path to the digital elevation model, produced by the process_pair
            or process_triplet function
        do_offset (optional, default: False): boolean flag to decide wether the
            x, y coordinates of points in the ply file will be translated or
            not (translated to be close to 0, to avoid precision loss due to
            huge numbers)
    """
    # output files
    common.run('mkdir -p %s' % out_dir)
    crop = '%s/roi_ref.tif' % out_dir
    crop_color = '%s/roi_color_ref.tif' % out_dir
    cloud = '%s/cloud.ply' % out_dir

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    z = global_params.subsampling_factor
    assert(z > 0 and z == np.floor(z))
    if (z != 1):
        x = z * np.floor(float(x) / z)
        y = z * np.floor(float(y) / z)
        w = z * np.ceil(float(w) / z)
        h = z * np.ceil(float(h) / z)

    # build the matrix of the zoom + translation transformation
    A = common.matrix_translation(-x, -y)
    f = 1.0/z
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = common.tmpfile('.txt')
    np.savetxt(trans, A)

    # compute offset
    if do_offset:
        r = rpc_model.RPCModel(rpc)
        lat = r.latOff
        lon = r.lonOff
        off_x, off_y = geographiclib.geodetic_to_utm(lat, lon)[0:2]
    else:
        off_x, off_y = 0, 0

    # crop the ROI and zoom
    if z == 1:
        common.image_crop_TIFF(img, x, y, w, h, crop)
    else:
        # gdal is used for the zoom because it handles BigTIFF files, and
        # before the zoom out the image may be that big
        tmp_crop = common.image_crop_TIFF(img, x, y, w, h)
        common.image_zoom_gdal(tmp_crop, z, crop, w, h)

    # colorize, then generate point cloud
    try:
        with open(clr):
            triangulation.colorize(crop, clr, x, y, z, crop_color)
    except IOError:
        print 'no color image available for this dataset.'
        crop_color = common.image_qauto(crop)

    triangulation.compute_point_cloud(crop_color, dem, rpc, trans, cloud,
            off_x, off_y)

    # cleanup
    if global_params.clean_tmp:
        while common.garbage:
            common.run('rm ' + common.garbage.pop())


if __name__ == '__main__':

    if len(sys.argv) == 2:
      config_file  = sys.argv[1]
    else:
      print """
      Incorrect syntax, use:
        > %s config.json

        Launches the s2p pipeline. All the parameters, paths to input and
        output files, are defined in the json configuration file.
      """ % sys.argv[0]
      sys.exit(1)

    # read the json configuration file
    f = open(config_file)
    cfg = json.load(f)
    f.close()

    # fill the global_params module
    global_params.subsampling_factor = cfg['subsampling_factor']
    global_params.subsampling_factor_registration = cfg['subsampling_factor_registration']
    global_params.sift_match_thresh = cfg['sift_match_thresh']
    global_params.disp_range_extra_margin = cfg['disp_range_extra_margin']
    global_params.n_gcp_per_axis = cfg['n_gcp_per_axis']
    global_params.epipolar_thresh = cfg['epipolar_thresh']
    global_params.matching_algorithm = cfg['matching_algorithm']
    global_params.use_pleiades_unsharpening = cfg['use_pleiades_unsharpening']
    global_params.debug = cfg['debug']
    if "disp_range_method" in cfg:
        global_params.disp_range_method = str(cfg['disp_range_method'])
    if "disp_range_srtm_low_margin" in cfg:
        global_params.disp_range_srtm_low_margin = float(cfg['disp_range_srtm_low_margin'])
    if "disp_range_srtm_high_margin" in cfg:
        global_params.disp_range_srtm_high_margin = float(cfg['disp_range_srtm_high_margin'])

    if "temporary_dir" in cfg:
        global_params.temporary_dir = str(cfg['temporary_dir'])
    if "tile_size" in cfg:
        global_params.tile_size = cfg['tile_size']
    if "max_nb_threads" in cfg:
        global_params.max_nb_threads = cfg['max_nb_threads']
    if "clean_tmp" in cfg:
        global_params.clean_tmp = cfg['clean_tmp']
    if "fusion_thresh" in cfg:
        global_params.fusion_thresh = cfg['fusion_thresh']

    # other params to be read in the json
    if "offset_ply" in cfg:
        do_offset = cfg['offset_ply']
    else:
        do_offset = False

    # roi definition and output path
    x = cfg['roi']['x']
    y = cfg['roi']['y']
    w = cfg['roi']['w']
    h = cfg['roi']['h']
    out_dir = str(cfg['out_dir'])
    img1 = cfg['images'][0]['img']
    rpc1 = cfg['images'][0]['rpc']
    clr1 = cfg['images'][0]['clr']
    img2 = cfg['images'][1]['img']
    rpc2 = cfg['images'][1]['rpc']

    # create output directory for the experiment, and store a copy the json
    # config file there
    common.run('mkdir -p %s' % out_dir)
    common.run('cp %s %s/config.json' % (config_file, out_dir))

    # dem generation
    if len(cfg['images']) == 2:
        dem = process_pair(out_dir, img1, rpc1, img2, rpc2, x, y, w, h)
    else:
        img3 = cfg['images'][2]['img']
        rpc3 = cfg['images'][2]['rpc']
        dem = process_triplet(out_dir, img1, rpc1, img2, rpc2, img3, rpc3, x,
                y, w, h, global_params.fusion_thresh)

    # point cloud generation
    generate_cloud(out_dir, img1, rpc1, clr1, x, y, w, h, dem, do_offset)
