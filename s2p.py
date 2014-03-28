#!/usr/bin/env python

# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>
# Copyright (C) 2013, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>

import multiprocessing
import sys
import json
import numpy as np
import os.path
import copy

from python import tee
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
from python.config import cfg



def process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
        w=None, h=None, A_global=None, prv1=None, cld_msk=None, roi_msk=None):
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
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # redirect stdout and stderr to log file
    if not cfg['debug']:
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
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir

    if os.path.isfile(dem) and cfg['skip_existing']:
        print "Tile %d, %d, %d, %d already generated, skipping" % (x, y, w, h)
        if not cfg['debug']:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            fout.close()
        return dem


    ## select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        x, y, w, h = common.get_roi_coordinates(rpc1, prv1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor
    z = cfg['subsampling_factor']
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

    ## rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
        rpc2, x, y, w, h, rect1, rect2, A, m)

    ## save the subsampling factor, the sift matches, the pointing
    # correction matrix, the rectifying homographies and the disparity bounds
    # ATTENTION if subsampling_factor is set the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact
    np.savetxt(subsampling, np.array([z]))
    np.savetxt(pointing, A)
    np.savetxt(sift_matches, m)
    np.savetxt(H_ref, H1)
    np.savetxt(H_sec, H2)
    np.savetxt(disp_min_max, np.array([disp_min, disp_max]))


    if cfg['disp_range_method'] in ["auto_srtm", "wider_sift_srtm"]:
        # Read models
        rpci1 = rpc_model.RPCModel(rpc1)
        rpci2 = rpc_model.RPCModel(rpc2)
        srtm_disp_min, srtm_disp_max = rpc_utils.rough_disparity_range_estimation(rpci1, rpci2, x, y, w, h,
            H1, H2, A, cfg['disp_range_srtm_low_margin'], cfg['disp_range_srtm_high_margin'])

        if cfg['disp_range_method'] == "auto_srtm":
            disp_min = srtm_disp_min
            disp_max = srtm_disp_max
            print "Auto srtm disp range: [%s, %s]" % (disp_min, disp_max)
        elif cfg['disp_range_method'] == "wider_sift_srtm":
            disp_min = min(srtm_disp_min, disp_min)
            disp_max = max(srtm_disp_max, disp_max)
            print "Wider sift srtm disp range: [%s, %s]" % (disp_min, disp_max)
    else:
         print "Auto sift disp range:  [%s, %s]" % (disp_min, disp_max)

    ## block-matching
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        cfg['matching_algorithm'], disp_min, disp_max)

    ## update mask with cloud mask and roi mask
    if cld_msk is not None:
        triangulation.update_mask(mask, H1, cld_msk, True)
    if roi_msk is not None:
        triangulation.update_mask(mask, H1, roi_msk, False, cfg['msk_erosion'])

    ## triangulation
    if A_global is not None:
        A = A_global
    triangulation.compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height,
        rpc_err, A)
    triangulation.transfer_map(height, H1, x, y, w, h, z, dem)

    ## save json file with all the parameters needed to reproduce this tile
    tile_cfg = copy.deepcopy(cfg)
    tile_cfg['roi']['x'] = x
    tile_cfg['roi']['y'] = y
    tile_cfg['roi']['w'] = w
    tile_cfg['roi']['h'] = h
    f = open(config, 'w')
    json.dump(tile_cfg, f, indent=2)
    f.close()

    # close logs
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return dem


def safe_process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x=None,
        y=None, w=None, h=None, A_global=None, prv1=None, cld_msk=None,
        roi_msk=None):
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
            w, h, A_global, prv1, cld_msk, roi_msk)
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
        h=None, tw=None, th=None, ov=None, cld_msk=None, roi_msk=None):
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
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.

    Returns:
        path to the digital elevation model (dem), resampled on the grid of the
        reference image.
    """
    # duplicate stdout and stderr to log file
    log = tee.Tee('%s/stdout.log' % out_dir, 'w')

    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    z = cfg['subsampling_factor']
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
        if w <= z * cfg['tile_size']:
            tw = w
        else:
            tw = z * cfg['tile_size']
            #TODO: modify tiles size to be close do a divisor of w
            #while (np.ceil((w - ov) / (tw - ov)) - .2 > (w - ov) / (tw - ov)):
            #    tw += 1
        if h <= z * cfg['tile_size']:
            th = h
        else:
            th = z * cfg['tile_size']
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
        matrix_file = '%s/pointing_global.txt' % out_dir

        if not os.path.isfile(matrix_file) or not cfg['skip_existing']:
            p = multiprocessing.Process(target=pointing_accuracy.compute_correction,
                    args=(img1, rpc1, img2, rpc2, x, y, w, h, out_dict))
            p.start()
            p.join()
            if 'correction_matrix' in out_dict:
                A = out_dict['correction_matrix']
                np.savetxt(matrix_file, A)
            else:
                print """WARNING: global correction matrix not found. The
            estimation process seems to have failed. No global correction will
            be applied."""

    # create pool with less workers than available cores
    max_processes = multiprocessing.cpu_count()
    if cfg['max_nb_threads'] > 0:
        max_processes = min(max_processes, cfg['max_nb_threads'])

    print 'Creating pool with %d processes\n' % max_processes
    pool = multiprocessing.Pool(max_processes)

    # process the tiles
    # don't parallellize if in debug mode
    tiles = []
    for j in np.arange(y, y + h - ov, th - ov):
        for i in np.arange(x, x + w - ov, tw - ov):
            tile_dir = '%s/tile_%d_%d_%d_%d' % (out_dir, i, j, tw, th)
            tiles.append('%s/dem.tif' % tile_dir)
            if cfg['debug']:
                process_pair_single_tile(tile_dir, img1, rpc1, img2, rpc2, i,
                    j, tw, th, A, None, cld_msk, roi_msk)
            else:
                pool.apply_async(safe_process_pair_single_tile, args=(tile_dir,
                    img1, rpc1, img2, rpc2, i, j, tw, th, A, None, cld_msk,
                    roi_msk))

    # wait for all the processes to terminate
    pool.close()
    pool.join()


    # Check if all tiles were computed
    for j in np.arange(y, y + h - ov, th - ov):
        for i in np.arange(x, x + w - ov, tw - ov):
            tile_dir = '%s/tile_%d_%d_%d_%d' % (out_dir, i, j, tw, th)
            dem = '%s/dem.tif' % tile_dir

            if not os.path.isfile(dem):
                print "WARNING: Tile %d %d %d %d failed. Retrying..." % (i, j,
                        tw, th)
                process_pair_single_tile(tile_dir, img1, rpc1, img2, rpc2, i,
                        j, tw, th, A, None, cld_msk, roi_msk)

    # tiles composition
    out = '%s/dem.tif' % out_dir
    if not os.path.isfile(out) or not cfg['skip_existing']:
        print "Mosaic method: %s" % cfg['mosaic_method']
        if cfg['mosaic_method'] == 'gdal':
            tile_composer.mosaic_gdal(out, w/z, h/z, tiles, tw/z, th/z, ov/z)
        else:
            tile_composer.mosaic(out, w/z, h/z, tiles, tw/z, th/z, ov/z)

    # cleanup
    if cfg['clean_tmp']:
        while common.garbage:
            common.run('rm ' + common.garbage.pop())

    return out


def process_triplet(out_dir, img1, rpc1, img2, rpc2, img3, rpc3, x=None,
        y=None, w=None, h=None, thresh=3, tile_w=None, tile_h=None,
        overlap=None, prv1=None, cld_msk=None, roi_msk=None):
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
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.

    Returns:
        path to the digital elevaton model (dem), resampled on the grid of the
        reference image.
    """
    # duplicate stdout and stderr to log file
    log = tee.Tee('%s/stdout.log' % out_dir, 'w')

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
            tile_w, tile_h, overlap, cld1, roi1)

    out_dir_right = '%s/right' % out_dir
    dem_right = process_pair(out_dir_right, img1, rpc1, img3, rpc3, x, y, w, h,
            tile_w, tile_h, overlap, cld1, roi1)

    # merge the two digital elevation models
    dem = '%s/dem.tif' % out_dir
    fusion.merge(dem_left, dem_right, thresh, dem)

    # cleanup
    if cfg['clean_tmp']:
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
    z = cfg['subsampling_factor']
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
    except (IOError, TypeError):
        print 'no color image available for this dataset.'
        crop_color = common.image_qauto(crop)

    triangulation.compute_point_cloud(crop_color, dem, rpc, trans, cloud,
            off_x, off_y)

    # cleanup
    if cfg['clean_tmp']:
        while common.garbage:
            common.run('rm ' + common.garbage.pop())


def check_parameters(usr_cfg):
    """
    Checks that the provided dictionary defines all the mandatory
    arguments, and warns about unknown optional arguments.

    Args:
        usr_cfg: python dict read from the json input file

    Returns:
        nothing
    """

    ## verify that i/o files and roi are defined
    # i/o files
    if 'out_dir' not in usr_cfg:
        print "missing output dir: abort"
        sys.exit(1)
    if 'images' not in usr_cfg or len(usr_cfg['images']) < 2:
        print "missing input data paths: abort"
        sys.exit(1)
    if 'img' not in usr_cfg['images'][0] or 'rpc' not in usr_cfg['images'][0]:
        print "missing input data paths for image 0: abort"
        sys.exit(1)
    if 'img' not in usr_cfg['images'][1] or 'rpc' not in usr_cfg['images'][1]:
        print "missing input data paths for image 1: abort"
        sys.exit(1)

    # roi
    if ('full_img' not in usr_cfg) or (not usr_cfg['full_img']):
        if 'roi' not in usr_cfg or ('x' not in usr_cfg['roi']) or ('y' not in
                usr_cfg['roi']) or ('w' not in usr_cfg['roi']) or ('h' not in
                        usr_cfg['roi']):
            print "bad or missing roi definition: abort"
            sys.exit(1)

    ## warn about unknown optional parameters: these parameters have no default
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
        if not k in cfg:
            print """parameter %s unknown: you should remove it from the input
            json file. It will be ignored.""" % k


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
        sz = common.image_size_gdal(cfg['images'][0]['img'])
        cfg['roi'] = {}
        cfg['roi']['x'] = 0
        cfg['roi']['y'] = 0
        cfg['roi']['w'] = sz[0]
        cfg['roi']['h'] = sz[1]

    # create output directory for the experiment, and store a json dump of the
    # config.cfg dictionary there
    common.run('mkdir -p %s' % cfg['out_dir'])
    f = open('%s/config.json' % cfg['out_dir'], 'w')
    json.dump(cfg, f, indent=2)
    f.close()

    # dem generation
    if len(cfg['images']) == 2:
        dem = process_pair(cfg['out_dir'],
            cfg['images'][0]['img'], cfg['images'][0]['rpc'],
            cfg['images'][1]['img'], cfg['images'][1]['rpc'],
            cfg['roi']['x'], cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h'],
            None, None, None, cfg['images'][0]['cld'], cfg['images'][0]['roi'])
    else:
        dem = process_triplet(cfg['out_dir'],
            cfg['images'][0]['img'], cfg['images'][0]['rpc'],
            cfg['images'][1]['img'], cfg['images'][1]['rpc'],
            cfg['images'][2]['img'], cfg['images'][2]['rpc'],
            cfg['roi']['x'], cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h'],
            cfg['fusion_thresh'], None, None, None, None,
            cfg['images'][0]['cld'], cfg['images'][0]['roi'])

    # point cloud generation
    generate_cloud(cfg['out_dir'],
      cfg['images'][0]['img'], cfg['images'][0]['rpc'], cfg['images'][0]['clr'],
      cfg['roi']['x'], cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h'],
      dem, cfg['offset_ply'])
