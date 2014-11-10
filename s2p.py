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
import operator

from python import tee
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
from python.config import cfg


def process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None, A=None):
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
        A (optional, default None): pointing correction matrix. If None, it
            will be estimated by this function.

    Returns:
        nothing
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

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
    print 'tile %d %d running on process %s' % (
        x, y, multiprocessing.current_process())

    # select ROI
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

    # check if the ROI is completely masked (water, or outside the image domain)
    H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
    if masking.cloud_water_image_domain(cwid_msk, w, h, H, rpc1, roi_msk,
                                        cld_msk):
        print "Tile masked by water or outside definition domain, skip"
        open("%s/this_tile_is_masked.txt" % out_dir, 'a').close()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if not cfg['debug']:
            fout.close()
        return

    # correct pointing error
    # A is the correction matrix and m is the list of sift matches
    if A is None:
        A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, x,
                                                    y, w, h)
        if A is not None:
            np.savetxt(pointing, A)
        if m is not None:
            np.savetxt(sift_matches, m)
            np.savetxt(center, np.mean(m[:, 2:4], 0))
            visualisation.plot_matches_pleiades(img1, img2, rpc1, rpc2, m, x, y,
                                                w, h, sift_matches_plot)
    else:
        m = None

    # rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m)

    # block-matching
    if cfg['disp_min'] is not None: disp_min = cfg['disp_min']
    if cfg['disp_max'] is not None: disp_max = cfg['disp_max']
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # intersect mask with the cloud_water_image_domain mask
    ww, hh = common.image_size(rect1)
    masking.cloud_water_image_domain(cwid_msk, ww, hh, H1, rpc1, roi_msk,
                                     cld_msk)
    try:
        masking.intersection(mask, mask, cwid_msk)
        masking.erosion(mask, mask, cfg['msk_erosion'])
    except OSError:
        print "file %s not produced" % mask

    # save the subsampling factor, the rectifying homographies and the
    # disparity bounds.
    # ATTENTION if subsampling_factor is > 1 the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact
    np.savetxt(subsampling, np.array([z]))
    np.savetxt(H_ref, H1)
    np.savetxt(H_sec, H2)
    np.savetxt(disp_min_max, np.array([disp_min, disp_max]))

    # save json file with all the parameters needed to reproduce this tile
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

    return


def safe_process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x=None,
                                  y=None, w=None, h=None, prv1=None,
                                  cld_msk=None, roi_msk=None):
    """
    Safe call to process_pair_single_tile (all exceptions will be
    caught). Arguments are the same. This safe version is used when
    processing with multiprocessing module, which will silent
    exceptions and break downstream tasks.

    Returns:
    path to the height map, resampled on the grid of the reference image
    """
    try:
        process_pair_single_tile(out_dir, img1, rpc1, img2, rpc2, x, y, w, h,
                                 prv1, cld_msk, roi_msk)
    # Catch all possible exceptions here
    except:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stdout__
        print sys.exc_info()
        print "failed to generate tile %i %i %i %i\n" % (x, y, w, h)
    return


def process_pair(out_dir, img1, rpc1, img2, rpc2, x, y, w, h, tw=None, th=None,
                 ov=None, cld_msk=None, roi_msk=None):
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
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # duplicate stdout and stderr to log file
    log = tee.Tee('%s/stdout.log' % out_dir, 'w')

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
    if tw is None and th is None and ov is None:
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

    print 'tiles size is tw, th = (%d, %d)' % (tw, th)
    print 'number of tiles in each dimension is %d, %d' % (ntx, nty)
    print 'total number of tiles is %d' % (ntx * nty)

    # create pool with less workers than available cores
    max_processes = multiprocessing.cpu_count()
    if cfg['max_nb_threads'] > 0:
        max_processes = min(max_processes, cfg['max_nb_threads'])

    print 'Creating pool with %d processes\n' % max_processes
    pool = multiprocessing.Pool(max_processes)

    # process the tiles
    # don't parallellize if in debug mode
    tiles = []
    for row in np.arange(y, y + h - ov, th - ov):
        for col in np.arange(x, x + w - ov, tw - ov):
            tile_dir = '%s/tile_%06d_%06d_%04d_%04d' % (out_dir, col, row, tw,
                    th)
            tiles.append(tile_dir)

            # check if the tile is already done, or masked
            if os.path.isfile('%s/rectified_disp.tif' % tile_dir):
                if cfg['skip_existing']:
                    print "stereo on tile %d %d already done, skip" % (col, row)
                    continue
            if os.path.isfile('%s/this_tile_is_masked.txt' % tile_dir):
                print "tile %d %d already masked, skip" % (col, row)
                continue

            # process the tile
            if cfg['debug']:
                process_pair_single_tile(tile_dir, img1, rpc1, img2, rpc2, col,
                                         row, tw, th, None, cld_msk, roi_msk)
            else:
                pool.apply_async(safe_process_pair_single_tile,
                                 args=(tile_dir, img1, rpc1, img2, rpc2, col,
                                       row, tw, th, None, cld_msk, roi_msk))

    # wait for all the processes to terminate
    pool.close()
    pool.join()
    common.garbage_cleanup() 

    # compute global pointing correction
    A_global = pointing_accuracy.global_from_local(tiles)
    np.savetxt('%s/pointing.txt' % out_dir, A_global)

    # Check if all tiles were computed
    # The only cause of a tile failure is a lack of sift matches, which breaks
    # the pointing correction step. Thus it is enough to check if the pointing
    # correction matrix was computed.
    for i, row in enumerate(np.arange(y, y + h - ov, th - ov)):
        for j, col in enumerate(np.arange(x, x + w - ov, tw - ov)):
            tile_dir = '%s/tile_%06d_%06d_%04d_%04d' % (out_dir, col, row, tw,
                    th)
            if not os.path.isfile('%s/this_tile_is_masked.txt' % tile_dir):
                if not os.path.isfile('%s/pointing.txt' % tile_dir):
                    print "%s retrying pointing corr..." % tile_dir
                    # estimate pointing correction matrix from neighbors, if it
                    # fails use A_global, then rerun the disparity map
                    # computation
                    A = pointing_accuracy.from_next_tiles(tiles, ntx, nty, j, i)
                    if A is None:
                        A = A_global
                    process_pair_single_tile(tile_dir, img1, rpc1, img2, rpc2,
                                             col, row, tw, th, None, cld_msk,
                                             roi_msk, A)
    common.garbage_cleanup() 

    # triangulation
    pool = multiprocessing.Pool(max_processes)
    for row in np.arange(y, y + h - ov, th - ov):
        for col in np.arange(x, x + w - ov, tw - ov):
            tile = '%s/tile_%06d_%06d_%04d_%04d' % (out_dir, col, row, tw, th)
            H1 = '%s/H_ref.txt' % tile
            H2 = '%s/H_sec.txt' % tile
            disp = '%s/rectified_disp.tif' % tile
            mask = '%s/rectified_mask.png' % tile
            rpc_err = '%s/rpc_err.tif' % tile
            dem = '%s/dem.tif' % tile

            # check if the tile is already done, or masked
            if os.path.isfile('%s/dem.tif' % tile):
                if cfg['skip_existing']:
                    print "triangulation on tile %d %d is done, skip" % (col,
                                                                         row)
                    continue
            if os.path.isfile('%s/this_tile_is_masked.txt' % tile):
                print "tile %d %d already masked, skip" % (col, row)
                continue

            # process the tile
            if cfg['debug']:
                triangulation.compute_dem(dem, col, row, tw, th, z, rpc1, rpc2,
                                          H1, H2, disp, mask, rpc_err, A_global)
            else:
                pool.apply_async(triangulation.compute_dem, args=(dem, col, row,
                                                                  tw, th, z,
                                                                  rpc1, rpc2,
                                                                  H1, H2, disp,
                                                                  mask, rpc_err,
                                                                  A_global))
    pool.close()
    pool.join()
    common.garbage_cleanup() 

    # tiles composition
    out = '%s/dem.tif' % out_dir
    tmp = ['%s/dem.tif' % t for t in tiles]
    if not os.path.isfile(out) or not cfg['skip_existing']:
        print "Mosaic method: %s" % cfg['mosaic_method']
        if cfg['mosaic_method'] == 'gdal':
            tile_composer.mosaic_gdal(out, w/z, h/z, tmp, tw/z, th/z, ov/z)
        else:
            tile_composer.mosaic(out, w/z, h/z, tmp, tw/z, th/z, ov/z)
    common.garbage_cleanup() 

    return out


def process_triplet(out_dir, img1, rpc1, img2, rpc2, img3, rpc3, x=None, y=None,
                    w=None, h=None, thresh=3, tile_w=None, tile_h=None,
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
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # duplicate stdout and stderr to log file
    log = tee.Tee('%s/stdout.log' % out_dir, 'w')

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        x, y, w, h = common.get_roi_coordinates(rpc1, prv1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # process the two pairs
    out_dir_left = '%s/left' % out_dir
    dem_left = process_pair(out_dir_left, img1, rpc1, img2, rpc2, x, y, w, h,
                            tile_w, tile_h, overlap, cld_msk, roi_msk)

    out_dir_right = '%s/right' % out_dir
    dem_right = process_pair(out_dir_right, img1, rpc1, img3, rpc3, x, y, w, h,
                             tile_w, tile_h, overlap, cld_msk, roi_msk)

    # merge the two digital elevation models
    dem = '%s/dem.tif' % out_dir
    fusion.merge(dem_left, dem_right, thresh, dem)

    common.garbage_cleanup()
    return dem


def generate_cloud(out_dir, im1, rpc1, clr, im2, rpc2, x, y, w, h, dem,
                   do_offset=False):
    """
    Args:
        out_dir: output directory. The file cloud.ply will be written there
        im1:  path to the panchro reference image
        rpc1: path to the xml file containing rpc coefficients for the
            reference image
        clr:  path to the xs (multispectral, ie color) reference image
        im2:  path to the panchro secondary image
        rpc2: path to the xml file containing rpc coefficients for the
            secondary image
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
    crop_color = '%s/roi_color_ref.tif' % out_dir
    crop_ref = '%s/roi_ref.tif' % out_dir
    crop_sec = '%s/roi_sec.tif' % out_dir
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
    trans = '%s/trans.txt' % out_dir
    np.savetxt(trans, A)

    # compute offset
    if do_offset:
        r = rpc_model.RPCModel(rpc1)
        lat = r.latOff
        lon = r.lonOff
        off_x, off_y = geographiclib.geodetic_to_utm(lat, lon)[0:2]
    else:
        off_x, off_y = 0, 0

    # crop the ROI in ref and sec images, then zoom
    if cfg['full_img'] and z == 1:
        crop_ref = im1
        crop_sec = im2
    else:
        r1 = rpc_model.RPCModel(rpc1)
        r2 = rpc_model.RPCModel(rpc2)
        x2, y2, w2, h2 = rpc_utils.corresponding_roi(r1, r2, x, y, w, h)
        if z == 1:
            common.image_crop_TIFF(im1, x, y, w, h, crop_ref)
            common.image_crop_TIFF(im2, x2, y2, w2, h2, crop_sec)
        else:
            # gdal is used for the zoom because it handles BigTIFF files, and
            # before the zoom out the image may be that big
            tmp_crop = common.image_crop_TIFF(im1, x, y, w, h)
            common.image_zoom_gdal(tmp_crop, z, crop_ref, w, h)
            tmp_crop = common.image_crop_TIFF(im2, x2, y2, w2, h2)
            common.image_zoom_gdal(tmp_crop, z, crop_sec, w2, h2)

    # colorize, then generate point cloud
    if clr is not None:
        print 'colorizing...'
        triangulation.colorize(crop_ref, clr, x, y, z, crop_color)
    elif common.image_pix_dim_tiffinfo(crop_ref) == 4:
        print 'the image is pansharpened fusioned'

        # if the image is big, use gdal
        if reduce(operator.mul, common.image_size_tiffinfo(crop_ref)) > 1e8:
            crop_color = common.rgbi_to_rgb_gdal(crop_ref)
            crop_color = common.image_qauto_gdal(crop_color)
        else:
            crop_color = common.rgbi_to_rgb(crop_ref)
            crop_color = common.image_qauto(crop_color)
    else:
        print 'no color data'
        if reduce(operator.mul, common.image_size_tiffinfo(crop_ref)) > 1e8:
            crop_color = common.image_qauto_gdal(crop_ref)
        else:
            crop_color = common.image_qauto(crop_ref)

    triangulation.compute_point_cloud(crop_color, dem, rpc1, trans, cloud,
                                      off_x, off_y)

    common.garbage_cleanup()


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
    if 'img' not in usr_cfg['images'][0] or 'rpc' not in usr_cfg['images'][0]:
        print "missing input data paths for image 0: abort"
        sys.exit(1)
    if 'img' not in usr_cfg['images'][1] or 'rpc' not in usr_cfg['images'][1]:
        print "missing input data paths for image 1: abort"
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


def main(config_file):
    """
    Launches s2p with the parameters given by a json file.

    Args:
        config_file: path to the config json file
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

    # create output directory for the experiment, and store a json dump of the
    # config.cfg dictionary there
    common.run('mkdir -p %s' % cfg['out_dir'])
    f = open('%s/config.json' % cfg['out_dir'], 'w')
    json.dump(cfg, f, indent=2)
    f.close()

    # dem generation
    if len(cfg['images']) == 2:
        dem = process_pair(cfg['out_dir'], cfg['images'][0]['img'],
                           cfg['images'][0]['rpc'], cfg['images'][1]['img'],
                           cfg['images'][1]['rpc'], cfg['roi']['x'],
                           cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h'],
                           None, None, None, cfg['images'][0]['cld'],
                           cfg['images'][0]['roi'])
    else:
        dem = process_triplet(cfg['out_dir'], cfg['images'][0]['img'],
                              cfg['images'][0]['rpc'], cfg['images'][1]['img'],
                              cfg['images'][1]['rpc'], cfg['images'][2]['img'],
                              cfg['images'][2]['rpc'], cfg['roi']['x'],
                              cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h'],
                              cfg['fusion_thresh'], None, None, None, None,
                              cfg['images'][0]['cld'], cfg['images'][0]['roi'])

    # point cloud generation
    generate_cloud(cfg['out_dir'], cfg['images'][0]['img'],
                   cfg['images'][0]['rpc'], cfg['images'][0]['clr'],
                   cfg['images'][1]['img'], cfg['images'][1]['rpc'],
                   cfg['roi']['x'], cfg['roi']['y'], cfg['roi']['w'],
                   cfg['roi']['h'], dem, cfg['offset_ply'])


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
