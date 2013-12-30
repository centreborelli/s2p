#!/usr/bin/env python

# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>
# Copyright (C) 2013, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>

import multiprocessing
import sys
import json
import numpy as np

from python import common
from python import rpc_model
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
    # debug print
    print 'tile %d %d, running on process ' % (x, y), multiprocessing.current_process()

    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

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
    A = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, x, y, w,
            h)

    ## save the subsampling factor and
    # the pointing correction matrix
    np.savetxt(pointing, A)
    np.savetxt(subsampling, np.array([z]))

    # ATTENTION if subsampling_factor is set the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact

    ## rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
        rpc2, x, y, w, h, rect1, rect2, A)

    ## block-matching
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        global_params.matching_algorithm, disp_min, disp_max)

    ## triangulation
    if A_global is not None:
        A = A_global
    triangulation.compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height,
        rpc_err, A)
    triangulation.transfer_map(height, H1, x, y, w, h, z, dem)

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
        x = z * np.floor(x / z)
        y = z * np.floor(y / z)
        w = z * np.ceil(w / z)
        h = z * np.ceil(h / z)

    # TODO: automatically compute optimal size for tiles
    # TODO: impose the constraint that ntx*nty is inferior to or equal to a
    # multiple of the number of cores
    if tw is None and th is None and ov is None:
        ov = z * np.ceil(100 / z)
        if w <= 1000:
            tw = w
        else:
            tw = 1000
            #TODO: modify tiles size to be close do a divisor of w
            #while (np.ceil((w - ov) / (tw - ov)) - .2 > (w - ov) / (tw - ov)):
            #    tw += 1
        if h <= 1000:
            th = h
        else:
            th = 1000
            #TODO: modify tiles size to be close do a divisor of h
            #hhile (np.ceil((h - ov) / (th - ov)) - .2 > (h - ov) / (th - ov)):
            #    th += 1
    ntx = np.ceil((w - ov) / (tw - ov))
    nty = np.ceil((h - ov) / (th - ov))
    # ensure that the coordinates of each tile are multiples of the zoom factor
    if (z != 1):
        ov = z * np.floor(ov / z)
        tw = z * np.floor(tw / z)
        th = z * np.floor(th / z)
    print 'tiles size is tw, th = (%d, %d)' % (tw, th)
    print 'number of tiles in each dimension is %d, %d' % (ntx, nty)
    print 'total number of tiles is %d' % (ntx * nty)

    # if several tiles, compute global pointing correction (on the whole ROI)
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
        A = out_dict['correction_matrix']
    else:
        A = None

    # create pool with less workers than available cores
    PROCESSES = int(0.75 * multiprocessing.cpu_count())
    print 'Creating pool with %d processes\n' % PROCESSES
    pool = multiprocessing.Pool(PROCESSES)
    print 'pool = %s' % pool
    print

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
                pool.apply_async(process_pair_single_tile, args=(tile_dir,
                    img1, rpc1, img2, rpc2, i, j, tw, th, A))

    # wait for all the processes to terminate
    pool.close()
    pool.join()

    # tiles composition
    out = '%s/dem.tif' % out_dir
    tile_composer.mosaic(out, w, h, ov, tiles)

    # cleanup
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
    dem = '%s/dem_fusion.tif' % out_dir
    fusion.merge(dem_left, dem_right, thresh, dem)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    return dem


def generate_cloud(out_dir, img, rpc, clr, x, y, w, h, dem, merc_x=None,
        merc_y=None):
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
        merc_{x,y}: mercator coordinates of the point we want to use as
            origin in the local coordinate system of the computed cloud
    """
    # output files
    common.run('mkdir -p %s' % out_dir)
    crop = '%s/roi_ref.tif' % out_dir
    crop_color = '%s/roi_color_ref.tif' % out_dir
    cloud = '%s/cloud.ply' % out_dir

    # read the zoom value
    zoom = global_params.subsampling_factor

    # build the matrix of the zoom + translation transformation
    A = common.matrix_translation(-x, -y)
    f = 1.0/zoom
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = common.tmpfile('.txt')
    np.savetxt(trans, A)

    # compute mercator offset
    if merc_x is None:
        r = rpc_model.RPCModel(rpc)
        lat = r.firstLat
        lon = r.firstLon
        merc_x, merc_y = geographiclib.geodetic_to_mercator(lat, lon)

    # crop the ROI and zoom
    if zoom == 1:
        common.image_crop_TIFF(img, x, y, w, h, crop)
    else:
        tmp_crop = common.image_crop_TIFF(img, x, y, w, h)
        common.image_safe_zoom_fft(tmp_crop, zoom, crop)

    # colorize, then generate point cloud
    try:
        with open(clr):
            triangulation.colorize(crop, clr, x, y, zoom, crop_color)
    except IOError:
        print 'no color image available for this dataset.'
        crop_color = common.image_qauto(crop)

    triangulation.compute_point_cloud(crop_color, dem, rpc, trans, cloud,
            merc_x, merc_y)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())


if __name__ == '__main__':

    if len(sys.argv) == 2:
      config  = sys.argv[1]
    else:
      print """
      Incorrect syntax, use:
        > %s config.json

        Launches the s2p pipeline. All the parameters, paths to input and
        output files, are defined in the json configuration file.
      """ % sys.argv[0]
      sys.exit(1)

    # read the json configuration file
    f = open(config)
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

    # dem generation
    if len(cfg['images']) == 2:
        dem = process_pair(out_dir, img1, rpc1, img2, rpc2, x, y, w, h)
    else:
        img3 = cfg['images'][2]['img']
        rpc3 = cfg['images'][2]['rpc']
        dem = process_triplet(out_dir, img1, rpc1, img2, rpc2, img3, rpc3, x,
                y, w, h, thresh=3)

    # point cloud generation
    generate_cloud(out_dir, img1, rpc1, clr1, x, y, w, h, dem)
