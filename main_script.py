#!/usr/bin/env python

import numpy as np
from python import common
from python import rpc_model
from python import pointing_accuracy
from python import rectification
from python import block_matching
from python import triangulation
from python import fusion
from shutil import copyfile

# import the global parameters module
# it permits to pass values between different modules
try:
    from python import global_params

    # zoom factor (determines outputs' resolution)
    global_params.subsampling_factor=1

    # zoom factor used when searching for sift matches
    global_params.subsampling_factor_registration=1

    # matching algorithm: 'tvl1', 'msmw', 'hirschmuller08',
    # hirschmuller08_laplacian'
    global_params.matching_algorithm='hirschmuller08'
#    global_params.matching_algorithm='msmw'

except ImportError:
    pass


def process_pair(img_name, exp_name, x=None, y=None, w=None, h=None,
    reference_image_id=1, secondary_image_id=2):
    """
    Computes a height map from a Pair of Pleiades images.

    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        reference_image_id: id (1, 2 or 3 for the tristereo datasets, and 1 or
            2 for the bistereo datasets) of the image used as the reference
            image of the pair
        secondary_image_id: id of the image used as the secondary image of the
            pair

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """

    # input files
    im1 = 'pleiades_data/images/%s/im%02d.tif' % (img_name, reference_image_id)
    im2 = 'pleiades_data/images/%s/im%02d.tif' % (img_name, secondary_image_id)
    rpc1 = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
    rpc2 = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, secondary_image_id)
    prev1 = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, reference_image_id)

    # output files
    rect1 = '/tmp/%s%d.tif' % (exp_name, reference_image_id)
    rect2 = '/tmp/%s%d.tif' % (exp_name, secondary_image_id)
    outrpc1 = '/tmp/%s_rpc%d.xml' % (exp_name, reference_image_id)
    outrpc2 = '/tmp/%s_rpc%d.xml' % (exp_name, secondary_image_id)
    disp    = '/tmp/%s_disp.pgm'   % (exp_name)
    mask    = '/tmp/%s_mask.png'   % (exp_name)
    height  = '/tmp/%s_height.tif' % (exp_name)
    rpc_err = '/tmp/%s_rpc_err.tif'% (exp_name)
    height_unrect  = '/tmp/%s_height_unrect.tif' % (exp_name)
    mask_unrect    = '/tmp/%s_mask_unrect.png'   % (exp_name)
    subsampling_file = '/tmp/%s_subsampling.txt' % (exp_name)
    pointing = '/tmp/%s_pointing_correction_%02d_%02d.txt' % (exp_name,
        reference_image_id, secondary_image_id)


    ## select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        x, y, w, h = common.get_roi_coordinates(rpc1, prev1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    ## correct pointing error - no subsampling!
    tmp = global_params.subsampling_factor_registration
    global_params.subsampling_factor_registration = 1

    r1 = rpc_model.RPCModel(rpc1)
    r2 = rpc_model.RPCModel(rpc2)
    m = pointing_accuracy.filtered_sift_matches_roi(im1, im2, r1, r2, x, y, w, h)
    A = pointing_accuracy.optimize_pair(im1, im2, r1, r2, None, m)

    global_params.subsampling_factor_registration = tmp

    ## copy the rpcs to the output directory, save the subsampling factor and
    # the pointing correction matrix
    copyfile(rpc1, outrpc1)
    copyfile(rpc2, outrpc2)
    np.savetxt(pointing, A)
    np.savetxt(subsampling_file, np.array([global_params.subsampling_factor]))

    # ATTENTION if subsampling_factor is set the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact

    ## rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1,
        rpc2, x, y, w, h, rect1, rect2, A)

    ## block-matching
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        global_params.matching_algorithm, disp_min, disp_max)

    ## triangulation
    triangulation.compute_height_map(rpc1, rpc2, H1, H2, disp, mask,
        height, rpc_err, A)
    try:
        zoom = global_params.subsampling_factor
    except NameError:
        zoom = 1
    tmp_crop = common.image_crop_TIFF(im1, x, y, w, h)
    ref_crop = common.image_safe_zoom_fft(tmp_crop, zoom)
    triangulation.transfer_map(height, ref_crop, H1, x, y, zoom, height_unrect)
    triangulation.transfer_map(mask, ref_crop, H1, x, y, zoom, mask_unrect)

    ## cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    ## display results
    print "v %s %s %s %s" % (rect1, rect2, disp, mask)

    return height_unrect


def process_triplet(img_name, exp_name, x=None, y=None, w=None, h=None,
    reference_image_id=2, left_image_id=1, right_image_id=3):
    """
    Computes a height map from three Pleiades images.

    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image.  (x, y) is the top-left corner, and (w, h) are the dimensions of
            the rectangle.
        reference_image_id: id (1, 2 or 3) of the image used as the reference
            image of the triplet
        left_image_id: id of the image used as the secondary image of the first
            pair.
        right_image_id: id of the image used as the secondary image of the
            second pair.

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        rpc = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
        prev = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, reference_image_id)
        x, y, w, h = common.get_roi_coordinates(rpc, prev)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # process the two pairs
    exp_name_left = '%s_left' % exp_name
    h_left = process_pair(img_name, exp_name_left, x, y, w, h,
        reference_image_id, left_image_id)

    exp_name_right = '%s_right' % exp_name
    h_right = process_pair(img_name, exp_name_right, x, y, w, h,
        reference_image_id, right_image_id)

    # merge the two height maps
    h = '/tmp/%s_height_merged.tif' % (exp_name)
    fusion.merge(h_left, h_right, 3, h)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    return h


def generate_cloud(img_name, exp_name, x, y, w, h, height_map,
    reference_image_id=1):
    """
    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the original
            panchro image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.
        height_map: path to the height_map, produced by the process_pair or
            process_triplet function
        reference_image_id: id (1, 2 or 3) of the image used as the reference
            image. The height map has been resampled on its grid.
    """

    rpc = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
    im = 'pleiades_data/images/%s/im%02d.tif' % (img_name, reference_image_id)
    im_color = 'pleiades_data/images/%s/im%02d_color.tif' % (img_name, reference_image_id)
    crop   = '/tmp/%s_roi_ref%02d.tif' % (exp_name, reference_image_id)
    crop_color = '/tmp/%s_roi_color_ref%02d.tif' % (exp_name, reference_image_id)
    cloud   = '/tmp/%s_cloud.ply'  % (exp_name)

    # read the zoom value
    zoom = global_params.subsampling_factor

    # build the matrix of the zoom + translation transformation
    A = common.matrix_translation(-x, -y)
    f = 1.0/zoom
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = common.tmpfile('.txt')
    np.savetxt(trans, A)

    # colorize, then generate point cloud
    tmp_crop = common.image_crop_TIFF(im, x, y, w, h)
    tmp_crop = common.image_safe_zoom_fft(tmp_crop, zoom)
    common.run('cp %s %s' % (tmp_crop, crop))
    try:
        with open(im_color):
            triangulation.colorize(crop, im_color, x, y, zoom, crop_color)
            triangulation.compute_point_cloud(crop_color, height_map, rpc,
                trans, cloud)
    except IOError:
        print 'no color image available for this dataset.'
        triangulation.compute_point_cloud(common.image_qauto(crop),
            height_map, rpc, trans, cloud)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    print "v %s %s %s" % (crop, crop_color, height_map)
    print "open %s" % (cloud)


if __name__ == '__main__':

    img_name = 'toulouse'
    exp_name = 'prison'
    x = 20380
    y = 18600
    w = 320
    h = 300

#    w = 700
#    h = 700
#
#    exp_name = 'aeroport'
#    x, y, w, h = 9197, 9630, 1892, 1847
#    img_name = 'calanques'
#    exp_name = 'collines'
#   x = 6600
#   y = 28800
#   w = 1000
#   h = 1000
#
#   img_name = 'cannes'
#   exp_name = 'theoule_sur_mer'
#   x = 5100
#   y = 32300
#   w = 1000
#   h = 1000
#
#    img_name = 'mera'
#    exp_name = 'crete'
#    x = 11127
#    y = 28545
#    w = 800
#    h = 800
#
#   img_name = 'new_york'
#   exp_name = 'manhattan'
#
#    img_name = 'ubaye'
#    exp_name = 'pic'
#    x = 11127
#    y = 13545
#    w = 800
#    h = 800
#
#   img_name = 'uy1'
#   exp_name = 'campo'
#   # FULL ROI
#   #x = 5500
#   #y = 13000
#   #w = 7000
#   #h = 10000
#   # portion inside ROI
#   x = 5500
#   y = 25000
#   w = 1500
#   h = 1500


#    img_name = 'montevideo'
#    exp_name = 'pza_independencia'
#    x, y, w, h = 13025, 26801, 2112, 1496
#    exp_name = 'fing_tvl1'
#    x, y, w, h = 19845, 29178, 1700, 1700

    # main call: STEREO PAIR
    height_map = process_pair(img_name, exp_name, x, y, w, h, 2, 1)
    generate_cloud(img_name, exp_name, x, y, w, h, height_map,
        reference_image_id=2)

    # main call: TRISTEREO
#    height_map = process_triplet(img_name, exp_name, x, y, w, h)
#    generate_cloud(img_name, exp_name, x, y, w, h, height_map,
#        reference_image_id=2)
