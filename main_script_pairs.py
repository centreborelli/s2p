#!/usr/bin/env python

import numpy as np
from python import common
from python import rectification
from python import block_matching
from python import triangulation



def main(img_name=None, exp_name=None, x=None, y=None, w=None, h=None,
    reference_image_id=1, secondary_image_id=2):


    ## Try to import the global parameters module
    #  it permits to pass values between different modules
    try:
       from python import global_params
    
       global_params.subsampling_factor=2
       global_params.subsampling_factor_registration=2

       # select matching algorithm: 'tvl1', 'msmw', 'hirschmuller08',
       # hirschmuller08_laplacian'
       global_params.matching_algorithm='hirschmuller08'

    except ImportError:
      pass


    # input files
    im1 = 'pleiades_data/images/%s/im%02d.tif' % (img_name, reference_image_id)
    im2 = 'pleiades_data/images/%s/im%02d.tif' % (img_name, secondary_image_id)
    rpc1 = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
    rpc2 = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, secondary_image_id)
    im1_color = 'pleiades_data/images/%s/im%02d_color.tif' % (img_name, reference_image_id)
    prev1 = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, reference_image_id)
    pointing = 'pleiades_data/images/%s/pointing_correction_%02d_%02d.txt' % (img_name, reference_image_id, secondary_image_id)

    # output files
    rect1 = '/tmp/%s%d.tif' % (exp_name, reference_image_id)
    rect2 = '/tmp/%s%d.tif' % (exp_name, secondary_image_id)
    hom1  = '/tmp/%s_hom%d.txt' % (exp_name, reference_image_id)
    hom2  = '/tmp/%s_hom%d.txt' % (exp_name, secondary_image_id)
    outrpc1 = '/tmp/%s_rpc%d.xml' % (exp_name, reference_image_id)
    outrpc2 = '/tmp/%s_rpc%d.xml' % (exp_name, secondary_image_id)
    crop1_color = '/tmp/%s%d_color.tif' % (exp_name, reference_image_id)
    disp    = '/tmp/%s_disp.pgm'   % (exp_name)
    mask    = '/tmp/%s_mask.png'   % (exp_name)
    cloud   = '/tmp/%s_cloud.ply'  % (exp_name)
    height  = '/tmp/%s_height.tif' % (exp_name)
    rpc_err = '/tmp/%s_rpc_err.tif'% (exp_name)
    height_unrect  = '/tmp/%s_height_unrect.tif' % (exp_name)
    mask_unrect    = '/tmp/%s_mask_unrect.png'   % (exp_name)
    subsampling_file = '/tmp/%s_subsampling.txt' % (exp_name)


    """
    Launches the s2p stereo pipeline on a pair of Pleiades images
    """
    ## 0. select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except (NameError,TypeError):
        x, y, w, h = common.get_roi_coordinates(rpc1, prev1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    ## 0.5 copy the rpcs to the output directory, and save the subsampling factor
    from shutil import copyfile
    copyfile(rpc1, outrpc1)
    copyfile(rpc2, outrpc2)
    np.savetxt(subsampling_file, np.array([global_params.subsampling_factor]))

    # ATTENTION if subsampling_factor is set the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact

    ## 1. rectification
    # If the pointing correction matrix is available, then use it. If not
    # proceed without correction
    try:
        with open(pointing):
            A = np.loadtxt(pointing)
            H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2,
                rpc1, rpc2, x, y, w, h, rect1, rect2, A)
    except IOError:
        H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1,
            rpc2, x, y, w, h, rect1, rect2)

    # save homographies to tmp files
    np.savetxt(hom1, H1)
    np.savetxt(hom2, H2)

    ## 2. block-matching
#    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
#        'hirschmuller08', disp_min, disp_max, extra_params='3')
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        global_params.matching_algorithm, disp_min, disp_max)


    ## 3. triangulation
    triangulation.compute_height_map(rpc1, rpc2, hom1, hom2, disp, mask, height,
        rpc_err)
    try:
        zoom = global_params.subsampling_factor
    except NameError:
        zoom = 1
    triangulation.transfer_height_map(height, mask, hom1, rpc1, x, y, w, h, zoom,
        height_unrect, mask_unrect)

    ## 4. colorize and generate point cloud
    crop1 = common.image_crop_TIFF(im1, x, y, w, h)
    trans1 = common.tmpfile('.txt')
    np.savetxt(trans1, common.matrix_translation(-x, -y))
    try:
        with open(im1_color):
            triangulation.colorize(crop1, im1_color, trans1, crop1_color)
            triangulation.compute_point_cloud(crop1_color, height_unrect, rpc1,
                trans1, cloud)
    except IOError:
        print 'no color image available for this dataset.'
        triangulation.compute_point_cloud(common.image_qauto(crop1),
            height_unrect, rpc1, trans1, cloud)

    ### cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())


    # display results
    print "v %s %s %s %s" % (rect1, rect2, disp, mask)
    print "v %s %s %s %s" % (crop1, crop1_color, height_unrect, mask_unrect)
    print "meshlab %s" % (cloud)



if __name__ == '__main__':

#   img_name = 'lenclio'
#   exp_name = 'tournon'
#   x = 15700
#   y = 16400
#   w = 1000
#   h = 1000
#
    img_name = 'toulouse'
    exp_name = 'prison'
    x = 18527
    y = 17865
    w = 2420
    h = 2464
#
#   img_name = 'calanques'
#   exp_name = 'collines'
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
#   img_name = 'mera'
#   exp_name = 'crete'
#   x = 11127
#   y = 28545
#   w = 1886
#   h = 1755
#
#   img_name = 'new_york'
#   exp_name = 'manhattan'
#
#   img_name = 'ubaye'
#   exp_name = 'pic'
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


#    img_name = 'new_york'
#    exp_name = 'manhattan'


#    # main call: STEREO PAIR WITHOUT ROI
#    main(img_name, exp_name)

#    # main call: STEREO PAIR WITH ROI
#    main(img_name, exp_name, x, y, w, h)

    # main call: TRISTEREO
    exp_name = 'prison_sgbm_21'
    main(img_name, exp_name, x, y, w, h, 2, 1)
    exp_name = 'prison_sgbm_23'
    main(img_name, exp_name, x, y, w, h, 2, 3)
