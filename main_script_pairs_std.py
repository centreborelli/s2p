#!/usr/bin/env python

import numpy as np
from python import common
from python import projective_model
from python import rectification_std as rectification
from python import block_matching
from python import triangulation


#### DIFFERENCES: 
# NO SRTM 
# NO RPC


def main(img_name=None, exp_name=None, x=None, y=None, w=None, h=None,
    reference_image_id=1, secondary_image_id=2):


    ## Try to import the global parameters module
    #  it permits to pass values between different modules
    try:
       from python import global_params
    
       global_params.subsampling_factor=4
       global_params.subsampling_factor_registration=4

       # select matching algorithm: 'tvl1', 'msmw', 'hirschmuller08',
       # hirschmuller08_laplacian'
       global_params.matching_algorithm='hirschmuller08'

    except ImportError:
      pass


    # input files
    im1 = 'data/%s/%04d.png' % (img_name, reference_image_id)
    im2 = 'data/%s/%04d.png' % (img_name, secondary_image_id)
    rpc1 = 'data/%s/%04d.png.P' % (img_name, reference_image_id)
    rpc2 = 'data/%s/%04d.png.P' % (img_name, secondary_image_id)
    im1_color = 'data/%s/%04d.png' % (img_name, reference_image_id)   
    prev1 = 'data/%s/%04d.png' % (img_name, reference_image_id)  ### GF: all the same image 
    pointing = 'data/%s/pointing_correction_%02d_%02d.txt' % (img_name, reference_image_id, secondary_image_id)

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
        x,y = 0,0
        w,h = common.image_size(im1)
    #    x, y, w, h = common.get_roi_coordinates(rpc1, prev1)    ### GF: ROI IS THE WHOLE IMAGE
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

    H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, 
        rpc1,rpc2, x, y, w, h, rect1, rect2)


    # save homographies to tmp files
    np.savetxt(hom1, H1)
    np.savetxt(hom2, H2)

    ## 2. block-matching
#    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
#        'hirschmuller08', disp_min, disp_max, extra_params='3')
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        global_params.matching_algorithm, disp_min, disp_max)

    print "MISSING TRIANGULATION FOR PROJECTIVE MATRICES: Hartley 12.2 or 12.5"
    exit(0)

    ## 3. triangulation
    triangulation.compute_height_map(rpc1, rpc2, hom1, hom2, disp, mask, height,
          rpc_err)                                                                       ### GF: THIS SHOULD BE REPLACED
    try:
        zoom = global_params.subsampling_factor
    except NameError:
        zoom = 1
    ref_crop = common.image_crop_TIFF(im1, x, y, w, h)
    triangulation.transfer_map(height, ref_crop, hom1, x, y, zoom, height_unrect)
    triangulation.transfer_map(mask, ref_crop, hom1, x, y, zoom, mask_unrect)

    ## 4. colorize and generate point cloud
    crop1 = common.image_crop_TIFF(im1, x, y, w, h)
    trans1 = common.tmpfile('.txt')
    np.savetxt(trans1, common.matrix_translation(-x, -y))
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


#    # main call: STEREO PAIR WITHOUT ROI
#    main(img_name, exp_name)

#    # main call: STEREO PAIR WITH ROI
#    main(img_name, exp_name, x, y, w, h)

    # main call: TRISTEREO
    img_name = 'fountain324'
#    exp_name = '32'
#    main(img_name, exp_name, reference_image_id=3, secondary_image_id=2)
    exp_name = '34'
    main(img_name, exp_name, reference_image_id=3, secondary_image_id=4)
