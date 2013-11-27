#!/usr/bin/env python

import numpy as np
from python import common
from python import projective_model
from python import rectification_std as rectification
from python import block_matching
from python import triangulation

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

    # don't use pleiades unsharpening filter
    global_params.use_pleiades_unsharpening = False

except ImportError:
    pass

#### DIFFERENCES: 
# NO SRTM 
# NO RPC

def compute_point_cloud(crop_colorized, heights, P, H, cloud):
    """
    Computes a color point cloud from a height map.

    Args:
        crop_colorized: path to the colorized rectified crop
        heights: height map. Its size is the same as the crop_color image
        P: path to file containing projection matrix 
        H: path to the file containing the coefficients of the rectifying
            homography
        cloud: path to the output points cloud (ply format)
    """
    common.run("colormesh_projective %s %s %s %s %s" % (crop_colorized, heights, P, H,
        cloud))
    return

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
        height_map: path to the height_map, produced by the process_pair of
            process_triplet function
        reference_image_id: id (1, 2 or 3) of the image used as the reference
            image. The height map has been resampled on its grid.
    """


    rpc = 'data/%s/%04d.png.P' % (img_name, reference_image_id)
    im = 'data/%s/%04d.png' % (img_name, reference_image_id)
    im_color = 'data/%s/%04d.png' % (img_name, reference_image_id)   
    crop   = '/tmp/%s_roi_ref%02d.tif' % (exp_name, reference_image_id)
    crop_color = '/tmp/%s_roi_color_ref%02d.tif' % (exp_name, reference_image_id)
    cloud   = '/tmp/%s_cloud.ply'  % (exp_name)



    # read the zoom value
    zoom = global_params.subsampling_factor

    # colorize, then generate point cloud
    tmp_crop = common.image_crop_TIFF(im, x, y, w, h)
    tmp_crop = common.image_safe_zoom_fft(tmp_crop, zoom)
    common.run('cp %s %s' % (tmp_crop, crop))
    A = common.matrix_translation(-x, -y)
    f = 1.0/zoom
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = common.tmpfile('.txt')
    np.savetxt(trans, A)
    
#    compute_point_cloud(common.image_qauto(crop),
#            height_map, rpc, trans, cloud)
    sz = common.image_size(crop)
    compute_point_cloud(common.image_qauto(crop),
          common.image_crop(height_map,0,0,sz[0],sz[1]) , rpc, trans, cloud)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    print "v %s %s %s" % (crop, crop_color, height_map)
    print "meshlab %s" % (cloud)





def main(img_name=None, exp_name=None, x=None, y=None, w=None, h=None,
    reference_image_id=1, secondary_image_id=2):


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


    ## 3. triangulation FOR PROJECTIVE MATRICES DLT algorithm Hartley chapter 12.2 or 12.5
#    from python import disp_to_h_projective as triangulate_proj
#    triangulate_proj.compute_height_map(rpc1,rpc2,hom1,hom2,disp,mask, height, rpc_err)
    common.run("disp_to_h_projective %s %s %s %s %s %s %s %s" % (rpc1, rpc2, hom1, hom2,
        disp, mask, height, rpc_err))

    try:
        zoom = global_params.subsampling_factor
    except NameError:
        zoom = 1
    ref_crop = common.image_crop_TIFF(im1, x, y, w, h)
    triangulation.transfer_map(height, ref_crop, H1, x, y, zoom, height_unrect)
    triangulation.transfer_map(mask, ref_crop, H1, x, y, zoom, mask_unrect)


    ## 4. colorize and generate point cloud
    print (img_name, exp_name, x, y, w, h, height_unrect)
    generate_cloud(img_name, exp_name, x, y, w, h, height_unrect, reference_image_id)



    ### cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())


#    # display results
#    print "v %s %s %s %s" % (rect1, rect2, disp, mask)
#    print "v %s %s %s %s" % (crop1, crop1_color, height_unrect, mask_unrect)
#    print "meshlab %s" % (cloud)



if __name__ == '__main__':


#    # main call: STEREO PAIR WITHOUT ROI
#    main(img_name, exp_name)

#    # main call: STEREO PAIR WITH ROI
#    main(img_name, exp_name, x, y, w, h)

    # main call: TRISTEREO
    img_name = 'fountain324'
    exp_name = '32'
    main(img_name, exp_name, reference_image_id=3, secondary_image_id=2)
    exp_name = '34'
    main(img_name, exp_name, reference_image_id=3, secondary_image_id=4)
    exp_name = '54'
    main(img_name, exp_name, reference_image_id=5, secondary_image_id=4)
    exp_name = '56'
    main(img_name, exp_name, reference_image_id=5, secondary_image_id=6)
#    exp_name = '35'
#    main(img_name, exp_name, reference_image_id=3, secondary_image_id=5)
