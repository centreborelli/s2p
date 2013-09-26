#!/usr/bin/env python

import numpy as np
from python import common
from python import rectification
from python import block_matching
from python import triangulation

img_name = 'lenclio'
exp_name = 'beu'
x = 15700
y = 16400
w = 1000
h = 1000

#img_name = 'toulouse'
#exp_name = 'blagnac'
#x = 6000
#y = 6000
#w = 1000
#h = 1000

#img_name = 'calanques'
#exp_name = 'collines'
#x = 6600
#y = 28800
#w = 1000
#h = 1000

#img_name = 'cannes'
#exp_name = 'theoule_sur_mer'
#x = 5100
#y = 32300
#w = 1000
#h = 1000

#img_name = 'mera'
#exp_name = 'crete'
#x = 15700
#y = 36400
#w = 1200
#h = 1000

img_name = 'uy1'
exp_name = 'campo'
# FULL ROI
#x = 4500
#y = 12000
#w = 8000
#h = 11000
# portion inside ROI
x = 5500
y = 25000
w = 1500
h = 1500

x = 7000
y = 25000
w = 2000
h = 2000



## Try to import the global parameters module
#  it permits to pass values between different modules
try:
   from python import global_params

   global_params.subsampling_factor=1
   global_params.subsampling_factor_registration=1

except ImportError:
  pass



im1 = 'pleiades_data/images/%s/im01.tif' % (img_name)
im2 = 'pleiades_data/images/%s/im02.tif' % (img_name)
im1_color = 'pleiades_data/images/%s/im01_color.tif' % (img_name)
rpc1 = 'pleiades_data/rpc/%s/rpc01.xml' % (img_name)
rpc2 = 'pleiades_data/rpc/%s/rpc02.xml' % (img_name)

rect1 = '/tmp/%s1.tif' % (exp_name)
rect2 = '/tmp/%s2.tif' % (exp_name)
hom1  = '/tmp/%s_hom1' % (exp_name)
hom2  = '/tmp/%s_hom2' % (exp_name)
rect1_color = '/tmp/%s1_color.tif' % (exp_name)
disp    = '/tmp/%s_disp.pgm'   % (exp_name)
mask    = '/tmp/%s_mask.png'   % (exp_name)
cloud   = '/tmp/%s_cloud.ply'  % (exp_name)
height  = '/tmp/%s_height.tif' % (exp_name)
rpc_err = '/tmp/%s_rpc_err.tif'% (exp_name)


def main():
    """
    Launches the s2p stereo pipeline on a pair of Pleiades images
    """

    # ATTENTION if subsampling_factor is set the rectified images will be smaller, 
    # and the homography matrices and disparity range will reflect this fact

    ## 1. rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1, rpc2,
        x, y, w, h, rect1, rect2)

    # save homographies to tmp files
    np.savetxt(hom1, H1)
    np.savetxt(hom2, H2)

    ## 2. block-matching
    #block_matching.compute_disparity_map(rect1, rect2, disp, mask,
    #    'hirschmuller02', disp_min, disp_max)
    #block_matching.compute_disparity_map(rect1, rect2, disp, mask,
    #    'hirschmuller02', disp_min, disp_max, '1 3')
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        'hirschmuller08', disp_min, disp_max)


    ## 3. triangulation
    triangulation.compute_height_map(rpc1, rpc2, hom1, hom2, disp, mask, height,
        rpc_err)

    ## 4. colorize and generate point cloud
    triangulation.colorize(rect1, im1_color, hom1, rect1_color)
    triangulation.compute_point_cloud(rect1_color, height, rpc1, hom1, cloud)

    # display results
    print "vflip %s %s %s %s %s %s" % (rect1, rect2, rect1_color, disp, mask, height)
    print "meshlab %s" % (cloud)

    ### cleanup 
    for name in common.garbage:
        common.run('rm ' + name)
    
if __name__ == '__main__': main()
