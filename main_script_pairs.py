#!/usr/bin/env python

import numpy as np
from python import common
from python import rectification
from python import block_matching
from python import triangulation

img_name = 'toulouse'
exp_name = 'blagnac'
x = 6000
y = 6000
w = 1000
h = 1000

#img_name = 'calanques'
#exp_name = 'collines'
#x = 6600
#y = 28800
#w = 1000
#h = 1000

#img_name = 'calanques'
#exp_name = 'collines2'
#x = 28800
#y = 6600
#w = 1000
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
w = 1000
h = 1000


im1  = 'pleiades_data/images/%s/im01.tif' % (img_name)
im2  = 'pleiades_data/images/%s/im02.tif' % (img_name)
im1_color  = 'pleiades_data/images/%s/im01_color.tif' % (img_name)
rpc1 = 'pleiades_data/rpc/%s/rpc01.xml' % (img_name)
rpc2 = 'pleiades_data/rpc/%s/rpc02.xml' % (img_name)
rect1 = '/tmp/%s1.tif' % (exp_name)
rect2 = '/tmp/%s2.tif' % (exp_name)
rect1_color = '/tmp/%s1_color.tif' % (exp_name)
disp_range  = '/tmp/%s_disp_range'% (exp_name)
disp  = '/tmp/%s_disp.pgm'% (exp_name)
mask  = '/tmp/%s_disp.pgm.mask.png'% (exp_name)
cloud = '/tmp/%s_cloud.ply'% (exp_name)
height = '/tmp/%s_height.tif'% (exp_name)
rpc_err = '/tmp/%s_rpc_err.tif'% (exp_name)


## 1. rectification
H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1, rpc2,
    x, y, w, h, rect1, rect2)

# save homographies and disp_range to tmp files
np.savetxt('/tmp/h1', H1)
np.savetxt('/tmp/h2', H2)
np.savetxt(disp_range, [disp_min, disp_max])


## 2. block-matching
# TODO: ATTENTION 
#   mask = '/tmp/%s_disp.pgm.mask.png is hardcoded for the hirshmuller method, 
#   the compute_disparity_map function shoould receive it like:
#   block_matching.compute_disparity_map(rect1, rect2, disp_range, disp, mask,  'hirshmuller')

#block_matching.compute_disparity_map(rect1, rect2, disp_range, disp,
#    'hirshmuller')
#block_matching.compute_disparity_map(rect1, rect2, disp_range, disp,
#    'hirshmuller', '1 3')
block_matching.compute_disparity_map(rect1, rect2, disp_range, disp,
    'hirshmuller08')


## 3. triangulation
triangulation.compute_height_map(rpc1, rpc2, '/tmp/h1', '/tmp/h2', disp, mask,
    height, rpc_err)

## 4. colorize and generate point cloud
triangulation.colorize(rect1, im1_color, '/tmp/h1', rect1_color)
triangulation.compute_point_cloud(rect1_color, height, rpc1, '/tmp/h1', cloud)

# display results
print "vflip %s %s %s %s %s %s" % (rect1, rect2, rect1_color, disp, mask, height)
print "meshlab %s" % (cloud)
