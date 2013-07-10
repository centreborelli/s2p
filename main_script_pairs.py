#!/usr/bin/env python

import numpy as np
from python import common
from python import rectification
from python import block_matching
from python import triangulation

img_name = 'toulouse'
exp_name = 'blagnac'

im1  = 'pleiades_data/images/%s/im01.tif' % (img_name)
im2  = 'pleiades_data/images/%s/im02.tif' % (img_name)
rpc1 = 'pleiades_data/rpc/%s/rpc01.xml' % (img_name)
rpc2 = 'pleiades_data/rpc/%s/rpc02.xml' % (img_name)
x = 6000
y = 6000
w = 1000
h = 1000
rect1 = '/tmp/%s1.tif' % (exp_name)
rect2 = '/tmp/%s2.tif' % (exp_name)
disp_range  = '/tmp/disp_range'
disp  = '/tmp/disp.pgm'
mask  = '/tmp/disp.pgm.mask.png'
cloud = '/tmp/cloud.ply'


## 1. rectification
H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1, rpc2,
                                                       x, y, w, h, rect1, rect2)
print disp_min, disp_max

# save homographies and disp_range to tmp files
np.savetxt('/tmp/h1', H1)
np.savetxt('/tmp/h2', H2)
np.savetxt('/tmp/disp_range', [disp_min, disp_max])


## 2. block-matching
block_matching(rect1, rect2, disp_range, disp, 'hirshmuller')


## 3. triangulation
triangulation(rpc1, '/tmp/h1', rpc2, '/tmp/h2', disp, mask, im1, cloud)
