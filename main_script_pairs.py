#!/usr/bin/env python

import numpy as np
from python import rectification
from python import rpc_model
from python import common

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

# 1. rectification
H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1, rpc2,
                                                       x, y, w, h, rect1, rect2)
print disp_min, disp_max
#disp_min = -55
#disp_max = 0


# 2. block-matching
bm_binary = 'bin/stereoHirschmuller2002/subpix.sh'
disp = '/tmp/disp.pgm'
common.run("%s %s %s %s %d %d" %(bm_binary, rect1, rect2, disp, disp_min, disp_max))

# 3. triangulation
#np.savetxt('/tmp/h1', H1)
#np.savetxt('/tmp/h2', H2)
#common.run("disp_to_h %s /tmp/h1 %s /tmp/h2 %s %s.png %s /tmp/{l_out.ply, l_height.tif, l_RPCerr.tif}" % (rpc1, rpc2, disp, disp, im1))
