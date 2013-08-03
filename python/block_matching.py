#!/usr/bin/env python

import os
import numpy as np
from python import common


# define paths of various bm binaries
hirshmuller = '%s/../3rdparty/stereo_hirschmuller_2002/subpix.sh' % (
    os.path.dirname( __file__))


def compute_disparity_map(im1, im2, disp_range, out_disp, algo, extra_params=''):
    # read the disp_range file
    disp = np.loadtxt(disp_range)
    disp_min = disp[0]
    disp_max = disp[1]

    # call the block_matching binary
    if (algo == 'hirshmuller'):
        bm_binary = hirshmuller
        common.run("%s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            disp_min, disp_max, extra_params))

    elif (algo == 'msmw'):
        bm_binary = iip_stereo_correlation_multi_win2
        common.run("%s -i 1 -n 4 -p 4 -W 9 -x 5 -y 5 -r 1 -d 1 -t 1 -s 0 -b 0 -o 0.25 -f 0 -P 1 -m %d -M %d %s %s %s %s" %(bm_binary, 
            disp_min, disp_max,im1, im2, out_disp, out_mask ))
