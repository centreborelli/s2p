#!/usr/bin/env python

import os
import numpy as np
from python import common


# define paths of various bm binaries
hirshmuller = '%s/../3rdparty/stereo_hirschmuller_2002/subpix.sh' % (
    os.path.dirname( __file__))

hirshmuller08 = '%s/../3rdparty/stereo_hirschmuller_2008/callSGBM.sh' % (
    os.path.dirname( __file__))


def compute_disparity_map(im1, im2, disp_range, out_disp, algo, extra_params=''):
    """
    Runs a block-matching binary on a pair of stereo-rectified images.

    Args:
        im1, im2: rectified stereo pair
        disp_range: path to a text file containing the disparity range
        out_disp: path to the output file (diparity map)
        algo: string used to indicate the desired binary
        extra_params: optional string with algorithm-dependent parameters
    """
    # read the disp_range file
    disp = np.loadtxt(disp_range)
    disp_min = disp[0]
    disp_max = disp[1]

    # call the block_matching binary
    if (algo == 'hirshmuller'):
        bm_binary = hirshmuller
        common.run("%s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            disp_min, disp_max, extra_params))
        # extra_params: LoG(0) regionRadius(3)
        #    LoG: Laplacian of Gaussian preprocess 1:enabled 0:disabled
        #    regionRadius: radius of the window

    if (algo == 'hirshmuller08'):
        bm_binary = hirshmuller08
        common.run("%s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            disp_min, disp_max, extra_params))
        # extra_params: regionRadius(3) P1(default) P2(default) LRdiff(1)
        #    regionRadius: radius of the window
        #    P1,P2 : regularization parameters 
        #    LRdiff: maximum difference between left and right disparity maps

    elif (algo == 'msmw'):
        bm_binary = iip_stereo_correlation_multi_win2
        common.run("%s -i 1 -n 4 -p 4 -W 9 -x 5 -y 5 -r 1 -d 1 -t 1 -s 0 -b 0 -o 0.25 -f 0 -P 1 -m %d -M %d %s %s %s %s" %(bm_binary, 
            disp_min, disp_max,im1, im2, out_disp, out_mask ))
