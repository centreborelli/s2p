#!/usr/bin/env python

import os
import numpy as np
from python import common


# define paths of various bm binaries
hirschmuller02 = '%s/../3rdparty/stereo_hirschmuller_2002/subpix.sh' % (
    os.path.dirname( __file__))

hirschmuller08 = '%s/../3rdparty/stereo_hirschmuller_2008/callSGBM.sh' % (
    os.path.dirname( __file__))

hirschmuller08_laplacian = '%s/../3rdparty/stereo_hirschmuller_2008/callSGBM_lap.sh' % (
    os.path.dirname( __file__))

msmw = '%s/../3rdparty/msmw/build/libstereo/iip_stereo_correlation_multi_win2' % (
    os.path.dirname( __file__))

def compute_disparity_map(im1, im2, out_disp, out_mask, algo, disp_min, disp_max, extra_params=''):
    """
    Runs a block-matching binary on a pair of stereo-rectified images.

    Args:
        im1, im2: rectified stereo pair
        out_disp: path to the output diparity map
        out_mask: path to the output rejection mask
        algo: string used to indicate the desired binary. Currently it can be
            one among 'hirschmuller02', 'hirschmuller08', 'hirschmuller08_laplacian' and 'msmw'
        disp_min : smallest disparity to consider
        disp_max : biggest disparity to consider
        extra_params: optional string with algorithm-dependent parameters
    """
    # call the block_matching binary
    if (algo == 'hirschmuller02'):
        bm_binary = hirschmuller02
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: LoG(0) regionRadius(3)
        #    LoG: Laplacian of Gaussian preprocess 1:enabled 0:disabled
        #    regionRadius: radius of the window

    if (algo == 'hirschmuller08'):
        bm_binary = hirschmuller08
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: regionRadius(3) P1(default) P2(default) LRdiff(1)
        #    regionRadius: radius of the window
        #    P1, P2 : regularization parameters
        #    LRdiff: maximum difference between left and right disparity maps

    if (algo == 'hirschmuller08_laplacian'):
        bm_binary = hirschmuller08
        common.run("%s %s %s %s %s %d %d %s" %(bm_binary, im1, im2, out_disp,
            out_mask, disp_min, disp_max, extra_params))
        # extra_params: regionRadius(3) P1(default) P2(default) LRdiff(1)
        #    regionRadius: radius of the window
        #    P1, P2 : regularization parameters
        #    LRdiff: maximum difference between left and right disparity maps

    elif (algo == 'msmw'):
        bm_binary = msmw
        common.run("%s -i 1 -n 4 -p 4 -W 5 -x 9 -y 9 -r 1 -d 1 -t -1 -s 0 -b 0 -o 0.25 -f 0 -P 32 -m %d -M %d %s %s %s %s" %(bm_binary,
            disp_min, disp_max, im1, im2, out_disp, out_mask))
