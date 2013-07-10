#!/usr/bin/env python

import os
import numpy as np
from python import common


# define paths of various bm binaries
hirshmuller = '%s/../stereoHirschmuller2002/subpix.sh' % (os.path.dirname(
                                                                __file__))


def compute_disparity_map(im1, im2, disp_range, out_disp, algo):
    # read the disp_range file
    disp = np.loadtxt(disp_range)
    disp_min = disp[0]
    disp_max = disp[1]

    # call the block_matching binary
    if (algo == 'hirshmuller'):
        bm_binary = hirshmuller
        common.run("%s %s %s %s %d %d" %(bm_binary, im1, im2, out_disp,
            disp_min, disp_max))
