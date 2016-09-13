# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>


import os
import numpy as np

from config import cfg
from python import common
from python import block_matching
from python import rectification
from python import masking


def cargarse_basura(inputf, outputf):
    se=5
    tmp1 = outputf + '1.tif'
    tmp2 = outputf + '2.tif'
    tmpM = outputf + 'M.tif'
    common.run('morphoop %s min %d %s'%(inputf,se,tmpM))
    common.run('morphoop %s max %d %s'%(inputf,se,tmp1))
    common.run('morphoop %s max %d %s'%(inputf,se,tmpM))
    common.run('morphoop %s min %d %s'%(inputf,se,tmp2))
    common.run('plambda %s %s %s "x y - fabs %d > nan z if" -o %s'%(tmp1,tmp2,inputf,5,tmpM))
    common.run('remove_small_cc %s %s %d %d'%(tmpM,outputf,200,5))
    common.run('rm -f %s %s %s'%(tmp1,tmp2,tmpM))


def rectify(out_dir, A_global, img1, rpc1, img2, rpc2, x, y, w, h, margin=0):
    """
    Perform stereo rectification of a given ROI in a pair of images.

    Args:
        out_dir: path to the output directory
        A_global:
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        margin (optional): horizontal margin added on the left an right sides of
            the rectified images
    """
    try:
        A = np.loadtxt(os.path.join(out_dir, 'pointing.txt'))
    except IOError:
        A = A_global
    try:
        m = np.loadtxt(os.path.join(out_dir, 'sift_matches.txt'))
    except IOError:
        m = None

    # rectification
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m,
                                                            margin=margin)
    np.savetxt(os.path.join(out_dir, 'H_ref.txt'), H1, fmt='%12.6f')
    np.savetxt(os.path.join(out_dir, 'H_sec.txt'), H2, fmt='%12.6f')
    np.savetxt(os.path.join(out_dir, 'disp_min_max.txt'), [disp_min, disp_max],
                            fmt='%3.1f')


def disparity(out_dir):
    """
    Computes a disparity map from a Pair of Pleiades images, without tiling

    Args:
        out_dir: path to the output directory
    """
    # disparity (block-matching)
    rect1 = os.path.join(out_dir, 'rectified_ref.tif')
    rect2 = os.path.join(out_dir, 'rectified_sec.tif')
    disp = os.path.join(out_dir, 'rectified_disp.tif')
    mask = os.path.join(out_dir, 'rectified_mask.png')
    disp_min, disp_max = np.loadtxt(os.path.join(out_dir, 'disp_min_max.txt'))
    if cfg['disp_min'] is not None: disp_min = cfg['disp_min']
    if cfg['disp_max'] is not None: disp_max = cfg['disp_max']
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # add margin around masked pixels
    masking.erosion(mask, mask, cfg['msk_erosion'])
