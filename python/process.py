# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import numpy as np
from config import cfg
import os
import sys
import multiprocessing

from python import common
from python import geographiclib
from python import triangulation
from python import fusion
from python import block_matching
from python import rectification
from python import masking


def color_crop_ref(tile_info, clr=None):
    """
    Colorizations of a crop_ref (for a given tile)

    Args:
        tile_info: a dictionary that provides all you need to process a tile
        clr (optional): if crop_ref is a pan image, will perform the pansharpening with the color image clr

        If clr is None then:
            case 1: tile is an RGBI image, so removes I channel, and perform rescaling of the remaining channels
            case 2: tile is already an RGB image, so just perform rescaling

        Note that if rescaling is already performed, then the file applied_minmax.txt exists:
            if applied_minmax.txt exists and cfg['skip_existing'] is True, rescaling won't be performed again
            if applied_minmax.txt exists and is different from global_minmax, rescaling will be compulsorily performed (can occur if a new tile is added)
    """
    # get info
    x, y, w, h = tile_info['coordinates']
    tile_dir = tile_info['directory']
    z = cfg['subsampling_factor']

    # paths
    global_minmax = cfg['out_dir'] + '/global_minmax.txt'
    applied_minmax = tile_dir + '/applied_minmax.txt'

    global_minmax_arr = np.loadtxt(global_minmax)

    # crop ref
    crop_ref = os.path.join(tile_dir, 'roi_ref_crop.tif')
    common.image_crop_tif(cfg['images'][0]['img'], x, y, w, h, crop_ref)

    if cfg['color_ply']:

        doProcess = False
        if not os.path.exists(applied_minmax):
            doProcess = True
            applied_minmax_arr = global_minmax_arr
        else:
            applied_minmax_arr = np.loadtxt(applied_minmax)

            if (applied_minmax_arr[0] != global_minmax_arr[0]) or (applied_minmax_arr[1] != global_minmax_arr[1]):
                doProcess = True
                applied_minmax_arr = global_minmax_arr

        if not doProcess and cfg['skip_existing']:
            print 'Rescaling of tile %s already done, skip' % tile_dir
        else:

            crop_color = tile_dir + '/roi_color_ref.tif'
            if clr is not None:
                triangulation.colorize(crop_ref, clr, x, y, z, crop_color,
                                       applied_minmax_arr[0],
                                       applied_minmax_arr[1])
            else:  # use of image_rescaleintensities
                np.savetxt(applied_minmax, applied_minmax_arr)

                if common.image_pix_dim_tiffinfo(crop_ref) == 4:
                    print 'the image is pansharpened fusioned'
                    tmp = common.rgbi_to_rgb(crop_ref, out=None, tilewise=True)
                    #common.image_qauto(tmp, crop_color, tilewise=False)
                    common.image_rescaleintensities(tmp, crop_color,
                                                    applied_minmax_arr[0],
                                                    applied_minmax_arr[1])
                else:
                    print 'no color data'
                    #common.image_qauto(crop_ref, crop_color, tilewise=False)
                    common.image_rescaleintensities(crop_ref, crop_color,
                                                    applied_minmax_arr[0],
                                                    applied_minmax_arr[1])


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


def rectify(out_dir, A_global, img1, rpc1, img2, rpc2, x=None, y=None,
            w=None, h=None, prv1=None):
    """
    Computes rectifications, without tiling

    Args:
        out_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image.

    Returns:
        nothing
    """
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir

    A, m = None, None

    if os.path.isfile('%s/pointing.txt' % out_dir):
        A = np.loadtxt('%s/pointing.txt' % out_dir)
    else:
        A = A_global
    if os.path.isfile('%s/sift_matches.txt' % out_dir):
        m = np.loadtxt('%s/sift_matches.txt' % out_dir)

    # rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m)

    z = cfg['subsampling_factor']
    np.savetxt(subsampling, np.array([z]), fmt='%.1f')
    np.savetxt(H_ref, H1, fmt='%12.6f')
    np.savetxt(H_sec, H2, fmt='%12.6f')
    np.savetxt(disp_min_max, np.array([disp_min, disp_max]), fmt='%3.1f')


def disparity(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
              w=None, h=None, prv1=None):
    """
    Computes a disparity map from a Pair of Pleiades images, without tiling

    Args:
        out_dir: path to the output directory
        img1: path to the reference image.
        rpc1: paths to the xml file containing the rpc coefficients of the
            reference image
        img2: path to the secondary image.
        rpc2: paths to the xml file containing the rpc coefficients of the
            secondary image
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        prv1 (optional): path to a preview of the reference image
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image
        wat_msk (optional): path to a tiff file containing a water mask.

    Returns:
        nothing
    """
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    cwid_msk_rect = '%s/rectified_cloud_water_image_domain_mask.png' % (out_dir)

    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir

    # disparity (block-matching)
    disp_min, disp_max = np.loadtxt(disp_min_max)

    if cfg['disp_min'] is not None:
        disp_min = cfg['disp_min']
    if cfg['disp_max'] is not None:
        disp_max = cfg['disp_max']
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
                                         cfg['matching_algorithm'], disp_min,
                                         disp_max)

    # intersect mask with the cloud_water_image_domain mask (recomputed here to
    # get to be sampled on the epipolar grid)
    ww, hh = common.image_size(rect1)
    H1 = np.loadtxt(H_ref)
    H_inv = np.array([[1, 0, x], [0, 1, y], [0, 0, 1]])
    common.image_apply_homography(cwid_msk_rect, cwid_msk, np.dot(H1,H_inv), ww, hh)

    try:
        masking.intersection(mask, mask, cwid_msk_rect)
        masking.erosion(mask, mask, cfg['msk_erosion'])
    except OSError:
        print "file %s not produced" % mask
