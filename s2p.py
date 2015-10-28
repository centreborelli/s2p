#!/usr/bin/env python

# s2p - Satellite Stereo Pipeline
# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@polytechnique.org>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import multiprocessing
import sys
import json
import numpy as np
import os.path
import copy
import glob
import time

from python import tee
from python import srtm
from python import common
from python import masking
from python import rpc_model
from python import rpc_utils
from python import geographiclib
from python import pointing_accuracy
from python import visualisation
from python import rectification
from python import block_matching
from python import triangulation
from python import tile_composer
from python import fusion
from python.config import cfg



def show_progress(a):
    show_progress.counter += 1
    if show_progress.counter > 1:
        print "Processed %d tiles" % show_progress.counter
    else:
        print "Processed 1 tile"



def getMinMaxFromExtract(tile_dir,tilesFullInfo):
    """
    Get min/max intensities of an extract ROI from the ref image.
    
    Args:
        fullInfo : all that you need to process a tile:
        tile_dir,pair_id,A_global,col,row,tw,th,ntx,nty,i,j,img1,rpc1,img2,rpc2
    """
    
    print "\nCrop ref image and compute min/max intensities..."

    #Get info
    col,row,tw,th,ov,i,j,images=tilesFullInfo[tile_dir]
    img1 = images[0]['img']
    
    # output files
    crop_ref = tile_dir + '/roi_ref.tif'
    local_minmax = tile_dir + '/local_minmax.txt'
    
    z = cfg['subsampling_factor']
    if z == 1:
        common.image_crop_TIFF(img1, col,row,tw,th, crop_ref)
    else:
        # gdal is used for the zoom because it handles BigTIFF files, and
        # before the zoom out the image may be that big
        tmp_crop = common.image_crop_TIFF(im1, col,row,tw,th)
        common.image_zoom_gdal(tmp_crop, z, crop_ref, tw,th)
		
    common.image_getminmax(crop_ref,local_minmax)
    
    
def colorCropRef(tile_dir,tilesFullInfo,clr=None):

    """
    Colorization of a crop_ref (for a given tile)
    
    Args:
        - fullInfo : all that you need to process a tile:
        tile_dir,pair_id,A_global,col,row,tw,th,ntx,nty,i,j,img1,rpc1,img2,rpc2
        
        - clr (optional) : if crop_ref is a pan image, will perform the pansharpening with the color image clr
            
            If clr is None then:
                case 1 : tile is an RGBI image, so get rid of I channel, and perform rescaling of the remaining channels
                case 2 : tile is already an RGB image, so just perform rescaling

                Note that if rescaling is already performed, then the file applied_minmax.txt exists :
                - if applied_minmax.txt exists and cfg['skip_existing'] is True, rescaling won't be performed again
                - if applied_minmax.txt exists and is different from global_minmax, rescaling will be compulsorily performed (can occur if a new tile is added)


    """

    # Get info
    col,row,tw,th,ov,i,j,images = tilesFullInfo[tile_dir]

    # Paths
    crop_ref = tile_dir + '/roi_ref.tif'
    global_minmax = cfg['out_dir']+'/global_minmax.txt'
    applied_minmax = tile_dir + '/applied_minmax.txt'
    
    global_minmax_arr = np.loadtxt(global_minmax)
            
    if cfg['color_ply']:
        
        crop_color = tile_dir + '/roi_color_ref.tif'
        if clr is not None:
            print 'colorizing...'
            triangulation.colorize(crop_ref, clr, col,row, z, crop_color)
        
        else: # use of image_rescaleintensities
        
            doProcess=False
            if not os.path.exists(applied_minmax):
                doProcess=True
                applied_minmax_arr = global_minmax_arr
            else:     
                applied_minmax_arr = np.loadtxt(applied_minmax)
		    
                if (applied_minmax_arr[0]!=global_minmax_arr[0]) or (applied_minmax_arr[1]!=global_minmax_arr[1]):
		            doProcess=True
		            applied_minmax_arr = global_minmax_arr
        
            if not doProcess and cfg['skip_existing']:
                print 'Rescaling of tile %s already done, skip' % tile_dir
            else:

                np.savetxt(applied_minmax,applied_minmax_arr)

                if common.image_pix_dim_tiffinfo(crop_ref) == 4:
                    print 'the image is pansharpened fusioned'
                    tmp = common.rgbi_to_rgb(crop_ref, out=None, tilewise=True)
                    #common.image_qauto(tmp, crop_color, tilewise=False)
                    common.image_rescaleintensities(tmp,crop_color,applied_minmax_arr[0],applied_minmax_arr[1])
                else:
                    print 'no color data'
                    #common.image_qauto(crop_ref, crop_color, tilewise=False)
                    common.image_rescaleintensities(crop_ref,crop_color,applied_minmax_arr[0],applied_minmax_arr[1])
	        
        

def generate_cloud(tile_dir,tilesFullInfo, do_offset=False):
    """
    Args:
        clr:  path to the xs (multispectral, ie color) reference image
        do_offset (optional, default: False): boolean flag to decide wether the
            x, y coordinates of points in the ply file will be translated or
            not (translated to be close to 0, to avoid precision loss due to
            huge numbers)
    """
    print "\nComputing point cloud..."
    
    # Get info
    col,row,tw,th,ov,i,j,images = tilesFullInfo[tile_dir]
    img1,rpc1 = images[0]['img'],images[0]['rpc']
    
    height_map = tile_dir + '/local_merged_height_map.tif'
    crop_color = tile_dir + '/roi_color_ref.tif'
    if not os.path.exists(crop_color):
        crop_color=''

    if cfg['full_img'] and z == 1:
        crop_ref = img1
    else: #Crop image has already been computed by getMinMaxFromExtract 
        crop_ref = tile_dir + '/roi_ref.tif'

    #Compute the homography transforming the coordinates system of the
    #original full size image into the coordinates system 
    #of the crop we are dealing with
    A = common.matrix_translation(-col, -row)
    z = cfg['subsampling_factor']
    f = 1.0/z
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = tile_dir + '/trans.txt'
    np.savetxt(trans, A)

    # compute coordinates (offsets) of the point we want to use as origin in the local coordinate system of the computed cloud
    if do_offset:
        r = rpc_model.RPCModel(rpc1)
        lat = r.latOff
        lon = r.lonOff
        off_x, off_y = geographiclib.geodetic_to_utm(lat, lon)[0:2]
    else:
        off_x, off_y = 0, 0



    # output 
    cloud = tile_dir + '/cloud.ply'
    
    triangulation.compute_point_cloud(cloud, height_map, rpc1, trans, crop_color,
                                      off_x, off_y)
                                                                            
    common.garbage_cleanup()



def generate_dsm(out, point_clouds_list, resolution):
    """
    Args:
        out: output geotiff file
        point_clouds_list: list of ply files
        resolution: in meters per pixel

    The point clouds are supposed to contain points in the same UTM zones.
    """
    if point_clouds_list:
        files = ' '.join(point_clouds_list)
        common.run("ls %s | plyflatten %f %s" % (files, resolution, out))


def crop_corresponding_areas(out_dir, images, roi, zoom=1):
    """
    Crops areas corresponding to the reference ROI in the secondary images.

    Args:
        out_dir:
        images: sequence of dicts containing the paths to input data
        roi: dictionary containing the ROI definition
        zoom: integer zoom out factor
    """
    rpc_ref = images[0]['rpc']
    for i, image in enumerate(images[1:]):
        x, y, w, h = rpc_utils.corresponding_roi(rpc_ref, image['rpc'],
                                                 roi['x'], roi['y'], roi['w'],
                                                 roi['h'])
        if zoom == 1:
            common.image_crop_TIFF(image['img'], x, y, w, h,
                                   '%s/roi_sec_%d.tif' % (out_dir, i))
        else:
            # gdal is used for the zoom because it handles BigTIFF files, and
            # before the zoom out the image may be that big
            tmp = common.image_crop_TIFF(image['img'], x, y, w, h)
            common.image_zoom_gdal(tmp, zoom, '%s/roi_sec_%d.tif' % (out_dir,
                                                                     i), w, h)


def check_parameters(usr_cfg):
    """
    Checks that the provided dictionary defines all the mandatory
    arguments, and warns about unknown optional arguments.

    Args:
        usr_cfg: python dict read from the json input file

    Returns:
        nothing
    """

    # verify that i/o files are defined
    if 'out_dir' not in usr_cfg:
        print "missing output dir: abort"
        sys.exit(1)
    if 'images' not in usr_cfg or len(usr_cfg['images']) < 2:
        print "missing input data paths: abort"
        sys.exit(1)
    for i in range(0,len(usr_cfg['images'])):
        if 'img' not in usr_cfg['images'][i] or 'rpc' not in usr_cfg['images'][i]:
            errMsg = 'missing input data paths for image ' + str(i) + ': abort'
            print errMsg
            sys.exit(1)           
	        

    # verify that roi or path to preview file are defined
    if ('full_img' not in usr_cfg) or (not usr_cfg['full_img']):
        if 'roi' not in usr_cfg or any(p not in usr_cfg['roi'] for p in ['x',
                                                                         'y',
                                                                         'w',
                                                                         'h']):
            if 'prv' not in usr_cfg['images'][0]:
                print """missing or incomplete roi definition, and no preview
                file is specified: abort"""
                sys.exit(1)

    # warn about unknown optional parameters: these parameters have no default
    # value in the global config.cfg dictionary, and thus they are not used
    # anywhere.  They may appear in the usr_cfg because of a typo.
    l = usr_cfg.keys()

    # remove mandatory parameters (they are not in config.cfg)
    l.remove('out_dir')
    l.remove('images')
    if 'roi' in l:
        l.remove('roi')

    # check
    for k in l:
        if k not in cfg:
            print """parameter %s unknown: you should remove it from the input
            json file. It will be ignored.""" % k





def pointing_correction(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None):
    """
    Computes pointing corrections, without tiling

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
            area contained in the full image.

    Returns:
        nothing
    """
    
    
    # create a directory for the experiment
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir


    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())

    # ensure that the coordinates of the ROI are multiples of the zoom factor
    z = cfg['subsampling_factor']
    x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)

    # check if the ROI is completely masked (water, or outside the image domain)
    H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
    if masking.cloud_water_image_domain(cwid_msk, w, h, H, rpc1, roi_msk,
                                        cld_msk):
        print "Tile masked by water or outside definition domain, skip"
        open("%s/this_tile_is_masked.txt" % out_dir, 'a').close()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        if not cfg['debug']:
            fout.close()
        return

    # correct pointing error
    # A is the correction matrix and m is the list of sift matches

    A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, x,
                                                    y, w, h)
    if A is not None:
        np.savetxt(pointing, A)
    if m is not None:
        np.savetxt(sift_matches, m)
        np.savetxt(center, np.mean(m[:, 2:4], 0))
        visualisation.plot_matches_pleiades(img1, img2, rpc1, rpc2, m, x, y,
                                                w, h, sift_matches_plot)
   
    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return


def rectify(out_dir, A_global, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None):
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
        prv1 (optional): path to a preview of the reference image
        cld_msk (optional): path to a gml file containing a cloud mask
        roi_msk (optional): path to a gml file containing a mask defining the
            area contained in the full image.

    Returns:
        nothing
    """
    
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir


    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())
                                                
    
    A,m=None,None

    if os.path.isfile('%s/pointing.txt' % out_dir): 
        A=np.loadtxt('%s/pointing.txt' % out_dir)
    else:
        A=A_global
    if os.path.isfile('%s/sift_matches.txt' % out_dir):
        m=np.loadtxt('%s/sift_matches.txt' % out_dir)

    # rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(img1, img2, rpc1,
                                                            rpc2, x, y, w, h,
                                                            rect1, rect2, A, m)
    
    z = cfg['subsampling_factor']                           
    np.savetxt(subsampling, np.array([z]))                                                        
    np.savetxt(H_ref, H1)
    np.savetxt(H_sec, H2)
    np.savetxt(disp_min_max, np.array([disp_min, disp_max]))



    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return


def disparity(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None):
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
            area contained in the full image.

    Returns:
        nothing
    """
    
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir


    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())
                                                

    # disparity (block-matching)
    
    disp_min,disp_max=np.loadtxt(disp_min_max)
    
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
    masking.cloud_water_image_domain(cwid_msk, ww, hh, H1, rpc1, roi_msk,
                                     cld_msk)
    try:
        masking.intersection(mask, mask, cwid_msk)
        masking.erosion(mask, mask, cfg['msk_erosion'])
    except OSError:
        print "file %s not produced" % mask


    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return




def triangulate(out_dir, img1, rpc1, img2, rpc2, x=None, y=None,
                             w=None, h=None, prv1=None, cld_msk=None,
                             roi_msk=None, A=None):
    """
    Computes triangulations, without tiling

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
            area contained in the full image.
        A (optional, default None): pointing correction matrix. 

    Returns:
        nothing
    """
    
    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    cwid_msk = '%s/cloud_water_image_domain_mask.png' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing.txt' % out_dir
    center = '%s/center_keypts_sec.txt' % out_dir
    sift_matches = '%s/sift_matches.txt' % out_dir
    sift_matches_plot = '%s/sift_matches_plot.png' % out_dir
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
    disp_min_max = '%s/disp_min_max.txt' % out_dir
    config = '%s/config.json' % out_dir
    rpc_err = '%s/rpc_err.tif' % out_dir
    height_map = '%s/height_map.tif' % out_dir

    # redirect stdout and stderr to log file
    if not cfg['debug']:
        fout = open('%s/stdout.log' % out_dir, 'w', 0)  # '0' for no buffering
        sys.stdout = fout
        sys.stderr = fout

    # debug print
    print 'tile %d %d running on process %s' % (x, y,
                                                multiprocessing.current_process())

    
    

    # triangulation
    z = cfg['subsampling_factor']
    triangulation.compute_dem(height_map, x, y, w, h, z,
                                              rpc1, rpc2, H_ref, H_sec, disp, mask,
                                              rpc_err, A)
                                                            

    # close logs
    common.garbage_cleanup()
    if not cfg['debug']:
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        fout.close()

    return

					 


def prepare_fullProcess(out_dir, images, x, y, w, h, tw=None, th=None,
                 ov=None, cld_msk=None, roi_msk=None):
              

    # create a directory for the experiment
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # duplicate stdout and stderr to log file
    tee.Tee('%s/stdout.log' % out_dir, 'w')

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        if prv1:
            x, y, w, h = common.get_roi_coordinates(img1, prv1)
        else:
            print 'Neither a ROI nor a preview file are defined. Aborting.'
            return

    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    z = cfg['subsampling_factor']
    x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)
    

    # TODO: automatically compute optimal size for tiles
    if tw is None and th is None and ov is None:
        ov = z * 100
        if w <= z * cfg['tile_size']:
            tw = w
        else:
            tw = z * cfg['tile_size']
        if h <= z * cfg['tile_size']:
            th = h
        else:
            th = z * cfg['tile_size']
    ntx = np.ceil(float(w - ov) / (tw - ov))
    nty = np.ceil(float(h - ov) / (th - ov))
    nt = ntx * nty

    print 'tiles size: (%d, %d)' % (tw, th)
    print 'total number of tiles: %d (%d x %d)' % (nt, ntx, nty)
    NbPairs = len(images)-1
    print 'total number of pairs: %d ' % NbPairs


    # Build tiles dicos
    tilesFullInfo={}
    pairedTilesPerPairId = {}
    for pair_id in range(1,len(images)) :
        pairedTilesPerPairId[pair_id]=[]
        for i, row in enumerate(np.arange(y, y + h - ov, th - ov)):
            for j, col in enumerate(np.arange(x, x + w - ov, tw - ov)):
                # ensure that the coordinates of the ROI are multiples of the zoom factor
                col, row, tw, th = common.round_roi_to_nearest_multiple(z, col, row, tw, th)
                tile_dir = '%s/tile_%d_%d_row_%d/col_%d/' % (out_dir, tw, th, row, col)
                paired_tile_dir = '%s/tile_%d_%d_row_%d/col_%d/pair_%d' % (out_dir, tw, th, row, col, pair_id)
                tilesFullInfo[tile_dir]=[col,row,tw,th,ov,i,j,images]
                pairedTilesPerPairId[pair_id].append(paired_tile_dir )



    return tilesFullInfo,pairedTilesPerPairId,NbPairs



def preprocess_tiles(out_dir,tilesFullInfo,pairedTilesPerPairId,NbPairs,cld_msk=None, roi_msk=None):

    # create pool with less workers than available cores
    nb_workers = multiprocessing.cpu_count()
    if cfg['max_nb_threads']:
        nb_workers = min(nb_workers, cfg['max_nb_threads'])
    pool = multiprocessing.Pool(nb_workers)
    
    # 1 - Pointing correction  
    results = []
    show_progress.counter = 0
    try:
        print 'Computing pointing correction...'
        for tile_dir in tilesFullInfo:
            for i in range(0,NbPairs):
            
                pair_id = i+1
                paired_tile_dir=tile_dir + 'pair_%d' % pair_id
                
                #tile_dir,pair_id,A_global,
                col,row,tw,th,ov,i,j,images=tilesFullInfo[tile_dir]
                img1,rpc1,img2,rpc2 = images[0]['img'],images[0]['rpc'],images[pair_id]['img'],images[pair_id]['rpc']


                # check if the tile is already done, or masked
                if os.path.isfile('%s/pointing.txt' % paired_tile_dir):
                    if cfg['skip_existing']:
                        print "pointing correction on tile %d %d (pair %d) already done, skip" % (col,row,pair_id)
                        #tiles.append(tile_dir)
                        continue
                    if os.path.isfile('%s/this_tile_is_masked.txt' % tile_dir):
                        print "tile %s already masked, skip" % tile_dir
                        #tiles.append(tile_dir)
                        continue
    
                
                if cfg['debug']:
                    pointing_correction(paired_tile_dir, img1, rpc1, img2, rpc2,
                                                 col, row, tw, th, None, cld_msk,
                                                 roi_msk)
                else:
                    p = pool.apply_async(pointing_correction,
                                             args=(paired_tile_dir, img1, rpc1, img2, rpc2,
                                                   col, row, tw, th, None, cld_msk,
                                                   roi_msk), callback=show_progress)
                    results.append(p)

        for r in results:
            try:
                r.get(3600)  # wait at most one hour per tile
            except multiprocessing.TimeoutError:
                print "Timeout while computing tile "+str(r)

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
        
    # 2 - Global pointing correction
    try:
        print 'Computing global pointing correction...'
        
        for i in range(0,NbPairs):
            
                pair_id = i+1
                tiles = pairedTilesPerPairId[pair_id]
                A_globalMat = pointing_accuracy.global_from_local(tiles)
                np.savetxt('%s/global_pointing_pair_%d.txt' % (out_dir,pair_id), A_globalMat)
            

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
        
    # 3 - Crop ref image into tile, and get min/max intensities values  
    results = []
    show_progress.counter = 0
    try:
        print 'Croping ref image into tiles, get min/max intensities values from each tile...'
        for tile_dir in tilesFullInfo:

			# check if the tile is already done, or masked
			if os.path.isfile('%s/local_minmax.txt' % tile_dir):
				if cfg['skip_existing']:
					print "Compute min/max intensities values on tile %s already done, skip" % (tile_dir)
					#tiles.append(tile_dir)
					continue
				if os.path.isfile('%s/this_tile_is_masked.txt' % (tile_dir)):
					print "tile %s already masked, skip" % (tile_dir)
					#tiles.append(tile_dir)
					continue

			
			if cfg['debug']:
				getMinMaxFromExtract(tile_dir,tilesFullInfo)
			else:
				p = pool.apply_async(getMinMaxFromExtract,
										 args=(tile_dir,tilesFullInfo), callback=show_progress)
				results.append(p)

        for r in results:
            try:
                r.get(3600)  # wait at most one hour per tile
            except multiprocessing.TimeoutError:
                print "Timeout while computing tile "+str(r)

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
        
    # 4 - Global Min/max
    try:
        print 'Computing global min max intensities...'
        
        minlist=[]
        maxlist=[]
        for tile_dir in tilesFullInfo:
            minmax = np.loadtxt(tile_dir + '/local_minmax.txt')
            minlist.append(minmax[0])
            maxlist.append(minmax[1])
            
        global_minmax=[min(minlist),max(maxlist)]
        
        np.savetxt(cfg['out_dir']+'/global_minmax.txt',global_minmax)			        

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
        
       


def mergeHeightMaps(height_maps,tile_dir,thresh,conservative,k=1,garbage=[]):

    local_merged_height_map = tile_dir +'/local_merged_height_map.tif'
    if os.path.isfile(local_merged_height_map) and cfg['skip_existing']:
        print 'final height map %s already done, skip' % local_merged_height_map
    else:
        list_height_maps=[]
        for i in range(len(height_maps)-1):
            height_map = tile_dir +'/height_map_'+str(i)+'_'+str(i+1)+'_'+str(k)+'.tif'
            fusion.merge(height_maps[i], height_maps[i+1], thresh, height_map,
                     conservative)
            list_height_maps.append(height_map)
            garbage.append(height_map)
        
        if len(list_height_maps) > 1:
            mergeHeightMaps(list_height_maps,tile_dir,thresh,conservative,k+1,garbage)
        else:
            common.run('cp %s %s' % (list_height_maps[0],local_merged_height_map))
            for imtemp in garbage:
                common.run('rm -f %s' % imtemp )


def process_tile(tile_dir,NbPairs,tilesFullInfo,cld_msk=None, roi_msk=None):

    height_maps = []

    #Get all info
    fullInfo=tilesFullInfo[tile_dir]
    col,row,tw,th,ov,i,j,images=fullInfo
    img1,rpc1 = images[0]['img'],images[0]['rpc']

    for i in range(0,NbPairs):
            
        pair_id = i+1
        paired_tile_dir=tile_dir + 'pair_%d' % pair_id

        img2,rpc2 = images[pair_id]['img'],images[pair_id]['rpc']
        A_global = '%s/global_pointing_pair_%d.txt' % (cfg['out_dir'],pair_id)

        if os.path.isfile('%s/this_tile_is_masked.txt' % tile_dir):
            print "tile %s already masked, skip" % paired_tile_dir

        else:
            print 'process tile %d %d...' % (col,row)
           
            height_maps.append('%s/height_map.tif' % paired_tile_dir)

            # Rectification
            if os.path.isfile('%s/rectified_ref.tif' % paired_tile_dir) and os.path.isfile('%s/rectified_sec.tif' % paired_tile_dir) and cfg['skip_existing']:
                print 'rectification on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
            else:
                print 'rectification : process tile %d %d (pair %d)...' % (col,row,pair_id)
                rectify(paired_tile_dir, A_global, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk)
            
            # Disparity
            if os.path.isfile('%s/rectified_disp.tif' % paired_tile_dir) and cfg['skip_existing']:
                print 'disparity on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
            else:
                print 'disparity : process tile %d %d (pair %d)...' % (col,row,pair_id)
                disparity(paired_tile_dir, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk)
                
            # Triangulate
            if os.path.isfile('%s/height_map.tif' % paired_tile_dir) and cfg['skip_existing']:
                print 'Triangulation on tile %d %d (pair %d) already done, skip' % (col, row, pair_id)
            else:
                print 'Triangulation : process tile %d %d (pair %d)...' % (col,row,pair_id)    
                triangulate(paired_tile_dir, img1, rpc1, img2, rpc2, col, row, tw, th, None, cld_msk, roi_msk, np.loadtxt(A_global))


    ## Finalization ##
                
    # merge the n height maps
    mergeHeightMaps(height_maps,tile_dir,cfg['fusion_thresh'],cfg['fusion_conservative'],1,[])    
    
    #Colors 
    colorCropRef(tile_dir,tilesFullInfo,None)
    
    #Generate cloud 
    generate_cloud(tile_dir,tilesFullInfo,cfg['offset_ply'])
    

def process_tiles(out_dir,tilesFullInfo,tilesLocPerPairId,NbPairs,cld_msk=None, roi_msk=None):
    
    # create pool with less workers than available cores
    nb_workers = multiprocessing.cpu_count()
    if cfg['max_nb_threads']:
        nb_workers = min(nb_workers, cfg['max_nb_threads'])
    pool = multiprocessing.Pool(nb_workers)
    
    
    #process_tiles
    results = []
    show_progress.counter = 0
    try:
        for tile_dir in tilesFullInfo:

            #fullInfoNpairs=[]
            #for i in range(0,NbPairs):
            #    pair_id = i+1
            #    tile_dir=tile + 'pair_%d' % pair_id
            #    fullInfo=tilesFullInfo[tile_dir]
            #    fullInfoNpairs.append(fullInfo)
            

            if cfg['debug']:
                process_tile(tile_dir,NbPairs,tilesFullInfo,cld_msk,roi_msk)
            else:
                p = pool.apply_async(process_tile,
                                         args=(tile_dir,NbPairs,tilesFullInfo, cld_msk,
                                               roi_msk), callback=show_progress)
                results.append(p)

        for r in results:
            try:
                r.get(3600)  # wait at most one hour per tile
            except multiprocessing.TimeoutError:
                print "Timeout while computing tile "+str(r)

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]
    


    
def tile_composition(height_map_path, tiles, x, y, w, h, tw=None, th=None, ov=None):


    if os.path.isfile(height_map_path) and cfg['skip_existing']:
        print 'tile composition for %s already done, skip' % height_map_path
    else:
        # ensure that the coordinates of the ROI are multiples of the zoom factor,
        # to avoid bad registration of tiles due to rounding problems.
        z = cfg['subsampling_factor']
        x, y, w, h = common.round_roi_to_nearest_multiple(z, x, y, w, h)
    
        # TODO: automatically compute optimal size for tiles
        if tw is None and th is None and ov is None:
            ov = z * 100
            if w <= z * cfg['tile_size']:
                tw = w
            else:
                tw = z * cfg['tile_size']
            if h <= z * cfg['tile_size']:
                th = h
            else:
                th = z * cfg['tile_size']
        
        tmp = ['%s/height_map.tif' % t for t in tiles]
        print "Mosaicing tiles with %s..." % cfg['mosaic_method']
        if cfg['mosaic_method'] == 'gdal':
            tile_composer.mosaic_gdal(height_map_path, w/z, h/z, tmp, tw/z, th/z, ov/z)
        else:
            tile_composer.mosaic(height_map_path, w/z, h/z, tmp, tw/z, th/z, ov/z)


def tiles_composition(out_dir,ensTiles,tilesLocPerPairId,NbPairs, x, y, w, h, tw=None, th=None, ov=None):

    # create pool with less workers than available cores
    nb_workers = multiprocessing.cpu_count()
    if cfg['max_nb_threads']:
        nb_workers = min(nb_workers, cfg['max_nb_threads'])
    pool = multiprocessing.Pool(nb_workers)

    # Tiles composition
    results = []
    show_progress.counter = 0
    out=[]
    try:

        for i in range(0,NbPairs):
        
            pair_id = i+1 
            height_map_path = '%s/height_map_pair_%d.tif' % (out_dir,pair_id)
            out.append(height_map_path)

            if cfg['debug']:
                    tile_composition(height_map_path,tilesLocPerPairId[pair_id],x, y, w, h, tw, th, ov)
            else:
                    p = pool.apply_async(tile_composition,
                                             args=(height_map_path,tilesLocPerPairId[pair_id],x, y, w, h, tw, th, ov), callback=show_progress)
                    results.append(p)

            for r in results:
                try:
                    r.get(3600)  # wait at most one hour per tile
                except multiprocessing.TimeoutError:
                    print "Timeout while computing tile "+str(r)

    except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)

    except common.RunFailure as e:
        print "FAILED call: ", e.args[0]["command"]
        print "\toutput: ", e.args[0]["output"]    
    
    return out
    
    

def main(config_file):
    """
    Launches s2p with the parameters given by a json file.

    Args:
        config_file: path to the config json file
    """
    # read the json configuration file
    f = open(config_file)
    user_cfg = json.load(f)
    f.close()

    # Check that all the mandatory arguments are defined, and warn about
    # 'unknown' params
    check_parameters(user_cfg)

    # fill the config module: updates the content of the config.cfg dictionary
    # with the content of the user_cfg dictionary
    cfg.update(user_cfg)

    # sets keys 'clr', 'cld' and 'roi' of the reference image to None if they
    # are not already defined. The default values of these optional arguments
    # can not be defined directly in the config.py module. They would be
    # overwritten by the previous update, because they are in a nested dict.
    cfg['images'][0].setdefault('clr')
    cfg['images'][0].setdefault('cld')
    cfg['images'][0].setdefault('roi')

    # update roi definition if the full_img flag is set to true
    if ('full_img' in cfg) and cfg['full_img']:
        sz = common.image_size_tiffinfo(cfg['images'][0]['img'])
        cfg['roi'] = {}
        cfg['roi']['x'] = 0
        cfg['roi']['y'] = 0
        cfg['roi']['w'] = sz[0]
        cfg['roi']['h'] = sz[1]

    # check that the roi is well defined
    if 'roi' not in cfg or any(p not in cfg['roi'] for p in ['x', 'y', 'w',
                                                             'h']):
        print "missing or incomplete ROI definition"
        print "ROI will be redefined by interactive selection"
        x, y, w, h = common.get_roi_coordinates(cfg['images'][0]['img'],
                                                cfg['images'][0]['prv'])
        cfg['roi'] = {}
        cfg['roi']['x'] = x
        cfg['roi']['y'] = y
        cfg['roi']['w'] = w
        cfg['roi']['h'] = h

    # check the zoom factor
    z = cfg['subsampling_factor']
    assert(z > 0 and z == np.floor(z))

    # create tmp dir and output directory for the experiment, and store a json
    # dump of the config.cfg dictionary there
    if not os.path.exists(cfg['temporary_dir']):
        os.makedirs(cfg['temporary_dir'])
    if not os.path.exists(os.path.join(cfg['temporary_dir'], 'meta')):
        os.makedirs(os.path.join(cfg['temporary_dir'], 'meta'))
    if not os.path.exists(cfg['out_dir']):
        os.makedirs(cfg['out_dir'])
    f = open('%s/config.json' % cfg['out_dir'], 'w')
    json.dump(cfg, f, indent=2)
    f.close()

    # measure total runtime
    t0 = time.time()

    # needed srtm tiles
    srtm_tiles = srtm.list_srtm_tiles(cfg['images'][0]['rpc'],
                                           *cfg['roi'].values())
    for s in srtm_tiles:
        srtm.get_srtm_tile(s, cfg['srtm_dir'])


    tilesFullInfo,pairedTilesPerPairId,NbPairs = prepare_fullProcess(cfg['out_dir'], cfg['images'], cfg['roi']['x'],
                           cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h'],
                           None, None, None, cfg['images'][0]['cld'],
                           cfg['images'][0]['roi'])
                           
                    
    preprocess_tiles(cfg['out_dir'],tilesFullInfo,pairedTilesPerPairId,NbPairs,
                        cfg['images'][0]['cld'], cfg['images'][0]['roi'])         
     
     
    process_tiles(cfg['out_dir'],tilesFullInfo,pairedTilesPerPairId,NbPairs,
                        cfg['images'][0]['cld'], cfg['images'][0]['roi'])                         


    #TODO
    # digital surface model
    #for i,height_map in enumerate(height_maps):
     #   pair_id = i+1
      #  out_dsm = '%s/dsm_pair_%d.tif' % (cfg['out_dir'] , pair_id)
       # cloud   = '%s/cloud_pair_%d.ply' % ( cfg['out_dir'] , pair_id )
        #common.run("ls %s | plyflatten %f %s" % (cloud, cfg['dsm_resolution'], out_dsm))
        

    # crop corresponding areas in the secondary images
    if not cfg['full_img']:
        crop_corresponding_areas(cfg['out_dir'], cfg['images'], cfg['roi'])


    # also copy the RPC's
    for i in range(len(cfg['images'])):
        from shutil import copy2
        copy2(cfg['images'][i]['rpc'], cfg['out_dir'])


    # runtime
    t = int(time.time() - t0)
    h = t/3600
    m = (t/60) % 60
    s = t % 60
    print "Total runtime: %dh:%dm:%ds" % (h, m, s)
    common.garbage_cleanup()


if __name__ == '__main__':

    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        print """
        Incorrect syntax, use:
          > %s config.json

          Launches the s2p pipeline. All the parameters, paths to input and
          output files, are defined in the json configuration file.
        """ % sys.argv[0]
        sys.exit(1)
