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


def color_crop_ref(tile_dir,tilesFullInfo,clr=None):

    """
    Colorizations of a crop_ref (for a given tile)
    
    Args:
        - tile_dir : a key for the dictionary tilesFullInfo; refers to a particular tile
        - tilesFullInfo : a dictionary that provides all you need to process a tile -> col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk
        
        - clr (optional) : if crop_ref is a pan image, will perform the pansharpening with the color image clr
            
            If clr is None then:
                case 1 : tile is an RGBI image, so removes I channel, and perform rescaling of the remaining channels
                case 2 : tile is already an RGB image, so just perform rescaling

                Note that if rescaling is already performed, then the file applied_minmax.txt exists :
                - if applied_minmax.txt exists and cfg['skip_existing'] is True, rescaling won't be performed again
                - if applied_minmax.txt exists and is different from global_minmax, rescaling will be compulsorily performed (can occur if a new tile is added)
    """

    # Get info
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    z = cfg['subsampling_factor']

    # Paths
    #crop_ref = tile_dir + '/roi_ref.tif'
    crop_ref = tile_dir + '/roi_ref_crop.tif'
    global_minmax = cfg['out_dir']+'/global_minmax.txt'
    applied_minmax = tile_dir + '/applied_minmax.txt'
    
    global_minmax_arr = np.loadtxt(global_minmax)
            
    if cfg['color_ply']:
        
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

            crop_color = tile_dir + '/roi_color_ref.tif'
            if clr is not None:
                triangulation.colorize(crop_ref, clr, col,row, z, crop_color, applied_minmax_arr[0],applied_minmax_arr[1])
            else: # use of image_rescaleintensities
            
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
        - tile_dir : a key for the dictionary tilesFullInfo; refers to a particular tile
        - tilesFullInfo : a dictionary that provides all you need to process a tile -> col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk
        - do_offset (optional, default: False): boolean flag to decide wether the
            x, y coordinates of points in the ply file will be translated or
            not (translated to be close to 0, to avoid precision loss due to
            huge numbers)
    """
    print "\nComputing point cloud..."
    
    # Get info
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    img1,rpc1 = images[0]['img'],images[0]['rpc']
    
    #height_map = tile_dir + '/local_merged_height_map.tif'
    height_map = tile_dir +'/local_merged_height_map_crop.tif'
    crop_color = tile_dir + '/roi_color_ref.tif'
    if not os.path.exists(crop_color):
        crop_color=''

    z = cfg['subsampling_factor']

    #Compute the homography transforming the coordinates system of the
    #original full size image into the coordinates system 
    #of the crop we are dealing with
    A = common.matrix_translation(-col*z, -row*z) # col and row have been divided by z inside 'finalize_tile' for convinience; col*z and row*z allow to get the initial values back.
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
    
    
def merge_height_maps(height_maps,tile_dir,thresh,conservative,k=1,garbage=[]):
    """
    Merges a list of height maps recursively, computed for one tile from N image pairs.
  
    Args :
         - height_maps : list of height map directories
         - tile_dir : directory of the tile from which to get a merged height map
         - thresh : threshold used for the fusion algorithm, in meters.
         - conservative (optional, default is False): if True, keep only the
            pixels where the two height map agree (fusion algorithm)
         - k : used to identify the current call of merge_height_maps (default = 1, first call)
         - garbage : a list used to remove temp data (default = [], first call)

    """

    #output file
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
            merge_height_maps(list_height_maps,tile_dir,thresh,conservative,k+1,garbage)
        else:
            common.run('cp %s %s' % (list_height_maps[0],local_merged_height_map))
            for imtemp in garbage:
                common.run('rm -f %s' % imtemp )    


def finalize_tile(tile_dir, height_maps, tilesFullInfo):
    """
    Finalizes the processing of a tile : merge the height maps from the N pairs, remove overlapping areas, 
    get the colors from a XS image and use it to color and generate a ply file (colorization is not mandatory)

    Args:
        - tile_dir : directory of the tile to be processed
        - height_maps : list of the height maps generated from N pairs
        - tilesFullInfo : a dictionary that provides all you need to process a tile for a given tile directory : col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk = tilesFullInfo[tile_dir]
    """

    # Get all info
    fullInfo=tilesFullInfo[tile_dir]
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk=fullInfo
    img1,rpc1 = images[0]['img'],images[0]['rpc']


    # merge the n height maps
    local_merged_height_map = tile_dir +'/local_merged_height_map.tif'
    if len(height_maps)>1:
        merge_height_maps(height_maps,tile_dir,cfg['fusion_thresh'],cfg['fusion_conservative'],1,[])
    else:    
        common.run('cp %s %s' % (height_maps[0],local_merged_height_map))

    # Remove overlapping areas
    #By tile
    local_merged_height_map = tile_dir +'/local_merged_height_map.tif'
    local_merged_height_map_crop = tile_dir +'/local_merged_height_map_crop.tif'
    crop_ref = tile_dir + '/roi_ref.tif'
    crop_ref_crop = tile_dir + '/roi_ref_crop.tif'


    dicoPos={}
    dicoPos['M']  = [ov/2,ov/2,-ov,-ov]
    dicoPos['L']  = [0,ov/2,-ov/2,-ov]
    dicoPos['R']  = [ov/2,ov/2,-ov/2,-ov]
    dicoPos['U']  = [ov/2,0,-ov,-ov/2]
    dicoPos['B']  = [ov/2,ov/2,-ov,-ov/2]
    dicoPos['UL'] = [0,0,-ov/2,-ov/2]
    dicoPos['UR'] = [ov/2,0,-ov/2,-ov/2]
    dicoPos['BR'] = [ov/2,ov/2,-ov/2,-ov/2]
    dicoPos['BL'] = [0,ov/2,-ov/2,-ov/2]
    dicoPos['Single'] = [0,0,0,0]  

    z = cfg['subsampling_factor']
    #tilesFullInfo[tile_dir] = dicoPos[pos][0]
    newcol,newrow,difftw,diffth = np.array(dicoPos[pos])/z
    info=tilesFullInfo[tile_dir]
    info[0] = info[0]/z + newcol
    info[1] = info[1]/z + newrow
    info[2] = info[2]/z + difftw
    info[3] = info[3]/z + diffth
    tilesFullInfo[tile_dir]=info   
   
    # z=1 beacause local_merged_height_map, crop_ref (and so forth) have already been zoomed. So don't zoom again to crop these images.
    common.cropImage(local_merged_height_map,local_merged_height_map_crop,newcol,newrow,info[2],info[3],1)
    common.cropImage(crop_ref,crop_ref_crop,newcol,newrow,info[2],info[3],1)   

        
    #By pair    
    for i in range(0,NbPairs):
        pair_id = i+1
            
        single_height_map = tile_dir + '/pair_%d/height_map.tif' % pair_id
        single_height_map_crop =tile_dir + '/pair_%d/height_map_crop.tif' % pair_id
          
        single_rpc_err = tile_dir + '/pair_%d/rpc_err.tif' % pair_id
        single_rpc_err_crop =tile_dir + '/pair_%d/rpc_err_crop.tif' % pair_id
            
        common.cropImage(single_height_map,single_height_map_crop,newcol,newrow,info[2],info[3],1)
        common.cropImage(single_rpc_err,single_rpc_err_crop,newcol,newrow,info[2],info[3],1)

    
    
    #if not cfg['debug']:
     #    common.run('rm -f %s' % (local_merged_height_map)) 
         #don't remove crop_ref !!
                   
    
    # Colors 
    color_crop_ref(tile_dir,tilesFullInfo,cfg['images'][0]['clr'])
    
    # Generate cloud 
    generate_cloud(tile_dir,tilesFullInfo,cfg['offset_ply'])
    
    
    
# ----------------------------------------------------------------------------------------------------
# ----------------------------------- rectify disparity triangulate ----------------------------------
# ----------------------------------------------------------------------------------------------------

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
    disp = '%s/rectified_disp.tif' % (out_dir)
    mask = '%s/rectified_mask.png' % (out_dir)
    H_ref = '%s/H_ref.txt' % out_dir
    H_sec = '%s/H_sec.txt' % out_dir
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

# ----------------------------------------------------------------------------------------------------
# ------------------------------- rectify disparity triangulate (end) --------------------------------
# ----------------------------------------------------------------------------------------------------		
