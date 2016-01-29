# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import multiprocessing
import numpy as np

from config import cfg
from python import common
from python import masking
from python import pointing_accuracy
from python import visualisation



def getMinMaxFromExtract(tile_info):
    """
    Gets min/max intensities of an extract ROI from the ref image.
    
    Args:
        tile_info : a list that provides all you need to process a tile: col,
            row, tw, th, ov, i, j, pos, x, y, w, h, images, NbPairs, cld_msk,
            roi_msk
    """
    
    print "\nCrop ref image and compute min/max intensities..."

    #Get info
    col,row,tw,th,ov,i,j,pos,x,y,w,h,images,NbPairs,cld_msk,roi_msk, tile_dir = tile_info
    img1 = images[0]['img']
    
    # output files
    crop_ref = tile_dir + '/roi_ref.tif'
    local_minmax = tile_dir + '/local_minmax.txt'

    z = cfg['subsampling_factor']
    common.cropImage(img1,crop_ref,col,row,tw,th,z)
    
    if os.path.isfile('%s/local_minmax.txt' % tile_dir) and cfg['skip_existing']:
        print "Compute min/max intensities values on tile %s already done, skip" % (tile_dir)
    else:
        if os.path.isfile('%s/this_tile_is_masked.txt' % (tile_dir)):
            print "tile %s already masked, skip" % (tile_dir)
		
        else:
            common.image_getminmax(crop_ref,local_minmax)



def pointing_correction(tile_info):
    """                        
    Computes pointing corrections
    
    Args: 
        tile_info : a list that provides all you need to process a tile: col,
            row, tw, th, ov, i, j, pos, x, y, w, h, images, NbPairs, cld_msk,
            roi_msk
    """     
    col, row, tw, th, ov, i, j, pos, x, y, w, h, images, NbPairs, cld_msk,
    roi_msk, tile_dir = tile_info

    for i in range(0,NbPairs):
        pair_id = i+1
        paired_tile_dir=tile_dir + 'pair_%d' % pair_id

        img1,rpc1,img2,rpc2 = images[0]['img'],images[0]['rpc'],images[pair_id]['img'],images[pair_id]['rpc']
        
        
        # create a directory for the experiment
        if not os.path.exists(paired_tile_dir):
            os.makedirs(paired_tile_dir)
        
        # redirect stdout and stderr to log file
        if not cfg['debug']:
            fout = open('%s/stdout.log' % paired_tile_dir, 'w', 0)  # '0' for no buffering
            sys.stdout = fout
            sys.stderr = fout
    
        # debug print
        print 'tile %d %d running on process %s' % (x, y,
                                                    multiprocessing.current_process())
    
        # output files
        cwid_msk = '%s/cloud_water_image_domain_mask.png' % (paired_tile_dir)
        pointing = '%s/pointing.txt' % paired_tile_dir
        center = '%s/center_keypts_sec.txt' % paired_tile_dir
        sift_matches = '%s/sift_matches.txt' % paired_tile_dir
        sift_matches_plot = '%s/sift_matches_plot.png' % paired_tile_dir


        # check if the tile is already done
        if os.path.isfile('%s/pointing.txt' % paired_tile_dir) and cfg['skip_existing']:
            print "pointing correction on tile %d %d (pair %d) already done, skip" % (col,row,pair_id)
        else:
        
            # check if the ROI is completely masked (water, or outside the image domain)
            H = np.array([[1, 0, -x], [0, 1, -y], [0, 0, 1]])
            if masking.cloud_water_image_domain(cwid_msk, tw, th, H, rpc1, roi_msk,cld_msk):
                print "Tile masked by water or outside definition domain, skip"
                open("%s/this_tile_is_masked.txt" % paired_tile_dir, 'a').close()
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
                if not cfg['debug']:
                    fout.close()
            else:
            
                # correct pointing error
                # A is the correction matrix and m is the list of sift matches
                A, m = pointing_accuracy.compute_correction(img1, rpc1, img2, rpc2, col,
                                                                row, tw, th)
                if A is not None:
                    np.savetxt(pointing, A)
                if m is not None:
                    np.savetxt(sift_matches, m)
                    np.savetxt(center, np.mean(m[:, 2:4], 0))
                    visualisation.plot_matches_pleiades(img1, img2, rpc1, rpc2, m, col, row, tw, th, sift_matches_plot)
       
        # close logs
        common.garbage_cleanup()
        if not cfg['debug']:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            fout.close()

    return
