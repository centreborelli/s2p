from python import common
from python import fusion
from python import pointing_accuracy
from python import global_params
import main_script

img_name = 'uy1'
exp_name = 'campo1'

canvas_x = 5500
canvas_y = 13000
canvas_w = 7000
canvas_h = 10000

canvas_x = 5500
canvas_y = 13000
canvas_w = 3500
canvas_h = 3500

tile_w=2000
tile_h=2000
overlap=100

exp_dir = '/tmp'

def process_pair(img_name, exp_name, x, y, w, h, tile_w=1000, tile_h=1000,
        overlap=100, reference_image_id=1, secondary_image_id=2):
    """
    Computes a height map from a Pair of Pleiades images, using tiles.

    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        tile_w, tile_h: dimensions of the tiles
        overlap: width of overlapping bands between tiles
        reference_image_id: id (1, 2 or 3 for the tristereo datasets, and 1 or
            2 for the bistereo datasets) of the image used as the reference
            image of the pair
        secondary_image_id: id of the image used as the secondary image of the
            pair

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # compute global pointing correction (on the whole ROI)
    A = pointing_accuracy.compute_correction(img_name, exp_name, x, y, w, h,
        reference_image_id, secondary_image_id)

    # generate the tiles
    for i in range(x, x + w, tile_w - overlap):
        for j in range(y, y + h, tile_h - overlap):
            print i, j
            tile_exp = '%s_%d_%d' % (exp_name, i, j)
            height = main_script.process_pair(img_name, tile_exp, i, j, tile_w,
                tile_h, reference_image_id, secondary_image_id, A)
            # weird usage of 'crop' with negative coordinates, to do
            # 'anti-crop' (ie pasting a small image onto a bigger image
            # of nans)
            # TODO: manage zoom
            common.run('crop %s %s/%s_height_to_compose.tif %d %d %d %d' % (height,
                exp_dir, tile_exp, x-i, y-j, w, h))

    # compose the tiles
    out = '%s/%s_height_full.tif' % (exp_dir, exp_name)
    common.run('veco med %s/%s_*_height_to_compose.tif | iion - %s' % (exp_dir,
        exp_name, out))

    return out


def process_triplet(img_name, exp_name, x=None, y=None, w=None, h=None,
        tile_w=1000, tile_h=1000, overlap=100, reference_image_id=2,
        left_image_id=1, right_image_id=3):
    """
    Computes a height map from three Pleiades images.

    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        tile_w, tile_h: dimensions of the tiles
        overlap: width of overlapping bands between tiles
        reference_image_id: id (1, 2 or 3) of the image used as the reference
            image of the triplet
        left_image_id: id of the image used as the secondary image of the first
            pair.
        right_image_id: id of the image used as the secondary image of the
            second pair.

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        rpc = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
        prev = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, reference_image_id)
        x, y, w, h = common.get_roi_coordinates(rpc, prev)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # process the two pairs
    exp_name_left = '%s_left' % exp_name
    h_left = process_pair(img_name, exp_name_left, x, y, w, h, tile_w, tile_h,
        overlap, reference_image_id, left_image_id)

    exp_name_right = '%s_right' % exp_name
    h_right = process_pair(img_name, exp_name_right, x, y, w, h, tile_w,
        tile_h, overlap, reference_image_id, right_image_id)

    # merge the two height maps
    h = '%s/%s_height_merged.tif' % (exp_dir, exp_name)
    fusion.merge(h_left, h_right, 3, h)

    return h

#      # project and combine
#
#      # height projected on the reference image
#      do_project_height  = 1
#      projected1_height  = '%s/%s_proj1_height.tif' % (exper_dir,tile_exp)
#      projected1_mask    = '%s/%s_proj1_mask.tif' % (exper_dir,tile_exp)
#
#
#      import numpy as np
#      # read subsampling factor... 
#      subsampling_factor=np.loadtxt(subsampling_file)
#
#      hom00 = common.tmpfile('.txt')
#      f = 1.0/subsampling_factor
#      h0 = np.dot( np.diag([f,f,1]), common.matrix_translation(-canvas_x, -canvas_y))
#      common.matrix_write(hom00, h0)
#      tmp_h = common.tmpfile('.tif')
#      tmp_m = common.tmpfile('.tif')
#      common.run('height_rpc_move %s %s %s %s %s %s %s %s %d %d'%(rpc1,hom1,height,mask,  rpc1,hom00, tmp_h, tmp_m, canvas_w/subsampling_factor, canvas_h/subsampling_factor)) 
#
#
#      # median filter all the pixels that disapear after a 3x3 clausure
#      common.run('morphoop %s median 3 %s '%(tmp_h,projected1_height)) 
#      common.run('morphoop %s max 3 - | morphoop - min 3 %s'%(tmp_m,projected1_mask)) 
#      common.run('plambda %s %s %s %s "a 0 >    c   b 0 > d nan if  if " | iion - %s'%(tmp_m, projected1_mask, tmp_h, projected1_height, projected1_height)) 




##### COMPOSE PREVIOUSLY GENERATED TILES
#
#names = ['campo1_4x_msmw_11400_13000', 'campo1_4x_msmw_11400_21850', 'campo1_4x_msmw_5500_18900', 'campo1_4x_msmw_8450_15950',
#'campo1_4x_msmw_11400_15950', 'campo1_4x_msmw_5500_13000', 'campo1_4x_msmw_5500_21850', 'campo1_4x_msmw_8450_18900',
#'campo1_4x_msmw_11400_18900', 'campo1_4x_msmw_5500_15950', 'campo1_4x_msmw_8450_13000', 'campo1_4x_msmw_8450_21850']
#
#basename='campo1_4x_msmw'
#exper_dir='/tmp/campo1_4x_msmw/'
#
#for tile_exp in names:
#      print tile_exp
#
#      # generated output files 
#      rect1 = '%s/%s1.tif' % (exper_dir,tile_exp)
#      rect2 = '%s/%s2.tif' % (exper_dir,tile_exp)
#      hom1  = '%s/%s_hom1' % (exper_dir,tile_exp)
#      hom2  = '%s/%s_hom2' % (exper_dir,tile_exp)
#      rpc1 = '%s/%s_rpc1.xml' % (exper_dir,tile_exp)
#      rpc2 = '%s/%s_rpc2.xml' % (exper_dir,tile_exp)
#      rect1_color = '%s/%s1_color.tif' % (exper_dir,tile_exp)
#      disp    = '%s/%s_disp.pgm'   % (exper_dir,tile_exp)
#      mask    = '%s/%s_mask.png'   % (exper_dir,tile_exp)
#      cloud   = '%s/%s_cloud.ply'  % (exper_dir,tile_exp)
#      height  = '%s/%s_height.tif' % (exper_dir,tile_exp)
#      rpc_err = '%s/%s_rpc_err.tif'% (exper_dir,tile_exp)
#      subsampling_file = '%s/%s_subsampling.txt' % (exper_dir,tile_exp)
#
#
#      # project and combine
#
#      # height projected on the reference image
#      do_project_height  = 1
#      projected1_height  = '%s/%s_proj1_height.tif' % (exper_dir,tile_exp)
#      projected1_mask    = '%s/%s_proj1_mask.tif' % (exper_dir,tile_exp)
#
#
#      import numpy as np
#      # read subsampling factor... 
#      subsampling_factor=np.loadtxt(subsampling_file)
#
#      hom00 = common.tmpfile('.txt')
#      f = 1.0/subsampling_factor
#      h0 = np.dot( np.diag([f,f,1]), common.matrix_translation(-canvas_x, -canvas_y))
#      common.matrix_write(hom00, h0)
#      tmp_h = common.tmpfile('.tif')
#      tmp_m = common.tmpfile('.tif')
#      common.run('height_rpc_move %s %s %s %s %s %s %s %s %d %d'%(rpc1,hom1,height,mask,  rpc1,hom00, tmp_h, tmp_m, canvas_w/subsampling_factor, canvas_h/subsampling_factor)) 




### cleanup
#while common.garbage:
#    common.run('rm ' + common.garbage.pop())
