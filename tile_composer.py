from main_script_pairs import main
from python import common

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



#### COMPUTE THE RESULTS FOR EACH TILE
for i in range(canvas_x,canvas_x+canvas_w,tile_w-overlap):
   for j in range(canvas_y,canvas_y+canvas_h,tile_h-overlap):
      tile_exp = "%s_%d_%d"%(exp_name,i,j)
      print tile_exp


      # run
      main(img_name, tile_exp, i, j, tile_w, tile_h)



#### COMPOSE THE GENERATED TILES
exper_dir='/tmp/'    # default for main 
for i in range(canvas_x,canvas_x+canvas_w,tile_w-overlap):
   for j in range(canvas_y,canvas_y+canvas_h,tile_h-overlap):
      tile_exp = "%s_%d_%d"%(exp_name,i,j)
      print tile_exp

      # generated output files 
      rect1 = '%s/%s1.tif' % (exper_dir,tile_exp)
      rect2 = '%s/%s2.tif' % (exper_dir,tile_exp)
      hom1  = '%s/%s_hom1' % (exper_dir,tile_exp)
      hom2  = '%s/%s_hom2' % (exper_dir,tile_exp)
      rpc1 = '%s/%s_rpc1.xml' % (exper_dir,tile_exp)
      rpc2 = '%s/%s_rpc2.xml' % (exper_dir,tile_exp)
      rect1_color = '%s/%s1_color.tif' % (exper_dir,tile_exp)
      disp    = '%s/%s_disp.pgm'   % (exper_dir,tile_exp)
      mask    = '%s/%s_mask.png'   % (exper_dir,tile_exp)
      cloud   = '%s/%s_cloud.ply'  % (exper_dir,tile_exp)
      height  = '%s/%s_height.tif' % (exper_dir,tile_exp)
      rpc_err = '%s/%s_rpc_err.tif'% (exper_dir,tile_exp)
      subsampling_file = '%s/%s_subsampling.txt' % (exper_dir,tile_exp)


      # project and combine

      # height projected on the reference image
      do_project_height  = 1
      projected1_height  = '%s/%s_proj1_height.tif' % (exper_dir,tile_exp)
      projected1_mask    = '%s/%s_proj1_mask.tif' % (exper_dir,tile_exp)


      import numpy as np
      # read subsampling factor... 
      subsampling_factor=np.loadtxt(subsampling_file)

      hom00 = common.tmpfile('.txt')
      f = 1.0/subsampling_factor
      h0 = np.dot( np.diag([f,f,1]), common.matrix_translation(-canvas_x, -canvas_y))
      common.matrix_write(hom00, h0)
      tmp_h = common.tmpfile('.tif')
      tmp_m = common.tmpfile('.tif')
      common.run('height_rpc_move %s %s %s %s %s %s %s %s %d %d'%(rpc1,hom1,height,mask,  rpc1,hom00, tmp_h, tmp_m, canvas_w/subsampling_factor, canvas_h/subsampling_factor)) 


      # median filter all the pixels that disapear after a 3x3 clausure
      common.run('morphoop %s median 3 %s '%(tmp_h,projected1_height)) 
      common.run('morphoop %s max 3 - | morphoop - min 3 %s'%(tmp_m,projected1_mask)) 
      common.run('plambda %s %s %s %s "a 0 >    c   b 0 > d nan if  if " | iion - %s'%(tmp_m, projected1_mask, tmp_h, projected1_height, projected1_height)) 
common.run('veco med %s/%s*_proj1_height.tif | iion - %s/%s_FULLHEIGHT.tif'%(exper_dir,exp_name,exper_dir,exp_name))




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
#
#
#      # median filter all the pixels that disapear after a 3x3 clausure
#      common.run('morphoop %s median 3 %s '%(tmp_h,projected1_height)) 
#      common.run('morphoop %s max 3 - | morphoop - min 3 %s'%(tmp_m,projected1_mask)) 
#      common.run('plambda %s %s %s %s "a 0 >    c   b 0 > d nan if  if " | iion - %s'%(tmp_m, projected1_mask, tmp_h, projected1_height, projected1_height)) 
#
#common.run('veco med %s/%s*_proj1_height.tif | iion - %s/%s_FULLHEIGHT.tif'%(exper_dir,basename,exper_dir,basename))





### cleanup
while common.garbage:
    common.run('rm ' + common.garbage.pop())
