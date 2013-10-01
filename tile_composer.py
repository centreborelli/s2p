from main_script_pairs2 import main
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

for i in range(canvas_x,canvas_x+canvas_w,tile_w-overlap):
   for j in range(canvas_y,canvas_y+canvas_h,tile_h-overlap):
      tile_exp = "%s_%d_%d"%(exp_name,i,j)
      print tile_exp




      # run
      main(img_name, tile_exp, i, j, tile_w, tile_h)





      # generated output files 
      rect1 = '/tmp/%s1.tif' % (tile_exp)
      rect2 = '/tmp/%s2.tif' % (tile_exp)
      hom1  = '/tmp/%s_hom1' % (tile_exp)
      hom2  = '/tmp/%s_hom2' % (tile_exp)
      rpc1 = '/tmp/%s_rpc1.xml' % (tile_exp)
      rpc2 = '/tmp/%s_rpc2.xml' % (tile_exp)
      rect1_color = '/tmp/%s1_color.tif' % (tile_exp)
      disp    = '/tmp/%s_disp.pgm'   % (tile_exp)
      mask    = '/tmp/%s_mask.png'   % (tile_exp)
      cloud   = '/tmp/%s_cloud.ply'  % (tile_exp)
      height  = '/tmp/%s_height.tif' % (tile_exp)
      rpc_err = '/tmp/%s_rpc_err.tif'% (tile_exp)
      subsampling_file = '/tmp/%s_subsampling.txt' % (tile_exp)




      # project and combine

      # height projected on the reference image
      do_project_height  = 1
      projected1_height  = '/tmp/%s_proj1_height.tif' % (tile_exp)
      projected1_mask    = '/tmp/%s_proj1_mask.tif' % (tile_exp)


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

