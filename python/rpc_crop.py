#!/usr/bin/env python

# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>
# Copyright (C) 2013, Enric Meinhardt Llopis <enric.meinhardt@cmla.ens-cachan.fr>

import os
import sys

import rpc_model
import rpc_utils
import common


def procedure1(poly, a11, a22, b1, b2):
   """
   given the polynomial of degree three of three variables poly(x,y,z)
   computes the polynomial coefficients implementing
   the variable change x = a11 x' + b1
                       y = a22 y' + b2
   here the variable z is height of the RPC model
   VERIFIED!
   """
   newpoly = [None] * 20
   newpoly[0]    =  poly[0] + b2*b2*b2* poly[11] + b1*b1* b2* poly[12] + b1* b2*b2* poly[14] + b1*b1*b1* poly[15] + b2* poly[1] + b1* poly[2] + b1* b2* poly[4] + b2*b2* poly[7] + b1*b1* poly[8]
   newpoly[1]    =  (3* a22* b2*b2* poly[11] + a22* b1*b1 *poly[12] + 2* a22* b1* b2* poly[14] + a22* poly[1] + a22* b1* poly[4] + 2* a22* b2* poly[7])
   newpoly[2]    =  (2* a11* b1* b2* poly[12] + a11* b2*b2* poly[14] + 3* a11* b1*b1* poly[15] + a11* poly[2] + a11* b2* poly[4] + 2* a11* b1* poly[8])
   newpoly[3]    =  (b1* b2* poly[10] + b2*b2* poly[17] + b1*b1* poly[18] + poly[3] + b2* poly[5] + b1* poly[6])
   newpoly[4]    =  (2* a11* a22* b1* poly[12] + 2* a11* a22* b2* poly[14] + a11* a22* poly[4])
   newpoly[5]    =  (a22* b1* poly[10] + 2* a22* b2* poly[17] + a22* poly[5])
   newpoly[6]    =  (a11* b2* poly[10] + 2* a11* b1* poly[18] + a11* poly[6])
   newpoly[7]    =  (3* a22*a22* b2* poly[11] + a22*a22* b1* poly[14] + a22*a22* poly[7])
   newpoly[8]    =  (a11*a11* b2* poly[12] + 3 *a11*a11* b1* poly[15] + a11*a11* poly[8])
   newpoly[9]    =  (b2* poly[13] + b1* poly[16] + poly[9])
   newpoly[10]   =  a11* a22* poly[10]
   newpoly[11]   =  a22*a22*a22* poly[11]
   newpoly[12]   =  a11*a11* a22* poly[12]
   newpoly[13]   =  a22* poly[13]
   newpoly[14]   =  a11* a22*a22* poly[14]
   newpoly[15]   =  a11*a11*a11* poly[15]
   newpoly[16]   =  a11* poly[16]
   newpoly[17]   =  a22*a22* poly[17]
   newpoly[18]   =  a11*a11* poly[18]
   newpoly[19]   =  poly[19]
   return newpoly


def poly_variable_change_in(polyNum, polyDen, a11, a22, b1, b2):
   """
   given the RPC polynomials polyNum(x,y,z)/polyDen(x,y,z)
   computes the polynomial coefficients implementing
   the variable change x = a11 x' + b1
                       y = a22 y' + b2
   VERIFIED!
   """
   #print a11,a22,b1,b2
   newNum = procedure1(polyNum,a11, a22, b1, b2)
   newDen = procedure1(polyDen,a11, a22, b1, b2)
   return newNum, newDen


def poly_variable_change_out(polyNum, polyDen, a11, b1):
   """
   given the RPC polynomials polyNum(x,y,z)/polyDen(x,y,z)
   computes the polynomial coefficients implementing
   the operation   a11*(polyNum(x,y,z)/polyDen(x,y,z)) + b1
   VERIFIED!
   """
   import numpy as np
   newNum = list(float(a11) * np.array(polyNum) + float(b1) * np.array(polyDen))
   newDen = polyDen
   return newNum, newDen



def rpc_apply_crop_to_rpc_model(rpc, x0, y0, w, h):
   import copy
   rpcout = copy.deepcopy(rpc)

   ## compute the scale and shift parameter for the normalized RPC
   a11 = (float(w)/2) / (rpc.colScale)
   a22 = (float(h)/2) / (rpc.linScale)
#   a11 = 1.0
#   a22 = 1.0
   b1  = float(x0)/float(rpc.colScale)
   b2  = float(y0)/float(rpc.linScale)
   ## apply the transformation to the direct polynomials, BEWARE!! The order of the variables is reversed
   rpcout.directLonNum, rpcout.directLonDen = poly_variable_change_in(rpc.directLonNum, rpc.directLonDen, a22,a11,b2,b1)
   rpcout.directLatNum, rpcout.directLatDen = poly_variable_change_in(rpc.directLatNum, rpc.directLatDen, a22,a11,b2,b1)


   # scale the RPC domain (I'm not sure its [-1,1]^2)
   #   # TODO correct RPC so that the validity domain is still the square [-1,1]^2
   rpcout.colScale= float(w)/2
   rpcout.linScale= float(h)/2
#   rpcout.colScale= rpc.colScale  ## keep it unchanged (it also works)
#   rpcout.linScale= rpc.linScale


   ## compute the scale and shift parameters for the normalized RPC
   b1 = float(x0)/float(rpcout.colScale)
   b2 = float(y0)/float(rpcout.linScale)
   a11 = float(rpc.colScale)/float(rpcout.colScale)
   a22 = float(rpc.linScale)/float(rpcout.linScale)
#   a11 = 1.0
#   a22 = 1.0
   ## apply the transform to the inverse polynomials
   rpcout.inverseColNum, rpcout.inverseColDen  =  poly_variable_change_out(rpcout.inverseColNum, rpcout.inverseColDen, a11, -b1)
   rpcout.inverseLinNum, rpcout.inverseLinDen  =  poly_variable_change_out(rpcout.inverseLinNum, rpcout.inverseLinDen, a22, -b2)

   return rpcout


def test_me(rpcfile):
   import numpy as np
   import rpc_model, rpc_crop
   reload(rpc_crop)
   r1 = rpc_model.RPCModel('RPC_PHR1A_P_201309231105136_SEN_756965101-001.XML')
   r2 = rpc_crop.rpc_apply_crop_to_rpc_model(r1, 10000,20000,2000,2000)

   #print "direct estimate error:"
   geo1 = np.array(r1.direct_estimate(11000,20100,10, return_normalized=False))
   geo2 = np.array(r2.direct_estimate(1000,100,10, return_normalized=False))
   print geo1 - geo2

   #print "inverse estimate error:"
   pix1 = np.array(r1.inverse_estimate(geo1[0], geo1[1], geo1[2]))
   pix2 = np.array(r2.inverse_estimate(geo1[0], geo1[1], geo1[2]))
   print pix1 - pix2

   r2.write('cropped.xml')


def insert_suffix_filename(filename, suffix):
    """
    Inserts a suffix in a filename, just before the extension.

    Args:
        filename: input file name
        suffix: suffix to be inserted

    Returns:
        The same file name, with the suffix inserted just before the extension
    """
    filename = os.path.basename(filename)
    basename  = filename.split('.')[0]
    extension = filename.split('.')[1]
    return "%s%s.%s" % (basename, suffix, extension)


def crop_rpc_and_image(out_dir, img, rpc, rpc_ref, x, y, w, h):
    """
    Args:
        out_dir: path to the output directory. The cropped image and rpc files
            will be written there.
        img: path to the input image
        rpc: path to the input rpc
        rpc_ref: path to the rpc file of the reference image
        x, y, w, h: 4 integers defining a rectangular ROI in the reference
            image
    """
    r = rpc_model.RPCModel(rpc)

    # recompute the roi if the input image is not the reference image
    if rpc_ref is not rpc:
        r_ref = rpc_model.RPCModel(rpc_ref)
        x, y, w, h = rpc_utils.corresponding_roi(r_ref, r, x, y, w, h)

    # output filenames
    out_img_file = os.path.join(out_dir, insert_suffix_filename(img, "_crop"))
    out_rpc_file = os.path.join(out_dir, insert_suffix_filename(rpc, "_crop"))

    # do the crop
    out_r = rpc_apply_crop_to_rpc_model(r, x, y, w, h)
    out_r.write(out_rpc_file)
    os.system('gdal_translate -srcwin %d %d %d %d "%s" "%s"' % (x, y, w, h, img, out_img_file))



def main():

    # verify input
    if len(sys.argv) in [10, 12]:
       out_dir = sys.argv[1]
       img1 = sys.argv[2]
       rpc1 = sys.argv[3]
       img2 = sys.argv[4]
       rpc2 = sys.argv[5]
       if len(sys.argv) == 10:
           x = float(sys.argv[6])
           y = float(sys.argv[7])
           w = float(sys.argv[8])
           h = float(sys.argv[9])
       else:
           img3 = sys.argv[6]
           rpc3 = sys.argv[7]
           x = float(sys.argv[8])
           y = float(sys.argv[9])
           w = float(sys.argv[10])
           h = float(sys.argv[11])
    else:
       print "Tool to crop an image and its RPC."
       print "Incorrect syntax, use:"
       print "  >  %s out_dir img1 rpc1 img2 rpc2 [img3 rpc3] x y w h" % sys.argv[0]
       exit(1)

    try:
       os.stat(img1)
       os.stat(rpc1)
       os.stat(img2)
       os.stat(rpc2)
    except OSError:
       exit(1)

    # create output dir
    os.system('mkdir -p %s' % out_dir)

    # do the crops
    crop_rpc_and_image(out_dir, img1, rpc1, rpc1, x, y, w, h)
    crop_rpc_and_image(out_dir, img2, rpc2, rpc1, x, y, w, h)
    if 'img3' in locals() and 'rpc3' in locals():
        crop_rpc_and_image(out_dir, img3, rpc3, rpc1, x, y, w, h)

    # clean tmp files
    while common.garbage:
        common.run('rm %s' % common.garbage.pop())


if __name__ == '__main__': main()
