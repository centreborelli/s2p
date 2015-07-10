#!/usr/bin/env python
# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>


import numpy as np


def triangulate_DLT(P1, P2, x1, x2):
    """
    Triangulate the points x1, x2 using the DLT method (Hartley p312)

    Args:
        P1, P2: 3x4 numpy arrays representing the camera matrices
        x1: array-like containing the 2 coordinates of a pixel in image 1
        x2: array-like containing the 2 coordinates of a pixel in image 2

    Returns:
        X: array-like containing the 4 homogeneous coordinates of the
            triangulated point
        e1, e2: reprojection errors in pixels
    """
    # build homogeneous system
    A  = np.zeros((4,4))
    A[0,:] = x1[0] * P1[2,:] - P1[0,:]
    A[1,:] = x1[1] * P1[2,:] - P1[1,:]
    A[2,:] = x2[0] * P2[2,:] - P2[0,:]
    A[3,:] = x2[1] * P2[2,:] - P2[1,:]

    # solve: the returned s is sorder in descending order and V is transposed
    U, s, V  = np.linalg.svd(A)
    X = V[3,:]/V[3,3]

    # compute the distance of the projected point
    # can do the same for P2
    px1 = np.dot(P1, X)
    px2 = np.dot(P2, X)
    e1 = np.linalg.norm(px1[:2]/px1[2] - x1)
    e2 = np.linalg.norm(px2[:2]/px2[2] - x2)
    return X, e1, e2



def disp_to_h_projective(P1,P2,H1,H2,disp,mask):
    height   = disp.copy()* np.nan
    proj_err = disp.copy()* np.nan

    sz = disp.shape
    w,h = sz[1],sz[0]

    # inverse homographies
    iH1 = np.linalg.inv(H1)
    iH2 = np.linalg.inv(H2)

    # needed fpor the distance to the plane
    M  = P1[:,:3]
    signDetM = np.sign ( np.linalg.det(M) )

    for j in range(h):
       # print something while waiting
       import sys
       str = "%d %d\r"%(j,h)
       sys.stdout.write(str)
       sys.stdout.flush()

       for i in range(w):
          if np.isfinite( disp[j,i] ) and mask[j,i] > 0:

             X1 = iH1.dot(np.array([i          ,j,1]))
             X2 = iH2.dot(np.array([i+disp[j,i],j,1]))

             X3D, error = triangulate_DLT(P1,P2,X1,X2)[:2]
             proj_err[j,i] = error

             # distance to the projective plane. Equation 6.15 (Hartley)
             height[j,i] = P1[2,:].dot(X3D) * signDetM / X3D[3]
    return height, proj_err



def compute_height_map(P1, P2, H1, H2, f_disp, f_mask, fout_height, fout_proj_err):
    import piio
    import common
    P1 = common.matrix_read(P1,3,4)
    P2 = common.matrix_read(P2,3,4)
    H1 = common.matrix_read(H1,3,3)
    H2 = common.matrix_read(H2,3,3)
    disp = piio.read(f_disp)[:,:,0]
    mask = piio.read(f_mask)[:,:,0]

    height,proj_err = disp_to_h_projective(P1,P2,H1,H2,disp,mask)

    piio.write(fout_height  ,height)
    piio.write(fout_proj_err,proj_err)




if __name__ == '__main__':
    """
    disp to h projective
    """
    import sys

    # check parameters
    if len(sys.argv) > 7:
        v=sys.argv
        compute_height_map(v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8])
    else:
        print "Incorrect syntax, use:"
        print '  > ' + sys.argv[0] + "P1 P2 H1 H2 disp mask height proj_err"
        sys.exit(1)

