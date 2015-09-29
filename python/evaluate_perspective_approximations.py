# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import rpc_model
import rpc_utils
import estimation
import evaluation


def single_reprojection_error(rpc, proj, col, row, h):
    """
    Compute the distance, in the image plane, between two projections of the
    same artificial ground control point. The ground control point is obtained
    thanks to the RPC localisation function.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        proj: instance of the projective_model.ProjModel class
        col, row: coordinates of the pixel defining the ground control point
        h: altitude of the ground control point

    Returns:
        A distance (in pixels) measured in the image plane
    """

    # geodetic coordinates of the ground control point
    lon, lat, alt = rpc.direct_estimate(col, row, h)

    # find back col and row using the rpc, then the projective model
    col_rpc, row_rpc, alt_rpc = rpc.inverse_estimate(lon, lat, alt)
    col_proj, row_proj, alt_proj = proj.inverse_estimate(lon, lat, alt)

    # compute the distance
    return ((col_proj - col_rpc)**2 +  (row_proj - row_rpc)**2) ** 0.5


def camera_error(rpc, x0, y0, w, h, n_learn, n_test):
    """
    Measure the accuracy of RPC's projective approximation.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        x0, y0, w, h: four integers definig a rectangular region of interest
            (ROI) in the first view. (x0, y0) is the top-left corner, and (w,
            h) are the dimensions of the rectangle. The keypoints used for
            estimation and evaluation are located in that ROI.
        n_learn: number of points sampled in each of the 3 directions to
            generate 3D-->2D correspondences constraints for projective matrix
            estimation
        n_test: number of points sampled in each of the 3 directions to
            generate 3D points used to test the accuracy of the projective
            matrix
    Returns:
        A float value measuring the accuracy of the projective model. Smaller
        is better.

        This value is the biggest euclidean distance, among all the 'test
        points', between a point projected on the image with the approximated
        matrix and the corresponding projected point with the rpc.
    """

    # step 1: compute the projective model that best fits this range
    X, x = rpc_utils.world_to_image_correspondences_from_rpc(rpc,
                                                        x0, y0, w, h, n_learn)
    P = estimation.camera_matrix(X, x)

    # step 2: compute the error made by the projective model
    X, x = rpc_utils.world_to_image_correspondences_from_rpc(rpc,
                                                        x0, y0, w, h, n_test)
    err = evaluation.camera_matrix(X, x)
    return err


def fundamental_error(rpc1, rpc2, x, y, w, h, n_learn, n_test):
    """
    Measure the error made when imposing a fundamental matrix fitting on
    Pleiades data.

    Args:
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig a rectangular region of interest
            (ROI) in the first view. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle. In the first view, the keypoints
            used for estimation and evaluation are located in that ROI.
        n_learn: number of points sampled in each of the 3 directions to
            generate image correspondences constraints for fundamental matrix
            estimation.
        n_test: number of points sampled in each of the 3 directions to
            generate correspondences used to test the accuracy of the fundamental
            matrix.

    Returns:
        A float value measuring the accuracy of the fundamental matrix. Smaller
        is better.
    """

    # step 1: estimate the fundamental matrix
    matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, n_learn)
    F = estimation.fundamental_matrix(matches)

    # step 2: compute the residual error
    matches = rpc_utils.matches_from_rpc(rpc1, rpc2, x, y, w, h, n_test)
    err = evaluation.fundamental_matrix(F, matches)
    return err


def plot_projective_error(n_learn, n_test, rpc1, rpc2=None):
    """
    plot a 3D graphic showing projective error, depending on position and size.

    Args:
        n_learn: number of points sampled in each of the 3 directions to
            generate 3D-->2D correspondences constraints for projective matrix
            estimation
        n_test: number of points sampled in each of the 3 directions to
            generate 3D points used to test the accuracy of the projective
            matrix
        rpc1/2: instance of the rpc_model.RPCModel class
           If there is only one rpc, then compute the error associated to the
           approximated camera matrix. If two rpc are provided, then compute
           the epipolar error (associated to the estimated fundamental matrix)

    Returns:
        numpy 2D array containing all the computed errors
    """
    # center of the image
    p0 = 0.5 * np.array([rpc1.lastCol, rpc1.lastRow])
    p0 = np.rint(p0).astype(int) # round and convert to int

    # params
    r_max = 10000.0
    r_step = 1000.0
    w_min = 1000.0
    w_max = 4000.0
    w_step = 500.0
    angles = np.pi * np.arange(8)
    w_nb = np.ceil((w_max - w_min)/w_step)
    r_nb = np.ceil(r_max/r_step)
    err = np.zeros((w_nb, r_nb))

    # loop on radius and roi sizes
    # for each pair (w, r), find the maximal error we can obtain with the
    # function camera_error, or fundamental_error. For a fixed r, we
    # explore 8 positions on the corresponding circle and choose the one with
    # the highest error.
    for idx_w, w in enumerate(np.arange(w_min, w_max, w_step)):

        # first case: r=0
        col, row = p0[0], p0[1]
        x = col - 0.5*w
        y = row - 0.5*w
        if rpc2:
            err[idx_w, 0] = fundamental_error(rpc1, rpc2, x, y, w, w, n_learn,
                                                                         n_test)
        else:
            err[idx_w, 0] = camera_error(rpc1, x, y, w, w, n_learn, n_test)


        # second case: r>0
        for idx_r, r in enumerate(np.arange(r_step, r_max, r_step)):
            err_on_circle = np.zeros(8)
            for idx_t, theta in enumerate(angles):
                # implements [p0+r*cos(theta), p0+r*sin(theta)], with roundoff
                v = np.array([np.cos(theta), np.sin(theta)])
                v = np.rint(p0 + r*v).astype(int)
                col, row = v[0], v[1]
                x = col - 0.5*w
                y = row - 0.5*w
                if rpc2:
                    err_on_circle[idx_t] = fundamental_error(rpc1, rpc2, x, y,
                                                          w, w, n_learn, n_test)
                else:
                    err_on_circle[idx_t] = camera_error(rpc1, x, y, w, w,
                                                                n_learn, n_test)
            err[idx_w, idx_r+1] = err_on_circle.max()

    plot_3d_surface(np.arange(0, r_max, r_step),
            np.arange(w_min, w_max, w_step), err)
    return err


def plot_3d_surface(x, y, z, xlabel=None, ylabel=None, zlabel=None,
        outfile=None):
    """
    Call the Axes3D.plot_trisurf function from matplotlib to plot a surface.
    Args:
        x, y: 1D arrays containing the x and y values
        z: 2D array containing the z values. BE CAREFUL, its shape has to
            be (len(y), len(x))
        {x,y,z}label: string containing label for the {x,y,z} axis
        outfile (optional): path to the pdf file where to save the plot

    Returns:
        nothing, but opens a window with the plot
    """
    # reshape the data
    x, y = np.meshgrid(x, y)
    x = x.ravel()
    y = y.ravel()
    z = z.ravel()

    # make the plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z, cmap = cm.jet, linewidth = 0.2)

    # put labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if zlabel is not None:
        ax.set_zlabel(zlabel)

    # save or show the plot
    if outfile is not None:
        pp = PdfPages(outfile)
        plt.savefig(pp, format='pdf', bbox_inches='tight')
        pp.close()
    else:
        plt.show()


def main():
    """
    launch the compute_projective_error method with default params
    """
    # default parameters
    rpc = rpc_model.RPCModel('../rpc_data/haiti/rpc01.xml')
    col, row, w, h = 15000, 20000, 2000, 1000

    # check parameters
    if len(sys.argv) > 4:
        col = int(sys.argv[1])
        row = int(sys.argv[2])
        w = int(sys.argv[3])
        h = int(sys.argv[4])
    elif len(sys.argv) > 1:
        print "Incorrect syntax, use:"
        print '  > ' + sys.argv[0] + " col row w h"
        sys.exit(1)

    err = camera_error(rpc, col, row, w, h, 5, 7)
    print 'projective error: ', err

if __name__ == '__main__': main()
