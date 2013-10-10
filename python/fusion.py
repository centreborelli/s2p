import numpy as np
import common
import piio

def register_heights(im1, im2):
    """
    Affine registration of heights.
    Args:
        im1: first height map
        im2: second height map, to be registered on the first one

    Returns
        path to the registered second height map
    """
    # we search the (a, b) vector that minimizes the following sum (over
    # all the pixels):
    \sum (im1[i] - (a*im2[i]+b))^2

    # first read the images and store them as numpy 1D arrays, removing all the
    # nans and inf
    i1 = piio.read(im1).ravel() #np.ravel() gives a 1D view
    i2 = piio.read(im2).ravel()
    ind = np.logical_and(np.isfinite(i1), np.isfinite(i2))
    h1 = i1[ind]
    h2 = i2[ind]

    # it is a least squares minimization problem
    A = np.vstack((h2, h2*0+1)).T
    b = h1
    z = np.linalg.lstsq(A, b)[0]
    a = z[0]
    b = z[1]

    # apply the affine transform and return the modified im2
    out = common.tmpfile('.tif')
    common.run('plambda %s "x %f * %f +" > %s' % (im2, a, b, out))
    return out


def merge(im1, im2, thresh, out):
    """
    Args:
        im1, im2: paths to the two input images
        thresh: distance threshold on the intensity values
        out: path to the output image

    This function merges two images. They are supposed to be two height maps,
    sampled on the same grid. If a pixel has a valid height (ie not NAN) in
    only one of the two maps, then we keep this height. When two heights are
    available, if they differ less than the threshold we take the mean, if not
    we discard the pixel (ie assign NAN to the output pixel).
    """
    # the following plambda expression implements:
    #if x == nan
    #    return y
    #if y == nan
    #    return x
    #if fabs(x - y) < t
    #    return (x+y)/2
    #else
    #    return nan
    common.run('plambda %s %s "x isnan y y isnan x x y - fabs %f < x y + 2 / nan if if if" > %s' % (im1, im2, thresh, out)
