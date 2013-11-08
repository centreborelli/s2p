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
    # remove high frequencies with a morphological zoom out
    im1_low_freq = common.image_zoom_out_morpho(im1, 4)
    im2_low_freq = common.image_zoom_out_morpho(im2, 4)

    # first read the images and store them as numpy 1D arrays, removing all the
    # nans and inf
    i1 = piio.read(im1_low_freq).ravel() #np.ravel() gives a 1D view
    i2 = piio.read(im2_low_freq).ravel()
    ind = np.logical_and(np.isfinite(i1), np.isfinite(i2))
    h1 = i1[ind]
    h2 = i2[ind]

    # for debug
    print np.shape(i1)
    print np.shape(h1)

#    # 1st option: affine
#    # we search the (u, v) vector that minimizes the following sum (over
#    # all the pixels):
#    #\sum (im1[i] - (u*im2[i]+v))^2
#    # it is a least squares minimization problem
#    A = np.vstack((h2, h2*0+1)).T
#    b = h1
#    z = np.linalg.lstsq(A, b)[0]
#    u = z[0]
#    v = z[1]
#
#    # apply the affine transform and return the modified im2
#    out = common.tmpfile('.tif')
#    common.run('plambda %s "x %f * %f +" > %s' % (im2, u, v, out))

    # 2nd option: translation only
    v = np.mean(h1 - h2)
    out = common.tmpfile('.tif')
    common.run('plambda %s "x %f +" > %s' % (im2, v, out))

    return out


def merge(im1, im2, thresh, out):
    """
    Args:
        im1, im2: paths to the two input images
        thresh: distance threshold on the intensity values
        out: path to the output image

    This function merges two images. They are supposed to be two height maps,
    sampled on the same grid. If a pixel has a valid height (ie not inf) in
    only one of the two maps, then we keep this height. When two heights are
    available, if they differ less than the threshold we take the mean, if not
    we discard the pixel (ie assign NAN to the output pixel).
    """
    # first register the second image on the first
    im2 = register_heights(im1, im2)

    # then merge
    # the following plambda expression implements:
    # if isfinite x
    #   if isfinite y
    #     if fabs(x - y) < t
    #       return (x+y)/2
    #     return nan
    #   return x
    # return y
    common.run("""
        plambda %s %s "x isfinite y isfinite x y - fabs %f < x y + 2 / nan if x
        if y if" > %s
        """ % ( im1, im2, thresh, out))
