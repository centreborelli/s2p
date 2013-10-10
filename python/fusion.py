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
