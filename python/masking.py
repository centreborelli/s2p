# Copyright (C) 2013, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2013, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2013, Julien Michel <julien.michel@cnes.fr>

import common


def cloud_water_image_domain(out, w, h, H, rpc, roi_gml=None, cld_gml=None):
    """
    Computes a mask for pixels masked by clouds, water, or out of image domain.

    Args:
        out: path to the output image file.
        w, h: (w, h) are the dimensions of the output image mask.
        H: 3x3 numpy array representing the homography that transforms the
            original full image into the rectified tile.
        rpc: paths to the xml file containing the rpc coefficients of the image.
            RPC model is used with SRTM data to derive the water mask.
        roi_gml (optional, default None): path to a gml file containing a mask
            defining the area contained in the full image.
        cld_gml (optional, default None): path to a gml file containing a mask
            defining the areas covered by clouds.

    Returns:
        True if the tile is completely masked, False otherwise.
    """
    # put the coefficients of the homography in a string
    hij = ' '.join(['%f' % x for x in H.flatten()])

    # image domain mask
    if roi_gml is None:  # initialize to 255
        common.run('plambda zero:%dx%d "x 255 +" -o %s' % (w, h, out))
    else:
        common.run('cldmask %d %d -h "%s" %s %s' % (w, h, hij, roi_gml, out))
        if common.is_image_black(out):  # if we are already out, return
            return True

    # cloud mask
    if cld_gml is not None:
        cld_msk = common.tmpfile('.png')
        common.run('cldmask %d %d -h "%s" %s %s' % (w, h, hij, cld_gml,
                                                    cld_msk))
        # cld msk has to be inverted.
        # TODO: add flag to the cldmask binary, to avoid using read/write the
        # msk one more time for this
        common.run('plambda %s "255 x -" -o %s' % (cld_msk, cld_msk))

        intersection(out, out, cld_msk)

    # water mask
    water_msk = common.tmpfile('.png')
    common.run('watermask %d %d -h "%s" %s %s' % (w, h, hij, rpc, water_msk))
    intersection(out, out, water_msk)

    return common.is_image_black(out)


def intersection(out, in1, in2):
    """
    Computes the intersection between two mask files

    Args:
        out: path to the ouput mask image file
        in1, in2: paths to the input mask image files

    The masks are binary. Pixels may have value 0 or 255, 0 being the rejection
    value. The output mask rejects any pixel that is rejected in one input mask,
    or in both input masks. As 0 is the rejection value, the intersection is
    equivalent to a pixelwise product.
    """
    common.run('plambda %s %s "x y 255 / *" -o %s' % (in1, in2, out))


def erosion(out, msk, radius):
    """
    Erodes the accepted regions (ie eliminates more pixels)

    Args:
        out: path to the ouput mask image file
        msk: path to the input mask image file
        radius (in pixels): size of the disk used for the erosion
    """
    if radius >= 2:
        common.run('morsi disk%d erosion %s %s' % (int(radius), msk, out))
