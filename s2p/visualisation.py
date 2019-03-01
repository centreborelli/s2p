# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

from __future__ import print_function
import numpy as np
import os
import rasterio

from s2p import common
from s2p import sift
from s2p import estimation
from s2p import rpc_model
from s2p import rpc_utils
import s2p.pointing_accuracy


def plot_line(im, x1, y1, x2, y2, colour):
    """
    Plots a line on a rgb image stored as a numpy array

    Args:
        im: 3D numpy array containing the image values. It may be stored as
            uint8 or float32.
        x1, y1, x2, y2: integer coordinates of the line endpoints
        colour: list of length 3 giving the colour used for the plotted line
            (ie [r, g, b])

    Returns:
        a copy of the input numpy array, with the plotted lines on it. It means
        that the intensities of pixels located on the plotted line are changed.
    """
    # colour points of the line pixel by pixel. Loop over x or y, depending on
    # the biggest dimension.
    if np.abs(x2 - x1) >= np.abs(y2 - y1):
        n = np.abs(x2 - x1)
        for i in range(int(n+1)):
            x = int(x1 + i * (x2 - x1) / n)
            y = int(np.round(y1 + i * (y2 - y1) / n))
            try:
                im[y, x] = colour
            except IndexError:
                pass
    else:
        n = np.abs(y2 - y1)
        for i in range(int(n+1)):
            y = int(y1 + i * (y2 - y1) / n)
            x = int(np.round(x1 + i * (x2 - x1) / n))
            try:
                im[y, x] = colour
            except IndexError:
                pass

    return im


def plot_matches_low_level(im1, im2, matches):
    """
    Displays two images side by side with matches highlighted

    Args:
        im1, im2: paths to the two input images
        matches: 2D numpy array of size 4xN containing a list of matches (a
            list of pairs of points, each pair being represented by x1, y1, x2,
            y2)

    Returns:
        path to the resulting image, to be displayed
    """
    # load images
    with rasterio.open(im1, 'r') as f:
        img1 = f.read().squeeze()
    with rasterio.open(im2, 'r') as f:
        img2 = f.read().squeeze()

    # transform single channel to 3-channels
    if img1.ndim < 3:
        img1 = np.dstack([img1] * 3)
    if img2.ndim < 3:
        img2= np.dstack([img2] * 3)

    # if images have more than 3 channels, keep only the first 3
    if img1.shape[2] > 3:
        img1 = img1[:, :, 0:3]
    if img2.shape[2] > 3:
        img2 = img2[:, :, 0:3]

    # build the output image
    h1, w1 = img1.shape[:2]
    h2, w2 = img2.shape[:2]
    w = w1 + w2
    h = max(h1, h2)
    out = np.zeros((h, w, 3), np.uint8)
    out[:h1, :w1] = img1
    out[:h2, w1:w] = img2

    # define colors, according to min/max intensity values
    out_min = min(np.nanmin(img1), np.nanmin(img2))
    out_max = max(np.nanmax(img1), np.nanmax(img2))
    green = [out_min, out_max, out_min]
    blue = [out_min, out_min, out_max]

    # plot the matches
    for i in range(len(matches)):
        x1 = matches[i, 0]
        y1 = matches[i, 1]
        x2 = matches[i, 2] + w1
        y2 = matches[i, 3]
        # convert endpoints to int (nn interpolation)
        x1, y1, x2, y2 = list(map(int, np.round([x1, y1, x2, y2])))
        plot_line(out, x1, y1, x2, y2, blue)
        try:
            out[y1, x1] = green
            out[y2, x2] = green
        except IndexError:
            pass
    # save the output image, and return its path
    outfile = common.tmpfile('.png')
    common.rasterio_write(outfile, out)
    return outfile


def plot_matches(im1, im2, matches, x=None, y=None, w=None, h=None,
                 outfile=None):
    """
    Plot matches on Pleiades images

    Args:
        im1, im2: paths to full Pleiades images
        matches: 2D numpy array of size 4xN containing a list of matches (a
            list of pairs of points, each pair being represented by x1, y1, x2,
            y2). The coordinates are given in the frame of the full images.
        x, y, w, h (optional, default is None): ROI in the reference image
        outfile (optional, default is None): path to the output file. If None,
            the file image is displayed using the pvflip viewer

    Returns:
        path to the displayed output
    """
    # if no matches, no plot
    if not matches.size:
        print("visualisation.plot_matches: nothing to plot")
        return

    # read rpcs
    rpc1 = rpc_utils.rpc_from_geotiff(im1)
    rpc2 = rpc_utils.rpc_from_geotiff(im2)

    # determine regions to crop in im1 and im2
    if x is not None:
        x1 = x
    else:
        x1 = np.min(matches[:, 0])

    if y is not None:
        y1 = y
    else:
        y1 = np.min(matches[:, 1])

    if w is not None:
        w1 = w
    else:
        w1 = np.max(matches[:, 0]) - x1

    if h is not None:
        h1 = h
    else:
        h1 = np.max(matches[:, 1]) - y1

    x2, y2, w2, h2 = rpc_utils.corresponding_roi(rpc1, rpc2, x1, y1, w1, h1)
    # x2 = np.min(matches[:, 2])
    # w2 = np.max(matches[:, 2]) - x2
    # y2 = np.min(matches[:, 3])
    # h2 = np.max(matches[:, 3]) - y2

    # # add 20 pixels offset and round. The image_crop_gdal function will round
    # # off the coordinates before it does the crops.
    # x1 -= 20; x1 = np.round(x1)
    # y1 -= 20; y1 = np.round(y1)
    # x2 -= 20; x2 = np.round(x2)
    # y2 -= 20; y2 = np.round(y2)
    # w1 += 40; w1 = np.round(w1)
    # h1 += 40; h1 = np.round(h1)
    # w2 += 40; w2 = np.round(w2)
    # h2 += 40; h2 = np.round(h2)

    # do the crops
    crop1 = common.image_qauto(common.image_crop_gdal(im1, x1, y1, w1, h1))
    crop2 = common.image_qauto(common.image_crop_gdal(im2, x2, y2, w2, h2))

    # compute matches coordinates in the cropped images
    pts1 = matches[:, 0:2] - [x1, y1]
    pts2 = matches[:, 2:4] - [x2, y2]

    # plot the matches on the two crops
    to_display = plot_matches_low_level(crop1, crop2, np.hstack((pts1, pts2)))
    if outfile is None:
        os.system('v %s &' % (to_display))
    else:
        common.run('cp %s %s' % (to_display, outfile))

    return
