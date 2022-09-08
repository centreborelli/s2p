# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>

import numpy as np
import rasterio

from s2p import common
from s2p import rpc_utils

from rasterio.windows import Window

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


def plot_matches_low_level(crop1, crop2, matches, outfile):
    """
    Displays two images side by side with matches highlighted

    Args:
        crop1, crop2 (np.array): two input images
        matches: 2D numpy array of size 4xN containing a list of matches (a
            list of pairs of points, each pair being represented by x1, y1, x2,
            y2)
        outfile (str): path where to write the resulting image, to be displayed
    """
    # transform single channel to 3-channels
    if crop1.shape[0] < 3:
        crop1 = np.concatenate([crop1] * 3, axis=0)
    if crop2.shape[0] < 3:
        crop2 = np.concatenate([crop2] * 3, axis=0)

    # if images have more than 3 channels, keep only the first 3
    if crop1.shape[0] > 3:
        crop1 = crop1[:3, :, :]
    if crop2.shape[0] > 3:
        crop2 = crop2[:3, :, :]

    # build the output image
    h1, w1 = crop1.shape[1:]
    h2, w2 = crop2.shape[1:]
    w = w1 + w2
    h = max(h1, h2)
    out = np.zeros((3, h, w), np.uint8)
    out[:, :h1, :w1] = crop1
    out[:, :h2, w1:w] = crop2
    out = np.transpose(out, [1, 2, 0])

    # define colors, according to min/max intensity values
    out_min = min(np.nanmin(crop1), np.nanmin(crop2))
    out_max = max(np.nanmax(crop1), np.nanmax(crop2))
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

    # save output image
    common.rasterio_write(outfile, out)


def plot_matches(img1, img2, rpc1, rpc2, matches, outfile, x, y, w, h):
    """
    Plot keypoint matches on images corresponding ROIs.

    Args:
        im1, im2: paths to full Pleiades images
        rpc1, rpc2: two instances of the rpcm.RPCModel class
        matches: 2D numpy array of size 4xN containing a list of matches (a
            list of pairs of points, each pair being represented by x1, y1, x2,
            y2). The coordinates are given in the frame of the full images.
        outfile: path to the output file
        x, y, w, h (ints): ROI in the reference image

    Returns:
        path to the displayed output
    """
    # if no matches, no plot
    if not matches.size:
        print("visualisation.plot_matches: nothing to plot")
        return

    x1, y1, w1, h1 = x, y, w, h
    x2, y2, w2, h2 = map(int, rpc_utils.corresponding_roi(rpc1, rpc2, x1, y1, w1, h1))

    # read the crops
    with rasterio.open(img1, "r") as f:
        crop1 = f.read(window=Window(x1, y1, w1, h1))
    with rasterio.open(img2, "r") as f:
        crop2 = f.read(window=Window(x2, y2, w2, h2))

    crop1 = common.linear_stretching_and_quantization_8bit(crop1)
    crop2 = common.linear_stretching_and_quantization_8bit(crop2)

    # compute matches coordinates in the cropped images
    pts1 = matches[:, :2] - [x1, y1]
    pts2 = matches[:, 2:] - [x2, y2]

    # plot the matches on the two crops
    plot_matches_low_level(crop1, crop2, np.hstack((pts1, pts2)), outfile)
