# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>

import os.path
import numpy as np
import piio
import common

def mosaic_gdal(fout, w, h, list_tiles, tw, th, ov):
    """
    Compose several tiles of the same size into a bigger image (using gdal vrt)

    Args:
        fout: path to the output image
        w, h: output image dimensions
        list_tiles: list containing paths to the input tiles
        tw, th: dimensions of a tile (they must all have the same dimensions)
        ov: overlap between tiles (in pixels)

    Returns:
        nothing
    """
    N = len(list_tiles)
    ntx = np.ceil(float(w - ov) / (tw - ov)).astype(int)
    nty = np.ceil(float(h - ov) / (th - ov)).astype(int)
    assert(ntx * nty == N)

    vrtfilename = fout+'.vrt'
    
    vrtfile = open(vrtfilename,'w')

    vrtfile.write("<VRTDataset rasterXSize=\"%i\" rasterYSize=\"%i\">\n" %(w,h))
    vrtfile.write("\t<VRTRasterBand dataType=\"Float32\" band=\"1\">\n")
    vrtfile.write("\t\t<ColorInterp>Gray</ColorInterp>\n")

    # loop over all the tiles
    for j in range(nty):
        for i in range(ntx):
            x0 = i * (tw - ov)
            y0 = j * (th - ov)
            x1 = min(x0 + tw -ov, w)
            y1 = min(y0 + th - ov, h)
            tile_fname = list_tiles[j * ntx + i]
            if os.path.isfile(tile_fname):
                vrtfile.write("\t\t<SimpleSource>\n")
                vrtfile.write("\t\t\t<SourceFilename relativeToVRT=\"1\">%s</SourceFilename>\n" %tile_fname)
                vrtfile.write("\t\t\t<SourceBand>1</SourceBand>\n")
                vrtfile.write("\t\t\t<SrcRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" %(0,0,x1-x0,y1-y0))
                vrtfile.write("\t\t\t<DstRect xOff=\"%i\" yOff=\"%i\" xSize=\"%i\" ySize=\"%i\"/>\n" %(x0,y0,x1-x0,y1-y0))
                vrtfile.write("\t\t</SimpleSource>\n")
    
    vrtfile.write("\t</VRTRasterBand>\n")
    vrtfile.write("</VRTDataset>\n")
    vrtfile.close()

    common.run('gdal_translate %s %s' %(vrtfilename,fout))

    return

def mosaic(fout, w, h, list_tiles, tw, th, ov):
    """
    Compose several tiles of the same size into a bigger image.

    Args:
        fout: path to the output image
        w, h: output image dimensions
        list_tiles: list containing paths to the input tiles
        tw, th: dimensions of a tile (they must all have the same dimensions)
        ov: overlap between tiles (in pixels)

    Returns:
        nothing
    """
    N = len(list_tiles)
    ntx = np.ceil(float(w - ov) / (tw - ov)).astype(int)
    nty = np.ceil(float(h - ov) / (th - ov)).astype(int)
    assert(ntx * nty == N)

    out = np.zeros([h, w])
    count = np.zeros([h, w])

    # loop over all the tiles
    for j in range(nty):
        for i in range(ntx):
            # top-left and bottom-right corners of the tile in the output full
            # image
            x0 = i * (tw - ov)
            y0 = j * (th - ov)
            x1 = min(x0 + tw, w)
            y1 = min(y0 + th, h)

            # read the tile with piio. If the tile has not been produced,
            # nothing needs to be done. The corresponding pixels will get the
            # value 'nan' in the output full image.
            tile_fname = list_tiles[j * ntx + i]
            if os.path.isfile(tile_fname):
                tile = piio.read(tile_fname)[:, :, 0]
                assert(np.shape(tile) == (th, tw))

                # count the pixels different from nan and inf
                ind = np.isfinite(tile)
                count[y0:y1, x0:x1] += ind[:y1 - y0, :x1 - x0]

                # replace nan and inf with zeros, then add the tile to the output.
                # ~ind is the negation of ind
                tile[~ind] = 0
                out[y0:y1, x0:x1] += tile[:y1 - y0, :x1 - x0]

    # put nan wher count is zero, and take the average where count is nonzero.
    ind = (count > 0)
    out[ind] = out[ind] / count[ind]
    out[~ind] = np.nan

    piio.write(fout, out)
