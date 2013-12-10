import piio
import numpy as np
import common

def mosaic(fout, w, h, ov, list_tiles):
    """
    Compose several tiles of the same size into a bigger image.

    Args:
        fout: path to the output image
        w, h: output image dimensions
        ov: overlap between tiles (in pixels)
        list_tiles: list containing paths to the input tiles

    Returns:
        nothing
    """
    N = len(list_tiles)
    tw, th = common.image_size(list_tiles[0])
    ntx = np.ceil(w / (tw - ov)).astype(int)
    nty = np.ceil(h / (th - ov)).astype(int)
    assert(ntx * nty == N)

    out = np.zeros([h, w])
    count = np.zeros([h, w])

    # loop over all the tiles
    for j in range(nty):
        for i in range(ntx):
            # top-left and bottom-right corners of the tile
            x0 = i * (tw - ov)
            y0 = j * (th - ov)
            x1 = min(x0 + tw, w)
            y1 = min(y0 + th, h)

            # read the tile with piio
            tile = piio.read(list_tiles[j * ntx + i])[:, :, 0]
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
