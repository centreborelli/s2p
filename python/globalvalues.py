# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>


import os
import numpy as np

from config import cfg
from python import pointing_accuracy


def pointing_correction(tiles):
    """
    Compute the global pointing corrections for each pair of images.

    Args:
        tiles: list of tile_info dictionaries
    """
    for i in xrange(1, len(cfg['images'])):
        global_point_file = os.path.join(cfg['out_dir'], 'global_pointing_pair_%d.txt' % i)
        if not (os.path.isfile(global_point_file) and cfg['skip_existing']):
            list_of_tiles = [os.path.join(t['directory'], 'pair_%d' % i) for t in
                             tiles]
            np.savetxt(global_point_file,
                       pointing_accuracy.global_from_local(list_of_tiles), fmt='%12.6f')
