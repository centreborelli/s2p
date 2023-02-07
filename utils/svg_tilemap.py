#!/usr/bin/env python

# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
import os.path
import json
import argparse
import numpy as np
import datetime
import subprocess


import s2p
from s2p.config import cfg
from s2p import common
from s2p import initialization


def write_svg_tilemap(filename, cfg, tiles):
    '''
    draw tiles boundaries with names in an svg file
    '''
    with open(filename,'w') as f:
        f.write('<?xml version="1.0" standalone="no"?>\
        <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" \
         "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">')
        f.write('<svg version="1.1" xmlns="http://www.w3.org/2000/svg" \
               xmlns:xlink="http://www.w3.org/1999/xlink" \
               width="1000px" height="1000px" viewBox="%d %d %d %d">'%(
           cfg['roi']['x'], cfg['roi']['y'], cfg['roi']['w'], cfg['roi']['h']))
        for t in tiles:
            x, y ,w ,h = t['coordinates']
            dir = os.path.abspath(t['dir']).split(os.path.abspath(cfg['out_dir']))[-1]
            try:
               common.image_qauto("%s/dsm.tif"%t['dir'], "%s/dsm.tif.png"%t['dir'])
            except subprocess.CalledProcessError:
               pass

            f.write('<polygon style="fill:white;stroke:black;stroke-width:2" \
                  points="%d %d, %d %d, %d %d, %d %d" />'%(
                  x,y, x+w,y, x+w,y+h, x,y+h))
            f.write('<image xlink:href="./%s/dsm.tif.png" \
                  x="%d" y="%d" width="%d" height="%d"/>'%(
                  dir, x,y, w, h))

            f.write('<a xlink:href="./%s/" target="_blank">'%dir)
            f.write('<g transform="translate(%d,%d)">\
                  <text x="%d" y="%d"  text-anchor="middle" \
                  style="fill: #00FF00; stroke: #00FF00; stroke-width: 0.5; font-size: 12px;" \
                  alignment-baseline="central">%s</text></g>'%(
                     x,y,w/2,h/2,dir))
            f.write('</a>')

        f.write('</svg>')
        f.close()

def main(user_cfg):
    """
    Recompute the s2p tile geometry for the config file
    and produce an svg representing the tiles

    Args:
        user_cfg: user config dictionary
    """
    common.reset_elapsed_time()
    initialization.build_cfg(user_cfg)

    tw, th = initialization.adjust_tile_size()
    tiles_txt = os.path.join(cfg['out_dir'],'tiles.txt')
    tiles = initialization.tiles_full_info(tw, th, tiles_txt)

    # generate svg tile map
    write_svg_tilemap(os.path.join(cfg['out_dir'],'tiles.svg'), cfg, tiles)

    print("\n\n    svg tilemap saved in: %s\n"%os.path.join(cfg['out_dir'],'tiles.svg'))


    # cleanup
    common.garbage_cleanup()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('generate svg tilemap for S2P'))
    parser.add_argument('config', metavar='config.json',
                        help=('path to a json file containing the paths to '
                              'input and output files and the s2p algorithm '
                              'parameters'))
    args = parser.parse_args()

    user_cfg = s2p.read_config_file(args.config)

    main(user_cfg)
