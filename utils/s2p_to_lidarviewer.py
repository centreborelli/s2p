#!/usr/bin/env python
# Copyright (C) 2017, Gabriele Facciolo <gfacciol@gmail.com>
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

import argparse, json, os
import sys

# This is needed to import from a sibling folder
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import s2p
from s2plib import common




## temp files
garbage = list()
def tmpfile(ext='', tmpdir='tmp'):
    """
    Creates a temporary file in the cfg['temporary_dir'] directory.

    Args:
        ext: desired file extension. The dot has to be included.

    Returns:
        absolute path to the created file

    The path of the created file is added to the garbage list to allow cleaning
    at the end of the pipeline.
    """
    import tempfile

    try:
        os.mkdir(tmpdir)
    except OSError:
        pass
    pass
    fd, out = tempfile.mkstemp(suffix=ext, prefix='s2p_',
                               dir=os.path.expandvars(tmpdir))
    os.close(fd)           # http://www.logilab.org/blogentry/17873
    garbage.append(out)
    return out



def read_tiles(tile_files):
    tiles = []
    tile_file_dir = os.path.dirname(tile_files)
    err_log = os.path.join(outdir, "%s_invalid_tiles.txt" % key)

    with open(tile_files, 'r') as f:
        readlines = list(map(str.strip,
                             f.readlines()))
        with open(err_log, 'w') as ferr:
            for el in readlines:
                t = os.path.dirname(os.path.join(tile_file_dir, el))
                dsm = os.path.join(t, 'dsm.tif')
                message = "ok"
                if os.path.exists(dsm) is False:
                    message = "no dsm"

                if message == "ok":
                    tiles.append(t)

    return tiles


def produce_lidarviewer(s2poutdir, output):
    """
    Produce a single multiscale point cloud for the whole processed region.

    Args:
        tiles: list of tiles dictionaries
    """

    tiles_file = os.path.join(s2poutdir, 'tiles.txt')
    
    # Read the tiles file
    tiles = s2p.read_tiles(tiles_file)
    print(str(len(tiles))+' tiles found')

    # collect all plys
    plys = [os.path.join(os.path.abspath(os.path.dirname(t)), 'cloud.ply') for t in tiles]

	
    nthreads = 4
    plys = ' '.join(plys)
    common.run("LidarPreprocessor -to %s.LidarO -tp %s.LidarP -nt %d %s -o %s" % (output,
                                                                           output,
                                                                           nthreads,
                                                                           plys,
                                                                           output))

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: lidarviewer generation tool'))
    
    parser.add_argument('s2pout',metavar='s2poutdir',
                        help=('path to the s2p output directory'))
    parser.add_argument('outdir',metavar='potreeoutdir', default='',nargs='?',
                        help=('path to output lidarviewer (default: current dir)'))
    args = parser.parse_args()

    try:
    	produce_lidarviewer(args.s2pout,args.outdir)
    except common.RunFailure:
        print ('------------------------------------------------------------------------------------')
        print ('\nYou must install LidarPreprocessor from https://github.com/KeckCAVES/LidarViewer\n')
        print ('------------------------------------------------------------------------------------')

