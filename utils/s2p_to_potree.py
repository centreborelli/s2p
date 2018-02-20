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



def plys_to_potree(input_plys, output, bin_dir='.'):
    """
    Compute a multi-scale representation of a large point cloud.

    The output file can be viewed with a web browser. This is useful for
    huge point clouds. The input is a list of ply files.

    If PotreeConverter is not available it doesn't fail.

    Args:
        output: path to the output folder
        input_plys: list of paths to ply files
    """
    import os.path
    ply2ascii       = os.path.join(bin_dir, 'plytool/ply2ascii')
    txt2las         = os.path.join(bin_dir, 'PotreeConverter/LAStools/bin/txt2las')
    PotreeConverter = os.path.join(bin_dir, 'PotreeConverter/build/PotreeConverter/PotreeConverter')

    #if (not os.path.exists(ply2ascii)) or (not os.path.exists(txt2las)) or (not os.path.exists(PotreeConverter)) :
    #    return  

    outdir = os.path.dirname(output)

    plys = ' '.join(input_plys)

    las = []
    trash = []

    for p in input_plys:
        # make ascii ply if needed
        ap = tmpfile('.ply', outdir)
        lp = tmpfile('.las', outdir)

        las.append(lp)
        common.run("%s < %s > %s" % (ply2ascii, p, ap)) 
        # convert ply to las because PotreeConverter is not able to read PLY
        common.run("%s -parse xyzRGB -verbose -i  %s -o %s 2>/dev/null" % (txt2las, ap, lp)) 

    # generate potree output
    listfile = tmpfile('.txt', outdir)
    ff = open(listfile, 'w')
    for item in las:
        ff.write("%s\n" % item)
    ff.close()

    common.run("mkdir -p %s" % output)
    resourcedir = os.path.join(bin_dir, 'PotreeConverter/PotreeConverter/resources/page_template')
    common.run("LC_ALL=C %s --list-of-files %s -o %s -p cloud --edl-enabled --material ELEVATION --overwrite --page-template %s" % (PotreeConverter, listfile, output, resourcedir) )

    # clean intermediate files
    for p in garbage:
        common.run("rm %s"%p)


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


def test_for_potree(basedir):
    import os.path
    ply2ascii       = os.path.join(basedir, 'plytool/ply2ascii')
    txt2las         = os.path.join(basedir, 'PotreeConverter/LAStools/bin/txt2las')
    PotreeConverter = os.path.join(basedir, 'PotreeConverter/build/PotreeConverter/PotreeConverter')
    print ('looking for:\n    %s\n    %s\n    %s'%(txt2las, ply2ascii, PotreeConverter) )

    if (not os.path.exists(ply2ascii)) or (not os.path.exists(txt2las)) or (not os.path.exists(PotreeConverter)) :
        print ('not found\n')
	raise common.RunFailure


def produce_potree(s2poutdir, potreeoutdir):
    """
    Produce a single multiscale point cloud for the whole processed region.

    Args:
        tiles: list of tiles dictionaries
    """

    basedir = os.path.dirname(os.path.abspath(__file__))
    test_for_potree(os.path.join(basedir,'PotreeConverter_PLY_toolchain/'))

    tiles_file = os.path.join(s2poutdir, 'tiles.txt')
    
    # Read the tiles file
    tiles = s2p.read_tiles(tiles_file)
    print(str(len(tiles))+' tiles found')


    def plyvertex(fname):
        with open(fname) as f:
            for x in f:
                if x.split()[0] == 'element' and x.split()[1] == 'vertex':
                    return int(x.split()[2])


    # collect all plys
    plys = []
    for t in tiles:
        clo = os.path.join(os.path.abspath(os.path.dirname(t)), 'cloud.ply')
        if os.path.isfile(clo):
            if plyvertex(clo) > 0 :
                plys.append(clo)
#    plys = [os.path.join(os.path.abspath(os.path.dirname(t)), 'cloud.ply') for t in tiles if os.path.isfile(os.path.join(os.path.abspath(os.path.dirname(t)), 'cloud.ply'))]

    # produce the potree point cloud
    plys_to_potree(plys, os.path.join(potreeoutdir, 'cloud.potree'), 
		os.path.join(basedir, 'PotreeConverter_PLY_toolchain/'))


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: potree generation tool'))
    
    parser.add_argument('s2pout',metavar='s2poutdir',
                        help=('path to the s2p output directory'))
    parser.add_argument('potreeoutdir',metavar='potreeoutdir', default='',nargs='?',
                        help=('path to output potree (default: current dir)'))
    args = parser.parse_args()

    try:
    	produce_potree(args.s2pout,args.potreeoutdir)
    except common.RunFailure:
        basedir = os.path.dirname(os.path.abspath(__file__))
        print ('You must download and compile PotreeConverter. Run the following commands:')
	print ('    > cd %s'%basedir)
	print ('    > git clone https://github.com/gfacciol/PotreeConverter_PLY_toolchain --recurse-submodules')
	print ('    > cd PotreeConverter_PLY_toolchain')
	print ('    > CC=gcc CXX=g++ make')
