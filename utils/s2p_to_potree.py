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

from __future__ import print_function
import os
import argparse
from codecs import open
import re

from bs4 import BeautifulSoup

import s2p
from s2p import common


def tmpfile(ext='', tmpdir='tmp'):
    """
    Creates a temporary file in the cfg['temporary_dir'] directory.

    Args:
        ext: desired file extension. The dot has to be included.

    Returns:
        absolute path to the created file
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
    return out


def plys_to_potree(input_plys, output, bin_dir='.', cloud_name="cloud"):
    """
    Compute a multi-scale representation of a large point cloud.

    The output file can be viewed with a web browser. This is useful for
    huge point clouds. The input is a list of ply files.

    If PotreeConverter is not available it doesn't fail.

    Args:
        output: path to the output folder
        input_plys: list of paths to ply files
    """
    PotreeConverter = os.path.join(bin_dir, 'PotreeConverter/build/PotreeConverter/PotreeConverter')

    outdir = os.path.dirname(output)

    # List ply files in text file
    listfile = tmpfile('.txt', outdir)
    with open(listfile, 'w') as f:
        for p in input_plys:
            f.write("%s\n" % p)

    # Run PotreeConverter
    common.run("mkdir -p %s" % output)
    resourcedir = os.path.join(bin_dir, 'PotreeConverter/PotreeConverter/resources/page_template')
    common.run("LC_ALL=C %s --list-of-files %s -o %s -p %s --edl-enabled --material RGB --overwrite --page-template %s" % (PotreeConverter, listfile, output, cloud_name, resourcedir))

    # Cleanup
    os.remove(listfile)


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
    PotreeConverter = os.path.join(basedir, 'PotreeConverter/build/PotreeConverter/PotreeConverter')
    print('looking for:\n    %s' % PotreeConverter)

    if not os.path.exists(PotreeConverter):
        print('not found\n')
        raise common.RunFailure


def produce_potree(s2p_outdirs_list, potreeoutdir):
    """
    Produce a single multiscale point cloud for the whole processed region.

    Args:
        s2poutdirs_list: list of s2p output directories
    """
    basedir = os.path.dirname(os.path.abspath(__file__))
    test_for_potree(os.path.join(basedir, 'PotreeConverter_PLY_toolchain/'))

    def plyvertex(fname):
        with open(fname, 'r', 'utf-8') as f:
            for x in f:
                if x.split()[0] == 'element' and x.split()[1] == 'vertex':
                    return int(x.split()[2])

    js_scripts = []
    regex = re.compile("Potree\.loadPointCloud\(.*\);", re.DOTALL)
    cloudoutdir = os.path.join(potreeoutdir, "cloud.potree")

    # Produce a "cloud_?.html" file for all given s2p outdirs
    for i, s2p_outdir in enumerate(s2p_outdirs_list):
        tiles = s2p.read_tiles(os.path.join(s2p_outdir, 'tiles.txt'))
        print(str(len(tiles))+' tiles found')

        # collect all plys
        plys = []
        for t in tiles:
            clo = os.path.join(os.path.abspath(os.path.dirname(t)), 'cloud.ply')
            if os.path.isfile(clo):
                if plyvertex(clo) > 0:
                    plys.append(clo)

        # produce the potree point cloud
        cloud_name = "cloud_{}".format(i)
        plys_to_potree(
            plys,
            cloudoutdir,
            os.path.join(basedir, 'PotreeConverter_PLY_toolchain/'),
            cloud_name,
        )

        # Gather the js script inside the HTML file that is relevant
        # to the point cloud
        cloud_html = os.path.join(cloudoutdir, "{}.html".format(cloud_name))
        with open(cloud_html) as f:
            soup = BeautifulSoup(f, features="lxml")
        script = soup.find_all("script")[-1]
        js_script = re.search(regex, script.text).group(0)
        js_scripts.append(js_script)
        os.remove(cloud_html)

    # The "main.html" file will contain a concatenation of all the js
    # scripts that were gathered in the loop above.

    # Use the last HTML file as a basis for the "main.html", and replace
    # its js script by all the js scripts
    main_html = os.path.join(cloudoutdir, "main.html")
    script.string = re.sub(regex, "\n".join(js_scripts), script.text)
    with open(main_html, "w") as f:
        f.write(soup.prettify())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: potree generation tool'))

    parser.add_argument('s2pout', nargs='+',
                        help=('path(s) to the s2p output directory(ies)'))
    parser.add_argument('--outdir', metavar='potree_outdir', default='.',
                        help=('path to output directory'))
    args = parser.parse_args()

    try:
        produce_potree(args.s2pout, args.outdir)
    except common.RunFailure:
        basedir = os.path.dirname(os.path.abspath(__file__))
        print('You must download and compile PotreeConverter. Run the following commands:')
        print('    > cd %s'%basedir)
        print('    > git clone https://github.com/gfacciol/PotreeConverter_PLY_toolchain --recurse-submodules')
        print('    > cd PotreeConverter_PLY_toolchain')
        print('    > CC=gcc CXX=g++ make')
