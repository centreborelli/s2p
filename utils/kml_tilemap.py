#!/usr/bin/env python

# Copyright (C) 2015, David Youssefi <david.youssefi@cnes.fr>

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

import os
import argparse
import sys
import re
import subprocess
import json
import collections
import utm
import simplekml
import gdal
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import s2p
from s2plib import rpc_model
from s2plib import common

def pix_2_latlon(gt, px, py, zone_number, northern):
    x = px * gt[1] + gt[0]
    y = py * gt[5] + gt[3]

    if zone_number is not None:
        lat, lon = utm.to_latlon(x, y, zone_number, northern=northern)
    else:
        lat, lon = y, x

    return lon, lat, 0

def read_tiles(tile_files, outdir, m, M, key):
    tiles = []
    tile_file_dir = os.path.dirname(tile_files)
    err_log = os.path.join(outdir, "%s_invalid_tiles.txt" % key)
    kml = simplekml.Kml()

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

                if message != "ok":
                    ferr.write(t)
                    ferr.write('\n')
                    write_tiles_polygon(t, kml, m, M, message, True)
                    continue
                if message == "ok":
                    tiles.append(t)

    kml.save(os.path.join(outdir, "%s_error.kml" % key))

    return tiles

def get_coordinates_with_img(img):
    ds = gdal.Open(img)
    gt = ds.GetGeoTransform()

    prj = ds.GetProjection()
    utm_zone = re.findall("UTM zone (\w*)", prj)
    if utm_zone != list():
        zone_number = int(utm_zone[0][:-1])
        northern = (utm_zone[0][-1] == 'N')
    else:
        zone_number = None
        northern = None

    size_x = ds.RasterXSize
    size_y = ds.RasterYSize

    roi = [[0, size_y],
           [size_x, size_y],
           [size_x, 0],
           [0, 0],
           [0, size_y]]

    return [pix_2_latlon(gt, x[0], x[1], zone_number, northern) for x in roi]

def get_coordinates_with_config(tile, m, M):
    tile_cfg = s2p.read_config_file(os.path.join(tile, "config.json"))

    x = tile_cfg['roi']['x']
    y = tile_cfg['roi']['y']
    w = tile_cfg['roi']['w']
    h = tile_cfg['roi']['h']

    rpcfile = tile_cfg['images'][0]['rpc']
    rpc = rpc_model.RPCModel(rpcfile)

    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    c = np.array([m, M,   m,   M,   m,   M,   m,   M])

    lon, lat, __ = rpc.direct_estimate(a, b, c)

    out = list(common.bounding_box2D(np.vstack([lon, lat]).T))

    out[2] += out[0]
    out[3] += out[1]

    latlon = [[out[0], out[3], 0],
              [out[2], out[3], 0],
              [out[2], out[1], 0],
              [out[0], out[1], 0],
              [out[0], out[3], 0]]

    return latlon

def gdal_to_ground_overlay(img, name, kml, outdir, max_h, min_h):
    latlon = get_coordinates_with_img(img)
    tmp = os.path.join("href", name+".tif")
    preview = os.path.join("href", name+".png")
    cmd = ["gdal_translate", "-scale", str(min_h), str(max_h), "1", "255",
           "-ot", "Byte", img, os.path.join(outdir, tmp)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    q = p.communicate()

    cmd = ["convert", os.path.join(outdir, tmp),
           "-transparent", "black", os.path.join(outdir, preview)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    q = p.communicate()

    os.remove(os.path.join(outdir, tmp))

    ground = kml.newgroundoverlay(name=name)
    ground.icon.href = preview
    ground.gxlatlonquad.coords = latlon
    return ground

def get_polygon_description(dico):
    descr = ['<![CDATA[',
             '<style> table { border-collapse: collapse; width: 100%;} th, td {text-align: left; padding: 8px; } tr:nth-child(even){background-color: #f2f2f2}} </style>',
             '<table border="1">']

    for key in dico:
        if dico[key]["style"] is not None:
            descr += ['<tr style="%s">' % dico[key]["style"]]
        else:
            descr += ['<tr>']

        descr += ['<th>%s</th>' % key,
                  '<th>%s</th>' % dico[key]["value"],
                  '</tr>']
    descr += ['</table>]]>\n']

    return "".join(descr)

def write_tiles_polygon(tile, kml, m=None, M=None, message=None, error_mode=False):
    dsm = os.path.join(tile, 'dsm.tif')
    dico = []

    green_style = "background-color: #4CAF50; color: white;"
    red_style = "background-color: #FF4C4C; color: white;"
    blue_style = "background-color: #9999EB; color:white;"

    if error_mode is True:
        latlon = get_coordinates_with_config(tile, m, M)
        color = simplekml.Color.red
        head_style = red_style
    else:
        latlon = get_coordinates_with_img(dsm)
        color = simplekml.Color.green
        head_style = green_style

    tile_cfg = s2p.read_config_file(os.path.join(tile, "config.json"))
    x = tile_cfg['roi']['x']
    y = tile_cfg['roi']['y']
    w = tile_cfg['roi']['w']
    h = tile_cfg['roi']['h']

    pair_info = [{} for i in range(len(tile_cfg['images'])-1)]

    for i in range(len(tile_cfg['images'])-1):
        pair_dir = os.path.join(tile,
                                "pair_%d" % (i+1))
        disp = os.path.join(pair_dir,
                            "rectified_disp.tif")

        if os.path.exists(disp) is False:
            pair_info[i]['status'] = "failure"
            if error_mode is False:
                color = simplekml.Color.blue
                pair_info[i]['style'] = red_style
                head_style = blue_style
            else:
                pair_info[i]['style'] = None
        else:
            pair_info[i]['status'] = "success"
            pair_info[i]['style'] = None

        disp_min_max = os.path.join(pair_dir,
                                    "disp_min_max.txt")

        if os.path.exists(disp_min_max) is True:
            disp_min, disp_max = np.loadtxt(disp_min_max)
        else:
            disp_min, disp_max = None, None

        pair_info[i]['disp_min'] = disp_min
        pair_info[i]['disp_max'] = disp_max

    dico += [('roi', {"value"  : "x : %s - y : %s - w : %s - h : %s" % (x, y, w, h),
                      "style" : head_style})]

    for i in range(len(tile_cfg['images'])-1):
        disp_min = pair_info[i]['disp_min']
        disp_max = pair_info[i]['disp_max']
        status = pair_info[i]['status']
        style = pair_info[i]['style']
        dico += [('pair_%d' % (i+1), {"value":' - '.join(['disp_min : %s' % disp_min,
                                                          'disp_max : %s' % disp_max,
                                                          'status : %s' % status]),
                                      "style":style})]

    if message != None:
        dico += [('message', {"value": message,
                              "style": None})]
    
    pol = kml.newpolygon(name=tile)
    pol.outerboundaryis = latlon
    pol.style.linestyle.color = color
    pol.style.linestyle.width = 5
    pol.style.polystyle.color = simplekml.Color.changealphaint(100, color)

    dico = collections.OrderedDict(dico)
    pol.description = get_polygon_description(dico)
    

def get_min_max(im):
    ds = gdal.Open(im)

    print "Compute %s statistics..." % im
    try:
        stats = ds.GetRasterBand(1).GetStatistics(0, 1)
        M = stats[2] + stats[3]
        m = stats[2] - stats[3]
    except RuntimeError:
        print "no valid pixels found in sampling"
        return 0, 0
    return m, M

def write_overlay(tiles, outdir, m, M, key):
    print ("Rescale the input pixels values")
    print ("from the range %f to %f to the range 0 to 255" % (m,
                                                              M))
    print ("Create ground overlay...")
    kml = simplekml.Kml()
    add = 1
    add_tot = len(tiles)
    for tile in tiles:
        sys.stdout.write("\r%d/%d"% (add, add_tot))
        sys.stdout.flush()
        add += 1
        dsm = os.path.join(tile, 'dsm.tif')
        gdal_to_ground_overlay(dsm, "_".join(dsm.split(os.sep)[-3:-1]),
                               kml, outdir, M, m)
    kml.save(os.path.join(outdir, "%s_ground_overlay.kml" % key))
    print
    print "kml saved"

def write_tiles_info(tiles, outdir, key):
    print ("Create tiles info...")
    add = 1
    add_tot = len(tiles)

    kml = simplekml.Kml()
    for tile in tiles:
        sys.stdout.write("\r%d/%d"% (add, add_tot))
        sys.stdout.flush()
        add += 1
        write_tiles_polygon(tile, kml)

    kml.save(os.path.join(outdir, "%s_info.kml" % key))
    print
    print "kml saved"


def main(tiles_file, outdir, key, with_overlay=False):
    # makedirs
    try:
        os.makedirs(os.path.join(outdir, "href"))
    except OSError:
        pass

    # get min max
    final_dsm = os.path.join(os.path.dirname(tiles_file), "dsm.tif")
    if os.path.exists(final_dsm) is True:
        m, M = get_min_max(final_dsm)
    else:
        m, M = 0, 0
    print "min : %s, max : %s" % (m, M)

    # Read the tiles file
    tiles = read_tiles(tiles_file, outdir, m, M, key)
    add_tot = len(tiles)
    print (str(add_tot)+' tiles found')

    # overlay
    if with_overlay is True:
        write_overlay(tiles, outdir, m, M, key)

    # tiles info
    write_tiles_info(tiles, outdir, key)

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=('S2P: kml tilemap'))
    PARSER.add_argument('--tiles', required=True,
                        help=('path to the tiles.txt file'))
    PARSER.add_argument('--outdir', required=True,
                        help=('path to the output directory.'))
    PARSER.add_argument('--ID', default='dsm',
                        help=('basename for output files'))
    PARSER.add_argument('--overlay', action='store_true', help=('compute overlay'))

    ARGS = PARSER.parse_args()

    main(ARGS.tiles, ARGS.outdir, ARGS.ID, ARGS.overlay)
