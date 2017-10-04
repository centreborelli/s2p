#!/usr/bin/env python

# Copyright (C) 2017, Julien Michel <julien.michel@cnes.fr>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published                                                                                                                                        # by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import tempfile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from s2plib import rpc_model
from s2plib import common

available_filters = ['near', 'bilinear',
                     'cubic', 'cubicspline', 'lanczos', 'average']

if __name__ == "__main__":

    if len(sys.argv) not in [7, 8]:
        print("Usage: {} im_in rpc_in im_out rpc_out filter scale_x [scale_y]".format(
            sys.argv[0]))
        print("Use scale_x (resp. scale_y) > 1 to zoom in, and scale_x (resp. scale_y) < 1 to zoom out")
        print("Filter should be one of the following gdal filters: {}".format(
            available_filters))
        sys.exit(1)
    in_img_file = sys.argv[1]
    in_rpc_file = sys.argv[2]
    out_img_file = sys.argv[3]
    out_rcp_file = sys.argv[4]
    filt = sys.argv[5]
    scale_x = float(sys.argv[6])

    # handle optional scale_y parameter
    if len(sys.argv) == 7:
        scale_y = scale_x
    else:
        scale_y = float(sys.argv[7])

    # check if filter is valid
    if filt not in available_filters:
        print("Unknown filter {}. Should be one of {}".format(
            filt, available_filters))
        sys.exit(1)

    # generate image
    print("Generating {} ...".format(out_img_file))

    # First get input image size
    sz = common.image_size_gdal(in_img_file)
    w = sz[0]
    h = sz[1]

    # Generate a temporary vrt file to have the proper geotransform
    fd, tmp_vrt = tempfile.mkstemp(suffix='.vrt',
                                   dir=os.path.dirname(out_img_file))

    os.close(fd)

    common.run('gdal_translate -of VRT -a_ullr 0 0 %d %d %s %s' %
               (w, h, in_img_file, tmp_vrt))

    common.run(('gdalwarp -co RPB=NO -co PROFILE=GeoTIFF -r %s -co "BIGTIFF=IF_NEEDED" -co "TILED=YES" -ovr NONE -overwrite -to SRC_METHOD=NO_GEOTRANSFORM -to DST_METHOD=NO_GEOTRANSFORM -tr'
                ' %d %d %s %s') % (filt, scale_x, scale_y, tmp_vrt, out_img_file))

    try:
        # Remove aux files if any
        os.remove(out_img_file + ".aux.xml")
    except OSError:
        pass

    # Clean tmp vrt file
    os.remove(tmp_vrt)

    print("Done")

    # generate rpc file
    print("Generating {} ...".format(out_rcp_file))

    r = rpc_model.RPCModel(in_rpc_file)
    r.linScale /= scale_y
    r.linOff /= scale_y
    r.colScale /= scale_x
    r.colOff /= scale_x
    r.write(out_rcp_file)

    print("Done")

    sys.exit(0)
