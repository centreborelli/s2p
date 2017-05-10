#!/usr/bin/env python
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>
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

def vrt_body_source(fname,band,src_x,src_y,src_w,src_h,dst_x,dst_y,dst_w,dst_h):
    """
    Generate a source section in vrt body.
    
    Args:
        fname: Relative path to the source image
        band: index of the band to use as source
        src_x, src_y, src_w, src_h: source window (cropped from source image)
        dst_x, dst_y, dst_w, dst_h: destination window (where crop will be pasted)
    """
    
    body ='\t\t<SimpleSource>\n'
    body+='\t\t\t<SourceFileName relativeToVRT=\'1\'>%s</SourceFileName>\n' % fname
    body+='\t\t\t<SourceBand>%i</SourceBand>\n' % band
    body+='\t\t\t<SrcRect xOff=\'%i\' yOff=\'%i\'' % (src_x, src_y)
    body+='xSize=\'%i\' ySize=\'%i\'/>\n' % (src_w, src_h)
    body+='\t\t\t<DstRect xOff=\'%i\' yOff=\'%i\'' % (dst_x, dst_y)
    body+='xSize=\'%i\' ySize=\'%i\'/>\n'  %(dst_w, dst_h)
    body+='\t\t</SimpleSource>\n'

    return body

def vrt_header(w,h,dataType='Float32'):
    """
    Generate vrt header for a monoband image
    
    Args:
        w,h: size of the corresponding raster
        dataType: Type of the raster (default is Float32)
    
    """
    header = '<VRTDataset rasterXSize=\"%i\" rasterYSize=\"%i\">\n' %(w,h)
    header+= '\t<VRTRasterBand dataType=\"%s\" band=\"1\">\n' %(dataType)
    header+= '\t\t<ColorInterp>Gray</ColorInterp>\n'

    return header

def vrt_footer():
    """
    Generate vrt footer
    """
    footer = '\t</VRTRasterBand>\n'
    footer+= '</VRTDataset>\n'

    return footer

def global_extent(tiles):
    """
    Compute the global raster extent from a list of tiles
    Args:
        tiles: list of config files loaded from json files
    Returns:
         (min_x,max_x,min_y,max_y) tuple
    """
    min_x = None
    max_x = None
    min_y = None
    max_y = None

    # First loop is to compute global extent 
    for tile in tiles:
        with open(tile,'r') as f:

            tile_cfg = json.load(f)

            x = tile_cfg['roi']['x']
            y = tile_cfg['roi']['y']
            w = tile_cfg['roi']['w']
            h = tile_cfg['roi']['h']
            
            if min_x is None or x < min_x:
                min_x = x
            if min_y is None or y < min_y:
                min_y = y
            if max_x is None or x + w > max_x:
                max_x = x + w
            if max_y is None or y + h > max_y:
                max_y = y + h
                
    return(min_x,max_x,min_y,max_y)

def write_row_vrts(tiles,sub_img,vrt_basename,min_x,max_x):
    """
    Write intermediate vrts (one per row)
    
    Args:
        tiles: list of config files loaded from json files
        sub_img: Relative path of the sub-image to mosaic (for ex. height_map.tif)
        vrt_basename: basename of the output vrt 
        min_x, max_x: col extent of the raster
    Returns:
        A dictionnary of vrt files with sections vrt_body, th and vrt_dir
    """
    vrt_row = {}

    # First loop, write all row vrts body section
    for tile in tiles:
        with open(tile,'r') as f:
            
            tile_cfg = json.load(f)
            
            x = tile_cfg['roi']['x']
            y = tile_cfg['roi']['y']
            w = tile_cfg['roi']['w']
            h = tile_cfg['roi']['h']
            
            tile_dir = os.path.dirname(tile)
            row_vrt_dir = os.path.dirname(tile_dir)
            tile_sub_img_dir = os.path.basename(tile_dir)

            vrt_row.setdefault(y,{'vrt_body':"",'vrt_dir':row_vrt_dir,"th": h})
            
            tile_sub_img = os.path.join(tile_sub_img_dir,sub_img)

            # Check if source image exists
            if not os.path.exists(os.path.join(row_vrt_dir,tile_sub_img)):
                print('Warning: '+tile_sub_img+' does not exist, skipping ...')
                continue
            
            vrt_row[y]['vrt_body']+=vrt_body_source(tile_sub_img,1,0,0,w,h,
                                                    x-min_x,0,w,h)

    # Second loop, write all row vrts
    # Do not use items()/iteritems() here because of python 2 and 3 compat
    for y in vrt_row:
        vrt_data = vrt_row[y]
        row_vrt_filename = os.path.join(vrt_data['vrt_dir'],vrt_basename)
        
        with  open(row_vrt_filename,'w') as row_vrt_file:
            # Write vrt header
            row_vrt_file.write(vrt_header(max_x-min_x,vrt_data['th']))

            # Write vrt body
            row_vrt_file.write(vrt_data['vrt_body'])

            # Write vrt footer
            row_vrt_file.write(vrt_footer())
            
    return vrt_row

def write_main_vrt(vrt_row,vrt_name,min_x,max_x,min_y,max_y):
    """
    Write the main vrt file
    
    Args:
        vrt_row: The vrt files dictionnary from write_row_vrts()
        vrt_name: The output vrt_name
        min_x,max_x,min_y,max_y: Extent of the raster 
    """
    vrt_basename = os.path.basename(vrt_name)
    vrt_dirname = os.path.dirname(vrt_name)
    
    with open(vrt_name,'w') as main_vrt_file:
    
        main_vrt_file.write(vrt_header(max_x-min_x,max_y-min_y))
        # Do not use items()/iteritems() here because of python 2 and 3 compat
        for y in vrt_row:
            vrt_data = vrt_row[y]
            relative_vrt_dir = os.path.relpath(vrt_data['vrt_dir'],vrt_dirname)
            row_vrt_filename = os.path.join(relative_vrt_dir,vrt_basename)

            vrt_body_src=vrt_body_source(row_vrt_filename,1,0,0,max_x-min_x,
                                         vrt_data['th'],0,y-min_y,max_x-min_x,
                                         vrt_data['th'])
            
            main_vrt_file.write(vrt_body_src)
                
        main_vrt_file.write(vrt_footer())

            
def main(tiles_file,outfile,sub_img):

    outfile_basename = os.path.basename(outfile)
    outfile_dirname  = os.path.dirname(outfile)
    
    output_format = outfile_basename[-3:]

    print('Output format is '+output_format)

    # If output format is tif, we need to generate a temporary vrt
    # with the same name
    vrt_basename = outfile_basename
    
    if output_format == 'tif':
        vrt_basename = vrt_basename[:-3]+'vrt'
    elif output_format !='vrt':
        print('Error: only vrt or tif extension is allowed for output image.')
        return
    
    vrt_name = os.path.join(outfile_dirname,vrt_basename)
    
    # Read the tiles file
    tiles = s2p.read_tiles(tiles_file)

    print(str(len(tiles))+' tiles found')

    # Compute the global extent of the output image
    (min_x,max_x,min_y,max_y) = global_extent(tiles)
    
    print('Global extent: [%i,%i]x[%i,%i]'%(min_x,max_x,min_y,max_y))

    # Now, write all row vrts
    print("Writing row vrt files "+vrt_basename)
    vrt_row = write_row_vrts(tiles,sub_img,vrt_basename,min_x,max_x)
    
    # Finally, write main vrt
    print('Writing '+vrt_name)
    write_main_vrt(vrt_row,vrt_name,min_x,max_x,min_y,max_y)

    # If Output format is tif, convert vrt file to tif
    if output_format == 'tif':
        print('Converting vrt to tif ...')
        common.run(('gdal_translate -ot Float32 -co TILED=YES -co'
                    ' BIGTIFF=IF_NEEDED %s %s'
                    %(common.shellquote(vrt_name),common.shellquote(outfile))))

        print('Removing temporary vrt files')
        # Do not use items()/iteritems() here because of python 2 and 3 compat
        for y in vrt_row:
            vrt_data = vrt_row[y]
            row_vrt_filename = os.path.join(vrt_data['vrt_dir'],vrt_basename)
            try:
                os.remove(row_vrt_filename)
            except OSError:
                pass
        try:
            os.remove(vrt_name)
        except OSError:
            pass
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: mosaic tool'))
    
    parser.add_argument('tiles',metavar='tiles.txt',
                        help=('path to the tiles.txt file'))
    parser.add_argument('outfile',metavar='out.tif',
                        help=('path to the output file.'
                              ' File extension can be .tif or .vrt'))
    parser.add_argument('sub_img',metavar='pair_1/height_map.tif',
                        help=('path to the sub-image to mosaic.'
                              ' Can be (but not limited to) height_map.tif,'
                              ' pair_n/height_map.tif, pair_n/rpc_err.tif,'
                              ' cloud_water_image_domain_mask.png.'
                              ' Note that rectified_* files CAN NOT be used.'))
    args = parser.parse_args()

    main(args.tiles,args.outfile,args.sub_img)
