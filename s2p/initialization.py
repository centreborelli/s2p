# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import sys
import json
import copy
from typing import Any, List, Optional, Tuple
import rasterio
import numpy as np
import numpy.typing as npt
import rpcm

from s2p import common
from s2p import geographiclib
from s2p import rpc_utils
from s2p import masking
from s2p import parallel
from s2p.config import cfg
from s2p.tile import Tile

# This function is here as a workaround to python bug #24313 When
# using python3, json does not know how to serialize numpy.int64 on
# some platform numpy also decides to go for int64 when numpy.arange
# is called. This results in our json not being serializable anymore
# Calling json.dump(..,default=workaround_json_int64) fixes this
# https://bugs.python.org/issue24313
def workaround_json_int64(o: Any) -> int:
    if isinstance(o, np.integer) : return int(o)
    raise TypeError


def dict_has_keys(d: dict, l: List[str]) -> bool:
    """
    Return True if the dict d contains all the keys of the input list l.
    """
    return all(k in d for k in l)


def check_parameters(d: dict) -> None:
    """
    Check that the provided dictionary defines all mandatory s2p arguments.

    Args:
        d: python dictionary
    """
    # verify that input files paths are defined
    if 'images' not in d or len(d['images']) < 2:
        print('ERROR: missing paths to input images')
        sys.exit(1)
    for img in d['images']:
        if not dict_has_keys(img, ['img']):
            print('ERROR: missing img paths for image', img)
            sys.exit(1)

    # read RPCs
    for img in d['images']:
        if 'rpc' in img:
            if isinstance(img['rpc'], str):  # path to an RPC file
                img['rpcm'] = rpcm.rpc_from_rpc_file(img['rpc'])
            elif isinstance(img['rpc'], dict):  # RPC dict in 'rpcm' format
                img['rpcm'] = rpcm.RPCModel(img['rpc'], dict_format='rpcm')
            else:
                raise NotImplementedError(
                    'rpc of type {} not supported'.format(type(img['rpc']))
                )
        else:
            img['rpcm'] = rpcm.rpc_from_geotiff(img['img'])

    # verify that an input ROI is defined
    if d.get("full_img"):
        with rasterio.open(d['images'][0]['img'], "r") as f:
            width = f.width
            height = f.height
        d['roi'] = {'x': 0, 'y': 0, 'w': width, 'h': height}
    elif 'roi' in d and dict_has_keys(d['roi'], ['x', 'y', 'w', 'h']):
        pass
    elif 'roi_geojson' in d:
        ll_poly = geographiclib.read_lon_lat_poly_from_geojson(d['roi_geojson'])
        use_srtm = d.get('use_srtm', False)
        exogenous_dem = d.get('exogenous_dem')
        exogenous_dem_geoid_mode = d.get('exogenous_dem_geoid_mode', False)
        d['roi'] = rpc_utils.roi_process(d['images'][0]['rpcm'], ll_poly,
                                         use_srtm=use_srtm,
                                         exogenous_dem=exogenous_dem,
                                         exogenous_dem_geoid_mode=exogenous_dem_geoid_mode)
    else:
        print('ERROR: missing or incomplete roi definition')
        sys.exit(1)

    # d['roi'] : all the values must be integers
    d['roi']['x'] = int(np.floor(d['roi']['x']))
    d['roi']['y'] = int(np.floor(d['roi']['y']))
    d['roi']['w'] = int(np.ceil(d['roi']['w']))
    d['roi']['h'] = int(np.ceil(d['roi']['h']))

    # warn about unknown parameters. The known parameters are those defined in
    # the global config.cfg dictionary, plus the mandatory 'images' and 'roi'
    for k in d.keys():
        if k not in ['images', 'roi', 'roi_geojson']:
            if k not in cfg:
                print('WARNING: ignoring unknown parameter {}.'.format(k))


def build_cfg(user_cfg: dict) -> None:
    """
    Populate a dictionary containing the s2p parameters from a user config file.

    This dictionary is contained in the global variable 'cfg' of the config
    module.

    Args:
        user_cfg: user config dictionary
    """
    # check that all the mandatory arguments are defined
    check_parameters(user_cfg)

    # fill the config module: updates the content of the config.cfg dictionary
    # with the content of the user_cfg dictionary
    cfg.update(user_cfg)

    # set keys 'clr', 'cld' and 'roi' of the reference image to None if they
    # are not already defined. The default values of these optional arguments
    # can not be defined directly in the config.py module. They would be
    # overwritten by the previous update, because they are in a nested dict.
    cfg['images'][0].setdefault('clr')
    cfg['images'][0].setdefault('cld')
    cfg['images'][0].setdefault('roi')
    cfg['images'][0].setdefault('wat')

    # make sure that input data have absolute paths
    for i in range(len(cfg['images'])):
        for d in ['clr', 'cld', 'roi', 'wat', 'img']:
            if d in cfg['images'][i] and cfg['images'][i][d] is not None and not os.path.isabs(cfg['images'][i][d]):
                cfg['images'][i][d] = os.path.abspath(cfg['images'][i][d])

    # get out_crs
    if 'out_crs' not in cfg or cfg['out_crs'] is None:
        x, y, w, h = [cfg['roi'][k] for k in ['x', 'y', 'w', 'h']]
        utm_zone = rpc_utils.utm_zone(cfg['images'][0]['rpcm'], x, y, w, h)
        epsg_code = geographiclib.epsg_code_from_utm_zone(utm_zone)
        cfg['out_crs'] = "epsg:{}".format(epsg_code)
        if cfg['out_geoid']:
            # Use the EGM96 geoid model for the output CRS if out_geoid is True
            cfg['out_crs'] += "+5773"
    geographiclib.pyproj_crs(cfg['out_crs'])

    # get image ground sampling distance
    cfg['gsd'] = rpc_utils.gsd_from_rpc(cfg['images'][0]['rpcm'])


def make_dirs() -> None:
    """
    Create directories needed to run s2p.
    """
    os.makedirs(cfg['out_dir'], exist_ok=True)
    os.makedirs(os.path.expandvars(cfg['temporary_dir']), exist_ok=True)

    # store a json dump of the config.cfg dictionary
    with open(os.path.join(cfg['out_dir'], 'config.json'), 'w') as f:
        cfg_copy = copy.deepcopy(cfg)
        cfg_copy['out_dir'] = '.'
        for img in cfg_copy['images']:
            img.pop('rpcm', None)
        json.dump(cfg_copy, f, indent=2, default=workaround_json_int64)


def adjust_tile_size() -> Tuple[int, int]:
    """
    Adjust the size of the tiles.
    """

    tile_w = min(cfg['roi']['w'], cfg['tile_size'])  # tile width
    ntx = int(np.round(float(cfg['roi']['w']) / tile_w))
    # ceil so that, if needed, the last tile is slightly smaller
    tile_w = int(np.ceil(float(cfg['roi']['w']) / ntx))

    tile_h = min(cfg['roi']['h'], cfg['tile_size'])  # tile height
    nty = int(np.round(float(cfg['roi']['h']) / tile_h))
    tile_h = int(np.ceil(float(cfg['roi']['h']) / nty))

    print('tile size: {} {}'.format(tile_w, tile_h))
    n = len(cfg['images'])
    if n == 2:
        print('total number of tiles: {} ({} x {})'.format(ntx * nty, ntx, nty))
    else:
        print('total number of tiles: {} ({} x {}) x {} pairs'.format(ntx*nty*(n-1),
                                                                      ntx, nty, n-1))
    return tile_w, tile_h


def compute_tiles_coordinates(
    rx: int, ry: int, rw: int, rh: int, tw: int, th: int
) -> Tuple[List[Tuple[int, int, int, int]], dict]:
    out = []
    neighborhood_dict = dict()

    for y in np.arange(ry, ry + rh, th):
        h = min(th, ry + rh - y)
        for x in np.arange(rx, rx + rw, tw):
            w = min(tw, rx + rw - x)

            out.append((x, y, w, h))

            # get coordinates of tiles from neighborhood
            out2 = []
            for y2 in [y - th, y, y + th]:
                h2 = min(th, ry + rh - y2)
                for x2 in [x - tw, x, x + tw]:
                    w2 = min(tw, rx + rw - x2)
                    if rx + rw > x2 >= rx:
                        if ry + rh > y2 >= ry:
                            out2.append((x2, y2, w2, h2))

            neighborhood_dict[str((x, y, w, h))] = out2

    return out, neighborhood_dict


def get_tile_dir(x: int, y: int, w: int, h: int) -> str:
    """
    Get the name of a tile directory
    """
    return os.path.join('tiles','row_{:07d}_height_{}'.format(y, h),
                        'col_{:07d}_width_{}'.format(x, w))


def create_tile(
    coords: Tuple[int, int, int, int], neighborhood_coords_dict: dict
) -> Tile:
    """
    Return a dictionary with the data of a tile.

    Args:
        coords (tuple): 4-tuple of ints giving the x, y, w, h coordinates of a
            tile, where x, y are the top-left corner coordinates and w, h the
            width and height
        neighborhood_coords_dict (dict): dictionary with the list of
            neighboring tiles of each tile. The keys of this dict are string
            identifying the tiles, and the values are lists of tuples of
            coordinates of neighboring tiles coordinates

    Returns:
        tile (dict): dictionary with the metadata of a tile
    """
    dir = os.path.join(cfg['out_dir'], get_tile_dir(*coords))
    json = os.path.join(get_tile_dir(*coords), 'config.json')

    neighborhood_dirs = list()
    key = str(coords)
    if 'neighborhood_dirs' in cfg:
        neighborhood_dirs = cfg['neighborhood_dirs']
    elif key in neighborhood_coords_dict:
        for coords2 in neighborhood_coords_dict[key]:
            neighborhood_dirs.append(os.path.join('../../..',
                                                  get_tile_dir(*coords2)))

    return Tile(coordinates=coords, dir=dir, json=json, neighborhood_dirs=neighborhood_dirs)


def rectangles_intersect(
    r: Tuple[int, int, int, int], s: Tuple[int, int, int, int]
) -> bool:
    """
    Check intersection of two rectangles parallel to the coordinate axis.

    Args:
        r (tuple): 4 floats that define the coordinates of the top-left corner,
            the width and the height of a rectangle
        s (tuple): 4 floats that define the coordinates of the top-left corner,
            the width and the height of a rectangle

    Return:
        bool telling if the rectangles intersect
    """
    rx, ry, rw, rh = r
    sx, sy, sw, sh = s

    # check if one rectangle is entirely above the other
    if ry + rh < sy or sy + sh < ry:
        return False

    # check if one rectangle is entirely left of the other
    if rx + rw < sx or sx + sw < rx:
        return False

    return True


def is_tile_all_nodata(path: str, window: rasterio.windows.Window) -> bool:
    """Check if pixels in a given window are all nodata.

    Parameters
    ----------
    path
        Path to the raster.
    window
        A rasterio.windows.Window object.

    Returns
    -------
        Return True if all pixels in the window are nodata.
        Return False if at least one pixel is non-nodata.
    """
    with rasterio.open(path, "r") as ds:
        arr = ds.read(window=window)

        # NOTE: Many satellite imagery providers use ds.nodata as the value of
        # nodata pixels. Pleiades and PNeo imagery use None as nodata in their
        # profile while putting 0 to nodata pixel in reality. Thus, we have to
        # check both ds.nodata and 0 here. I.e., if a window is full of nodata
        # or 0, then this window is discarded.
        if np.all(arr == 0) or np.all(arr == ds.nodata):
            return True
        else:
            return False


def is_this_tile_useful(
    x: int, y: int, w: int, h: int, images_sizes: List[Tuple[int, int]]
) -> Tuple[bool, Optional[npt.NDArray[np.bool_]]]:
    """
    Check if a tile contains valid pixels.

    Valid pixels must be found in the reference image plus at least one other image.

    Args:
        x, y, w, h (ints): 4 ints that define the coordinates of the top-left corner,
            the width and the height of a rectangular tile
        images_sizes (list): list of tuples with the height and width of the images

    Return:
        useful (bool): bool telling if the tile has to be processed
        mask (np.array): tile validity mask. Set to None if the tile is discarded
    """
    if is_tile_all_nodata(cfg["images"][0]["img"], rasterio.windows.Window(x, y, w, h)):
        return False, None
        
    # check if the tile is partly contained in at least one other image
    rpc = cfg['images'][0]['rpcm']
    for img, size in zip(cfg['images'][1:], images_sizes[1:]):
        coords = rpc_utils.corresponding_roi(cfg, rpc, img['rpcm'], x, y, w, h)
        if rectangles_intersect(coords, (0, 0, size[1], size[0])):
            break  # the tile is partly contained
    else:  # we've reached the end of the loop hence the tile is not contained
        return False, None

    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    mask = masking.image_tile_mask(x, y, w, h, roi_msk, cld_msk, wat_msk,
                                   images_sizes[0], cfg['border_margin'])
    if not mask.any():
        return False, None
    return True, mask


def tiles_full_info(tw: int, th: int, tiles_txt: str, create_masks: bool = False) -> List[Tile]:
    """
    List the tiles to process and prepare their output directories structures.

    Most of the time is spent discarding tiles that are masked by water
    (according to exogenous dem).

    Returns:
        a list of dictionaries. Each dictionary contains the image coordinates
        and the output directory path of a tile.
    """
    rpc = cfg['images'][0]['rpcm']
    roi_msk = cfg['images'][0]['roi']
    cld_msk = cfg['images'][0]['cld']
    wat_msk = cfg['images'][0]['wat']
    rx = cfg['roi']['x']
    ry = cfg['roi']['y']
    rw = cfg['roi']['w']
    rh = cfg['roi']['h']

    # list of dictionaries (one for each non-masked tile)
    tiles = []

    # list tiles coordinates
    tiles_coords, neighborhood_coords_dict = compute_tiles_coordinates(rx, ry, rw, rh, tw, th)

    if create_masks or not os.path.exists(tiles_txt):
        print('\ndiscarding masked tiles...')
        images_sizes = []
        for img in cfg['images']:
            with rasterio.open(img['img'], 'r') as f:
                images_sizes.append(f.shape)

        # compute all masks in parallel as numpy arrays
        tiles_usefulnesses = parallel.launch_calls(is_this_tile_useful,
                                                   tiles_coords,
                                                   cfg['max_processes'],
                                                   images_sizes,
                                                   tilewise=False,
                                                   timeout=cfg['timeout'])

        # discard useless tiles from neighborhood_coords_dict
        discarded_tiles = set(x for x, (b, _) in zip(tiles_coords, tiles_usefulnesses) if not b)
        for k, v in neighborhood_coords_dict.items():
            neighborhood_coords_dict[k] = list(set(v) - discarded_tiles)

        for coords, usefulness in zip(tiles_coords, tiles_usefulnesses):

            useful, mask = usefulness
            if not useful:
                continue

            tile = create_tile(coords, neighborhood_coords_dict)
            tiles.append(tile)

            # make tiles directories and store json configuration dumps
            os.makedirs(tile.dir, exist_ok=True)
            for i in range(1, len(cfg['images'])):
                os.makedirs(os.path.join(tile.dir, 'pair_{}'.format(i)), exist_ok=True)

            # save a json dump of the tile configuration
            tile_cfg = copy.deepcopy(cfg)
            x, y, w, h = tile.coordinates
            for img in tile_cfg['images']:
                img.pop('rpcm', None)
            tile_cfg['roi'] = {'x': x, 'y': y, 'w': w, 'h': h}
            tile_cfg['full_img'] = False
            tile_cfg['max_processes'] = 1
            tile_cfg['neighborhood_dirs'] = tile.neighborhood_dirs
            tile_cfg['out_dir'] = '../../..'

            with open(os.path.join(cfg['out_dir'], tile.json), 'w') as f:
                json.dump(tile_cfg, f, indent=2, default=workaround_json_int64)

            # save the mask
            common.rasterio_write(os.path.join(tile.dir, 'mask.png'),
                                  mask.astype(np.uint8))
    else:
        if len(tiles_coords) == 1:
            tiles.append(create_tile(tiles_coords[0], neighborhood_coords_dict))
        else:
            with open(tiles_txt, 'r') as f_tiles:
                for config_json in f_tiles:
                    with open(os.path.join(cfg['out_dir'],
                                           config_json.rstrip(os.linesep)), 'r') as f_config:
                        tile_cfg = json.load(f_config)
                        roi = tile_cfg['roi']
                        coords = roi['x'], roi['y'], roi['w'], roi['h']
                        tiles.append(create_tile(coords, neighborhood_coords_dict))

    return tiles