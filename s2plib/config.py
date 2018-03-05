# Copyright (C) 2013, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2013, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2013, Julien Michel <julien.michel@cnes.fr>

# This module contains a dictionary, cfg, containing all the parameters of the
# s2p pipeline. This dictionary is updated at runtime with parameters defined
# by the user in the config.json file. All the optional parameters (that the
# user is not forced to define in config.json) must be defined here, otherwise
# they won't have a default value.

import os.path

cfg = {}

# path to output directory
cfg['out_dir'] = "s2p_output"

# path to directory where (many) temporary files will be stored
cfg['temporary_dir'] = "s2p_tmp"

# temporary files are erased when s2p terminates. Switch to False to keep them
cfg['clean_tmp'] = True

# remove all generated files except from ply point clouds and tif raster dsm
cfg['clean_intermediate'] = False

# switch to True if you want to process the whole image
cfg['full_img'] = False

# s2p processes the images tile by tile. The tiles are squares cropped from the
# reference image. The width and height of the tiles are given by this param, in pixels.
cfg['tile_size'] = 800

# margins used to increase the footprint of the rectified tiles, to
# account for poor disparity estimation close to the borders
cfg['horizontal_margin'] = 50  # for regularity and occlusions
cfg['vertical_margin'] = 10  # for regularity

# max number of processes launched in parallel. None means the number of available cores
cfg['max_processes'] = None

# max number of OMP threads used by programs compiled with openMP
cfg['omp_num_threads'] = 1

# debug mode (more verbose logs and intermediate results saved)
cfg['debug'] = False

# don't rerun when the output file is already there
cfg['skip_existing'] = False

# resolution of the output digital surface model, in meters per pixel
cfg['dsm_resolution'] = 4

# radius to compute altitudes (and to interpolate the small holes)
cfg['dsm_radius'] = 0

# dsm_sigma controls the spread of the blob from each point for the dsm computation
# (dsm_resolution by default)
cfg['dsm_sigma'] = None

# sift threshold on the first over second best match ratio
cfg['sift_match_thresh'] = 0.6

# disp range expansion facto
cfg['disp_range_extra_margin'] = 0.2

# register the rectified images with a shear estimated from the rpc data
cfg['register_with_shear'] = False

# number of ground control points per axis in matches from rpc generation
cfg['n_gcp_per_axis'] = 5

# max distance allowed for a point to the epipolar line of its match
cfg['epipolar_thresh'] = 0.5

# triangulation mode : 'pairwise'or 'geometric'
cfg['triangulation_mode'] = 'pairwise'

# use global pointing for geometric triangulation
cfg['use_global_pointing_for_geometric_triangulation'] = False

# stereo matching algorithm: 'tvl1', 'msmw', 'hirschmuller08',
# hirschmuller08_laplacian', 'sgbm', 'mgm'
cfg['matching_algorithm'] = 'mgm'

# size of the Census NCC square windows used in mgm
cfg['census_ncc_win'] = 5

# set these params if you want to impose the disparity range manually (cfg['disp_range_method'] == 'fixed_pixel_range')
cfg['disp_min'] = None
cfg['disp_max'] = None

# set these params if you want to impose the altitude range manually (cfg['disp_range_method'] == 'fixed_altitude_range')
cfg['alt_min'] = None
cfg['alt_max'] = None


# radius for erosion of valid disparity areas. Ignored if less than 2
cfg['msk_erosion'] = 2

cfg['fusion_operator'] = 'average_if_close'

# threshold (in meters) used for the fusion of two dems in triplet processing
# It should be adapted to the zoom factor
cfg['fusion_thresh'] = 3

cfg['disable_srtm'] = False
cfg['rpc_alt_range_scale_factor'] = 1

# method to compute the disparity range: "sift", "srtm", "wider_sift_srtm", "fixed_pixel_range", "fixed_altitude_range"
cfg['disp_range_method'] = "wider_sift_srtm"
cfg['disp_range_srtm_low_margin'] = -10
cfg['disp_range_srtm_high_margin'] = +100

# url of the srtm database mirror
cfg['srtm_url'] = 'http://138.231.80.250:443/srtm/tiff'
cfg['srtm_url'] = 'ftp://xftp.jrc.it/pub/srtmV4/tiff'
cfg['srtm_url'] = 'http://data_public:GDdci@data.cgiar-csi.org/srtm/tiles/GeoTIFF'

# directory where to store the srtm tiles
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
cfg['srtm_dir'] = os.path.join(parent_dir, '.srtm')

# clean height maps outliers
cfg['cargarse_basura'] = True

# longitude/latitude bounding box
cfg['ll_bbx'] = ("-inf", "inf", "-inf", "inf")

# use srtm to generate a watermask
cfg['use_srtm_for_water'] = False
