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

# resolution of the output digital surface model, in meters per pixel
cfg['dsm_resolution'] = 4

# radius to compute altitudes (and to interpolate the small holes)
cfg['dsm_radius'] = 0

# dsm_sigma controls the spread of the blob from each point for the dsm computation
# (dsm_resolution by default)
cfg['dsm_sigma'] = None

# relative sift match threshold (else sift match threshold is absolute)
cfg['relative_sift_match_thresh'] = True

# if cfg['relative_sift_match_thresh'] is True :
# sift threshold on the first over second best match ratio
# else (absolute) a reasonable value is between 200 and 300 (128-vectors SIFT descriptors)
cfg['sift_match_thresh'] = 0.6

# disp range expansion facto
cfg['disp_range_extra_margin'] = 0.2

# register the rectified images with a shear estimated from the rpc data
cfg['register_with_shear'] = False

# number of ground control points per axis in matches from rpc generation
cfg['n_gcp_per_axis'] = 5

# max distance allowed for a point to the epipolar line of its match
cfg['epipolar_thresh'] = 0.5

# maximal pointing error, in pixels
cfg['max_pointing_error'] = 10

# triangulation mode : 'pairwise'or 'geometric'
cfg['triangulation_mode'] = 'pairwise'

# use global pointing for geometric triangulation
cfg['use_global_pointing_for_geometric_triangulation'] = False

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

cfg['rpc_alt_range_scale_factor'] = 1

# method to compute the disparity range: "sift", "exogenous", "wider_sift_exogenous", "fixed_pixel_range", "fixed_altitude_range"
cfg['disp_range_method'] = "wider_sift_exogenous"
cfg['disp_range_exogenous_low_margin'] = -10
cfg['disp_range_exogenous_high_margin'] = +100

# exogenous dem
cfg['exogenous_dem'] = None
cfg['exogenous_dem_geoid_mode'] = True

### stereo matching parameters

# stereo matching algorithm: 'tvl1', 'msmw', 'hirschmuller08',
# hirschmuller08_laplacian', 'sgbm', 'mgm', 'mgm_multi'
cfg['matching_algorithm'] = 'mgm'

# size of the Census NCC square windows used in mgm
cfg['census_ncc_win'] = 5

# MGM parameter: speckle filter minimum area (REMOVESMALLCC flag)
cfg['stereo_speckle_filter'] = 25

# MGM parameter: regularity (multiplies P1 and P2)
cfg['stereo_regularity_multiplier'] = 1.0

# clean height maps outliers
cfg['cargarse_basura'] = True

# longitude/latitude bounding box
cfg['ll_bbx'] = ("-inf", "inf", "-inf", "inf")
