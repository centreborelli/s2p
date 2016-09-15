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
cfg['out_dir'] = "output"

# path to directory where (many) temporary files will be stored
cfg['temporary_dir'] = "/tmp"

# temporary files are erased when s2p terminates. Switch to False to keep them
cfg['clean_tmp'] = True

# switch to True if you want to process the whole image
cfg['full_img']  = False

# s2p processes the images tile by tile. The tiles are squares cropped from the
# reference image. The width and height of the tiles are given by this param, in pixels.
cfg['tile_size']  = 800

# margins used to increase the footprint of the rectified tiles, to
# account for poor disparity estimation close to the borders
cfg['horizontal_margin'] = 50  # for regularity and occlusions
cfg['vertical_margin'] = 10  # for regularity

# max number of tiles processed in parallel. None means the number of available cores
cfg['max_nb_threads'] = None 

# max number of OMP threads used to process a tile during stereo matching
cfg['omp_num_threads'] = 1

# debug mode: no parallelisation, and more verbose logs
cfg['debug'] = False

# skip tiles for which a file "rectified_disp.tif" already exists
cfg['skip_existing'] = False

# switch to True to translate the planar coordinates of the output point cloud
# around (0, 0)
cfg['offset_ply'] = False

# switch to False to disable cloud colorization
cfg['color_ply'] = True

# resolution of the output digital surface model, in meters per pixel
cfg['dsm_resolution'] = 4

# number of tiles of the output digital surface model
cfg['dsm_nb_tiles'] = 5

# zoom out applied to input images
cfg['subsampling_factor'] = 1

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

# stereo matching algorithm: 'tvl1', 'msmw', 'hirschmuller08',
# hirschmuller08_laplacian', 'sgbm', 'mgm'
cfg['matching_algorithm'] = 'mgm'

# size of the Census NCC square windows used in mgm
cfg['census_ncc_win'] = 5

# blur pleiades images before stereo matching
cfg['use_pleiades_unsharpening'] = True

# this param is used only if the tiles are bigger than 2MPix, to save time on
# sift calls
cfg['pointing_correction_rois_mode'] = 'automatic'

# set these params if you want to impose the disparity range manually
cfg['disp_min'] = None
cfg['disp_max'] = None

# radius for erosion of valid disparity areas. Ignored if less than 2
cfg['msk_erosion'] = 2

cfg['fusion_operator'] = 'average_if_close'

# threshold (in meters) used for the fusion of two dems in triplet processing
# It should be adapted to the zoom factor
cfg['fusion_thresh'] = 3

# set to True to keep only the pixels were the two height maps have a valid
# height. If False, pixels with a valid height in only one of the two maps will
# be kept.
cfg['fusion_conservative'] = False

cfg['disable_srtm'] = False

# method to compute the disparity range: "sift", "srtm" or "wider_sift_srtm"
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

# DSM options
cfg['dsm_option'] = 'median'

# clean height maps outliers
cfg['cargarse_basura'] = True

# longitude/latitude bounding box
cfg['ll_bbx'] = ("-inf", "inf", "-inf", "inf")
