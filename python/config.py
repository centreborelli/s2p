# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>

# This module contains a dictionary, cfg, containing all the parameters of the
# s2p pipeline. This dictionary is updated at runtime with parameters defined
# by the user in the config.json file. All the optional parameters (that the
# user is not forced to define in config.json) must be defined here, otherwise
# they won't have a default value.

cfg = {}
cfg['subsampling_factor'] = 1
# register using subsampled image
cfg['subsampling_factor_registration'] = 1

# other parameters
cfg['sift_match_thresh'] = 0.6
cfg['disp_range_extra_margin'] = 0.2
# number of ground control points per axis in matches from rpc generation
cfg['n_gcp_per_axis'] = 5
cfg['epipolar_thresh'] = 0.5

# matching algorithm: 'tvl1', 'msmw', 'hirschmuller08', hirschmuller08_laplacian', 'sgbm'
cfg['matching_algorithm'] = 'hirschmuller08'

cfg['pointing_correction_rois_mode'] = 'automatic'

cfg['use_pleiades_unsharpening'] = True

# One of "sift", "srtm", "wider_sift_srtm"
cfg['disp_range_method'] = "sift"
cfg['disp_range_srtm_low_margin'] = -10
cfg['disp_range_srtm_high_margin'] = +100
cfg['disp_min'] = None
cfg['disp_max'] = None

# threshold (in meters) used for the fusion of two dems in triplet processing
# It should be adapted to the zoom factor
cfg['fusion_thresh'] = 3

cfg['temporary_dir'] = "/tmp"
cfg['full_img']  = False
cfg['tile_size']  = 1000
cfg['max_nb_threads'] = -1
cfg['clean_tmp'] = True
cfg['debug'] = False
cfg['skip_existing'] = False
cfg['mosaic_method'] = 'piio'
cfg['offset_ply'] = False
cfg['msk_erosion'] = 50
