# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>

subsampling_factor = 1
# register using subsampled image
subsampling_factor_registration = 1

# other parameters
sift_match_thresh = 0.4
disp_range_extra_margin = 0.2
# number of ground control points per axis in matches from rpc generation
n_gcp_per_axis = 5
epipolar_thresh = 0.5

# matching algorithm: 'tvl1', 'msmw', 'hirschmuller08', hirschmuller08_laplacian'
matching_algorithm = 'hirschmuller08'

pointing_correction_rois_mode = 'automatic'

use_pleiades_unsharpening = True

# One of "auto_sift", "auto_srtm", "wider_sift_srtm"
disp_range_method = "auto_sift"
disp_range_srtm_low_margin = -10
disp_range_srtm_high_margin = +100

# threshold (in meters) used for the fusion of two dems in triplet processing
# It should be adapted to the zoom factor
fusion_thresh = 3

temporary_dir = "/tmp"
tile_size  = 500
max_nb_threads = -1
clean_tmp = True 

retry = False
mosaic_method = 'piio'

# bunch not used yet
##############################
# handy collector
# http://code.activestate.com/recipes/52308/
class Bunch:
   def __init__(self, **kwds):
      self.__dict__.update(kwds)

