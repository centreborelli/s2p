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

use_pleiades_unsharpening = True



# bunch not used yet
##############################
# handy collector
# http://code.activestate.com/recipes/52308/
class Bunch:
   def __init__(self, **kwds):
      self.__dict__.update(kwds)

