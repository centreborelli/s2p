subsampling_factor = 4
# register using subsampled image
subsampling_factor_registration = 4

# other parameters
sift_match_thresh = 0.6
disp_range_extra_margin = 0.2
# number of ground control points per axis in matches from rpc generation
n_gcp_per_axis = 5
epipolar_thresh = 2


# matching parameters #  'tvl1','msmw', 'hirschmuller08', hirschmuller08_laplacian'
matching_algorithm = 'tvl1'



# bunch not used yet
##############################
# handy collector
# http://code.activestate.com/recipes/52308/
class Bunch:
   def __init__(self, **kwds):
      self.__dict__.update(kwds)

