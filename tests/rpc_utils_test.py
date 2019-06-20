# s2p (Satellite Stereo Pipeline) testing module

import os
import numpy as np

from s2p import rpc_model, rpc_utils
from tests_utils import data_path


def test_matches_from_rpc():
    """
    Test for rpc_utils.matches_from_rpc().
    """
    r1 = rpc_utils.rpc_from_geotiff(data_path(os.path.join('input_pair', 'img_01.tif')))
    r2 = rpc_utils.rpc_from_geotiff(data_path(os.path.join('input_pair', 'img_02.tif')))

    test_matches = rpc_utils.matches_from_rpc(r1, r2, 100, 100, 200, 200, 5)
    expected_matches = np.loadtxt(data_path(os.path.join('expected_output',
                                                         'units',
                                                         'unit_matches_from_rpc.txt')))

    np.testing.assert_equal(test_matches.shape[0], 125, verbose=True)
    np.testing.assert_allclose(test_matches, expected_matches, rtol=0.01,
                               atol=0.1, verbose=True)
