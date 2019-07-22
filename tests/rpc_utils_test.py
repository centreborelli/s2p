# s2p (Satellite Stereo Pipeline) testing module

import os
import numpy as np
import rpcm

from s2p import rpc_utils
from tests_utils import data_path


def test_matches_from_rpc():
    """
    Test for rpc_utils.matches_from_rpc().
    """
    r1 = rpcm.rpc_from_geotiff(data_path(os.path.join('input_pair', 'img_01.tif')))
    r2 = rpcm.rpc_from_geotiff(data_path(os.path.join('input_pair', 'img_02.tif')))

    test_matches = rpc_utils.matches_from_rpc(r1, r2, 100, 100, 200, 200, 5)
    expected_matches = np.loadtxt(data_path(os.path.join('expected_output',
                                                         'units',
                                                         'unit_matches_from_rpc.txt')))

    np.testing.assert_equal(test_matches.shape[0], 125, verbose=True)
    np.testing.assert_allclose(test_matches, expected_matches, rtol=0.01,
                               atol=0.1, verbose=True)


def test_roi_process():
    """
    Test for rpc_utils.roi_process().
    """
    rpc = rpcm.rpc_from_geotiff(data_path(os.path.join('input_pair',
                                                            'img_01.tif')))
    ll_poly = np.asarray([[55.649517, -21.231542],
                          [55.651502, -21.231542],
                          [55.651502, -21.229672],
                          [55.649517, -21.229672]])
    computed = [rpc_utils.roi_process(rpc, ll_poly)[k] for k in
                ['x', 'y', 'w', 'h']]
    expected = (271.48531909338635,
                1.5901905457030807,
                407.3786143153775,
                413.5301010405019)
    np.testing.assert_allclose(computed, expected, atol=1e-3)
