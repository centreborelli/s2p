# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os
import numpy as np

from s2plib import rpc_model, rpc_utils
from tests_utils import data_path


def test_matches_from_rpc():
    """
    Test for the rpc_model and rpc_utils modules
    """
    r1 = rpc_model.RPCModel(data_path(os.path.join('input_pair', 'rpc_01.xml')))
    r2 = rpc_model.RPCModel(data_path(os.path.join('input_pair', 'rpc_02.xml')))

    test_matches = rpc_utils.matches_from_rpc(r1, r2, 100, 100, 200, 200, 5)
    expected_matches = np.loadtxt(data_path(os.path.join('expected_output',
                                                         'units',
                                                         'unit_matches_from_rpc.txt')))

    np.testing.assert_equal(test_matches.shape[0], 125, verbose=True)
    np.testing.assert_allclose(test_matches, expected_matches, rtol=0.01,
                               atol=0.1, verbose=True)
