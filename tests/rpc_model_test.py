# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

from tests_utils import TestWithDefaultConfig, data_path

from s2plib import rpc_model, rpc_utils
import numpy as np

class TestRpcModel(TestWithDefaultConfig):
    """
    Test for the rpc_model and rpc_utils modules
    """
    def test_matches_from_rpc(self):

        rpc1 = rpc_model.RPCModel(data_path('testdata/input_pair/rpc_01.xml'))
        rpc2 = rpc_model.RPCModel(data_path('testdata/input_pair/rpc_02.xml'))

        test_matches = rpc_utils.matches_from_rpc(rpc1,rpc2,100,100,200,200,5)
        expected_matches = np.loadtxt(data_path('testdata/expected_output/units/unit_matches_from_rpc.txt'))

        np.testing.assert_equal(test_matches.shape[0],125,verbose=True)
        np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)