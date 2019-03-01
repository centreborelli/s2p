# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os
from s2plib import config, common
import unittest

def data_path(data_path):
    """
    Build an absolute data path from an input test datafile
    params:
        data_path Path to the input test data
    returns:
        absolute path to that data file
    """
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(here,data_path)

class TestWithDefaultConfig(unittest.TestCase):
    """
    Base class for all tests classes requiring to reset config between each test
    """
    def __init__(self, name):
        super(TestWithDefaultConfig,self).__init__(name)
        # Init saves default config
        self.test_default_cfg = config.cfg.copy()

    def setUp(self):
        """
        Reset config properly, also ensures that temporary_dir exists
        """
        config.cfg.clear()
        config.cfg.update(self.test_default_cfg)
        # Ensure default temporary dir exists
        if not os.path.isdir(config.cfg['temporary_dir']):
            os.mkdir(config.cfg['temporary_dir'])

    def tearDown(self):
        """
        Reset config and call garbage cleanup
        """
        config.cfg.clear()
        config.cfg.update(self.test_default_cfg)
        common.garbage_cleanup()
      