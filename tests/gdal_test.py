# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import unittest

class TestGdal(unittest.TestCase):
    """
    Test gdal version
    """
    def test_gdal_version(self):
        try:
            from osgeo import gdal
            version_num = int(gdal.VersionInfo('VERSION_NUM'))
            if (version_num < 2010000):
                raise AssertionError(("The version of GDAL should be at least 2.1.\n",
                                      "Recommended fix for Ubuntu 16.04:\n",
                                      "add-apt-repository -y ppa:ubuntugis/ppa\n",
                                      "apt-get update\n",
                                      "apt-get install gdal-bin libgdal-dev\n"))
        except ImportError:
            raise AssertionError('GDAL does not seem to be installed.')