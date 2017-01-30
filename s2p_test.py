#!/usr/bin/env python

# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>

from __future__ import print_function
import numpy as np
from osgeo import gdal

import s2p

s2p.main('test.json')
computed = gdal.Open('test/dsm.tif').ReadAsArray()
expected = gdal.Open('testdata/expected_output/dsm.tif').ReadAsArray()

# compare shapes
np.testing.assert_equal(computed.shape, expected.shape)

# compare number of valid pixels
n_computed = np.count_nonzero(np.isfinite(computed))
n_expected = np.count_nonzero(np.isfinite(expected))
np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100)

# check largest difference
print('99th percentile abs difference', 
      np.nanpercentile(np.abs(computed - expected), 99))
assert(np.nanpercentile(np.abs(computed - expected), 99) < 1)
