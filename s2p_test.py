import tifffile
import numpy as np

import s2p

s2p.main('test.json')
computed = tifffile.imread('test/dsm.tif')
expected = tifffile.imread('testdata/expected_output/dsm.tif')
np.testing.assert_allclose(computed, expected, atol=.5)
