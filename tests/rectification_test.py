import os
import numpy as np

import s2p
from tests_utils import data_path


def test_rectification_homographies():
    """
    Test for rectification.rectification_homographies().
    """
    matches = np.loadtxt(data_path(os.path.join('expected_output', 'units',
                                                'unit_matches_from_rpc.txt')))
    
    x, y, w, h = 100, 100, 200, 200
    H1, H2, F = s2p.rectification.rectification_homographies(matches, x, y, w, h)

    for variable, filename in zip([H1, H2, F], ['H1.txt', 'H2.txt', 'F.txt']):
        expected = np.loadtxt(data_path(os.path.join('expected_output', 'units', filename)))
        np.testing.assert_allclose(variable, expected, rtol=0.01, atol=1e-6, verbose=True)
