import os
import subprocess

import pytest

import s2p
from tests_utils import data_path


def test_compute_disparity_map_timeout(timeout=1):
    """
    Run a long call to compute_disparity_map to check that the timeout kills it.
    """
    img = data_path(os.path.join("input_pair", "img_01.tif"))
    disp = data_path(os.path.join("testoutput", "d.tif"))
    mask = data_path(os.path.join("testoutput", "m.tif"))

    with pytest.raises(subprocess.TimeoutExpired):
        s2p.block_matching.compute_disparity_map(img, img, disp, mask,
                                                 "mgm_multi", -100, 100,
                                                 timeout)
