import os
import subprocess

import pytest

import s2p
from s2p.config import cfg
from tests_utils import data_path


def test_compute_disparity_map_timeout(timeout=1):
    """
    Run a long call to compute_disparity_map to check that the timeout kills it.
    """
    img = data_path(os.path.join("input_pair", "img_01.tif"))
    disp = data_path(os.path.join("testoutput", "d.tif"))
    mask = data_path(os.path.join("testoutput", "m.tif"))

    with pytest.raises(subprocess.TimeoutExpired):
        s2p.block_matching.compute_disparity_map(cfg, img, img, disp, mask,
                                                 "mgm_multi", -100, 100,
                                                 timeout)


def test_compute_disparity_map_max_disp_range(max_disp_range=10):
    """
    Run a call to compute_disparity_map with a small max_disp_range
    to check that an error is raised.
    """
    img = data_path(os.path.join("input_pair", "img_01.tif"))
    disp = data_path(os.path.join("testoutput", "d.tif"))
    mask = data_path(os.path.join("testoutput", "m.tif"))

    with pytest.raises(s2p.block_matching.MaxDisparityRangeError):
        s2p.block_matching.compute_disparity_map(cfg, img, img, disp, mask,
                                                 "mgm_multi", -100, 100,
                                                 max_disp_range=max_disp_range)
