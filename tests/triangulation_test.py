import os
import shutil

import pytest

from s2p import disparity_to_ply, read_config_file
from s2p.config import get_default_config
from s2p.initialization import build_cfg
from s2p.ply import read_3d_point_cloud_from_ply
from s2p.tile import Tile
from tests_utils import data_path


@pytest.mark.parametrize("out_crs", [None, "epsg:32740", "epsg:32740+5773"])
def test_disparity_to_ply(tmp_path, out_crs):
    """
    Check that disparity_to_ply() functions correctly when given
    different out_crs parameters
    """
    cfg = get_default_config()

    # Setup test data
    tile_dir = str(tmp_path / "tile_dir")
    shutil.copytree(data_path("input_triangulation"), tile_dir)

    # Initialize s2p's state
    config_file = data_path(os.path.join("input_pair", "config.json"))
    test_cfg = read_config_file(config_file)
    test_cfg["out_crs"] = out_crs
    build_cfg(cfg, test_cfg)

    tile = Tile(coordinates=(500, 150, 350, 350), dir=tile_dir, neighborhood_dirs=[], json="")
    disparity_to_ply(cfg, tile)

    _, comments = read_3d_point_cloud_from_ply(os.path.join(tile_dir, "cloud.ply"))
    expected_crs = out_crs or "epsg:32740"
    assert comments[-1] == "projection: CRS {}".format(expected_crs)
