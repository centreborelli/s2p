import json
import os
import shutil
from unittest.mock import MagicMock

import rasterio
import rpcm

import pytest
import s2p
from s2p.config import get_default_config
from tests_utils import data_path


@pytest.fixture(name="mocks")
def fixture_mocks():
    """
    Mock functions whose call patterns we want to check
    """
    rpcm.rpc_from_geotiff = MagicMock(side_effect=rpcm.rpc_from_geotiff)
    rpcm.rpc_from_rpc_file = MagicMock(side_effect=rpcm.rpc_from_rpc_file)


@pytest.fixture(name="data")
def fixture_data(tmp_path):
    """
    Copy the test data to a temporary directory
    """
    img1 = data_path(os.path.join("input_pair", "img_01.tif"))
    with rasterio.open(img1) as f:
        rpc1 = rpcm.RPCModel(f.tags(ns="RPC"))
    shutil.copy(img1, tmp_path / "img_01.tif")

    img2 = data_path(os.path.join("input_pair", "img_02.tif"))
    with rasterio.open(img2) as f:
        rpc2 = rpcm.RPCModel(f.tags(ns="RPC"))
    shutil.copy(img2, tmp_path / "img_02.tif")

    config = data_path(os.path.join("input_pair", "config.json"))
    tmp_config = tmp_path / "config.json"
    shutil.copy(config, tmp_config)

    return tmp_config, tmp_path, rpc1, rpc2


def test_no_rpc(data, mocks):
    """
    Initialize s2p with no `rpc` key.
    The RPCs should be read from the geotiff tags.
    """

    tmp_config, _, _, _ = data
    user_cfg = s2p.read_config_file(tmp_config)
    cfg = get_default_config()
    s2p.initialization.build_cfg(cfg, user_cfg)

    rpcm.rpc_from_geotiff.assert_called()
    assert rpcm.rpc_from_geotiff.call_count == 2
    rpcm.rpc_from_rpc_file.assert_not_called()


def test_rpc_path(data, mocks):
    """
    Initialize s2p with `rpc` keys that are paths to text files.
    The RPCs should be loaded from the text files.
    """
    tmp_config, tmp_path, rpc1, rpc2 = data
    with open(tmp_config) as f:
        cfg = json.load(f)

    rpc1_path = str(tmp_path / "rpc1.txt")
    rpc1.write_to_file(rpc1_path)
    cfg["images"][0]["rpc"] = rpc1_path

    rpc2_path = str(tmp_path / "rpc2.txt")
    rpc2.write_to_file(rpc2_path)
    cfg["images"][1]["rpc"] = rpc2_path

    with open(tmp_config, "w") as f:
        json.dump(cfg, f)

    cfg = get_default_config()
    user_cfg = s2p.read_config_file(tmp_config)
    s2p.initialization.build_cfg(cfg, user_cfg)

    rpcm.rpc_from_geotiff.assert_not_called()
    rpcm.rpc_from_rpc_file.assert_called()
    assert rpcm.rpc_from_rpc_file.call_count == 2


def test_rpc_dict(data, mocks):
    """
    Initialize s2p with `rpc` keys that are dicts with the RPC contents.
    The RPCs should be loaded from the dicts.
    """
    tmp_config, _, rpc1, rpc2 = data
    with open(tmp_config) as f:
        cfg = json.load(f)

    cfg["images"][0]["rpc"] = rpc1.__dict__
    cfg["images"][1]["rpc"] = rpc2.__dict__

    with open(tmp_config, "w") as f:
        json.dump(cfg, f)

    cfg = get_default_config()
    user_cfg = s2p.read_config_file(tmp_config)
    s2p.initialization.build_cfg(cfg, user_cfg)

    rpcm.rpc_from_geotiff.assert_not_called()
    rpcm.rpc_from_rpc_file.assert_not_called()


def test_roi_geojson(data):
    tmp_config, _, _, _ = data
    cfg = get_default_config()
    user_cfg = s2p.read_config_file(tmp_config)

    user_cfg["roi_geojson"] = {
      "coordinates" : [
        [
          [
            55.64943405,
            -21.23207174
          ],
          [
            55.65212062,
            -21.23207174
          ],
          [
            55.65212062,
            -21.23460474
          ],
          [
            55.64943405,
            -21.23460474
          ],
          [
            55.64943405,
            -21.23207174
          ]
        ]
      ],
      "type" : "Polygon"
    }

    s2p.initialization.build_cfg(cfg, user_cfg)
    assert user_cfg["roi"] == {'x': 150, 'y': 150, 'w': 700, 'h': 700}
