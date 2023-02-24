# s2p (Satellite Stereo Pipeline) testing module

import os
import numpy as np
import rpcm
import pytest

from s2p import rpc_utils
from s2p.config import get_default_config
from tests_utils import data_path


def test_matches_from_rpc():
    """
    Test for rpc_utils.matches_from_rpc().
    """
    r1 = rpcm.rpc_from_geotiff(data_path(os.path.join('input_pair', 'img_01.tif')))
    r2 = rpcm.rpc_from_geotiff(data_path(os.path.join('input_pair', 'img_02.tif')))

    cfg = get_default_config()
    test_matches = rpc_utils.matches_from_rpc(cfg, r1, r2, 100, 100, 200, 200, 5)
    expected_matches = np.loadtxt(data_path(os.path.join('expected_output',
                                                         'units',
                                                         'unit_matches_from_rpc.txt')))

    np.testing.assert_equal(test_matches.shape[0], 125, verbose=True)
    np.testing.assert_allclose(test_matches, expected_matches, rtol=0.01,
                               atol=0.1, verbose=True)


@pytest.mark.parametrize(
    "use_srtm, exogenous_dem, exogenous_dem_geoid_mode, expected",
    [
        (
            False,
            None,
            True,
            (271.48531, 1.59019, 407.37861, 413.53010),
        ),
        (
            True,
            None,
            True,
            (353.49632, 296.69818, 408.16015, 413.54849),
        ),
        (
            False,
            data_path(os.path.join("expected_output", "pair", "dsm.tif")),
            True,
            (356.65154, 308.01931, 408.19018, 413.54920)
        ),
        (
            False,
            data_path(os.path.join("expected_output", "pair", "dsm.tif")),
            False,
            (356.46596, 307.35347, 408.18841, 413.54916),
        ),
    ],
)
def test_roi_process(use_srtm, exogenous_dem, exogenous_dem_geoid_mode, expected):
    """
    Test for rpc_utils.roi_process().
    """
    rpc = rpcm.rpc_from_geotiff(data_path(os.path.join("input_pair", "img_01.tif")))
    ll_poly = np.asarray(
        [
            [55.649517, -21.231542],
            [55.651502, -21.231542],
            [55.651502, -21.229672],
            [55.649517, -21.229672],
        ]
    )
    output = rpc_utils.roi_process(
        rpc,
        ll_poly,
        use_srtm=use_srtm,
        exogenous_dem=exogenous_dem,
        exogenous_dem_geoid_mode=exogenous_dem_geoid_mode,
    )
    computed = [output[k] for k in ["x", "y", "w", "h"]]
    np.testing.assert_allclose(computed, expected, atol=1e-3)
