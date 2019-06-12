# s2p (Satellite Stereo Pipeline) testing module

import math

import numpy as np
import rasterio

import s2p
from tests_utils import data_path


def test_plyflatten():
    # Test data
    f = data_path("input_ply/cloud.ply")
    raster, profile = s2p.rasterization.plyflatten_from_plyfiles_list([f], resolution=0.4)
    test_raster = raster[:, :, 0]  # keep only band with height

    # Expected data
    e = data_path("expected_output/plyflatten/dsm_40cm.tiff")
    with rasterio.open(e) as src:
        expected_raster = src.read(1)
        expected_crs = src.crs
        expected_transform = src.transform
        expected_is_tiled = src.is_tiled
        expected_nodata = src.nodata

    # Check that both rasters are equal pixel-wise within a tolerance
    assert np.allclose(test_raster, expected_raster, equal_nan=True)

    # Check that both images have the same CRS
    test_crs = rasterio.crs.CRS.from_proj4(profile['crs'])
    assert test_crs == expected_crs

    # Check that both images have the same transform
    test_transform = profile['transform']
    assert np.allclose(test_transform, expected_transform)

    test_is_tiled = profile['tiled']
    assert test_is_tiled == expected_is_tiled

    test_nodata = profile.get('nodata')
    if expected_nodata and math.isnan(expected_nodata):
        assert math.isnan(test_nodata)
    else:
        assert test_nodata == expected_nodata
