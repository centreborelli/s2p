# s2p.sift testing module
# Copyright (C) 2019, Carlo de Franchis (CMLA) <carlo.de-franchis@ens-paris-saclay.fr>
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import numpy as np
import rpcm

from s2p import sift
from s2p.config import get_default_config
from tests_utils import data_path


def test_image_keypoints():
    """
    Unit test for the function s2p.sift.image_keypoints.

    Right now it tests only the x, y, scale, orientation of keypoints, not the
    descriptors.
    """
    computed = sift.image_keypoints(data_path('input_triplet/img_02.tif'),
                                    100, 100, 200, 200)
    expected = np.loadtxt(data_path('expected_output/units/unit_image_keypoints.txt'))
    np.testing.assert_allclose(computed[:, :4], expected[:, :4], atol=1e-3)


def test_matching():
    """
    Unit test for the function s2p.sift.keypoints_match.
    """
    computed = sift.keypoints_match(np.loadtxt(data_path('units/sift1.txt')),
                                    np.loadtxt(data_path('units/sift2.txt')))
    expected = np.loadtxt(data_path('expected_output/units/unit_keypoints_match.txt'))
    np.testing.assert_allclose(computed, expected, rtol=0.01, atol=0.1,
                               verbose=True)


def test_matches_on_rpc_roi():
    """
    Unit test for the function sift.matches_on_rpc_roi.
    """
    img1 = data_path('input_triplet/img_01.tif')
    img2 = data_path('input_triplet/img_02.tif')
    rpc1 = rpcm.rpc_from_geotiff(img1)
    rpc2 = rpcm.rpc_from_geotiff(img2)
    cfg = get_default_config()
    computed = sift.matches_on_rpc_roi(
        cfg, img1, img2, rpc1, rpc2, 100, 100, 200, 200,
        method='relative', sift_thresh=0.6, epipolar_threshold=10
    )
    expected = np.loadtxt(data_path('expected_output/units/matches_on_rpc_roi.txt'))
    np.testing.assert_allclose(computed, expected, rtol=0.01, atol=0.1,
                               verbose=True)
