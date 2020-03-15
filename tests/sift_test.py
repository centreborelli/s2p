# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Carlo de Franchis (CMLA) <carlo.de-franchis@ens-paris-saclay.fr>
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os
import numpy as np
import rpcm

import s2p
from s2p import sift
from s2p import config
from tests_utils import data_path


def test_image_keypoints():
    """
    """
    kpts = sift.image_keypoints(data_path('input_triplet/img_02.tif'), 100,
                                100, 200, 200)

    ref_kpts = np.loadtxt(data_path('expected_output/units/unit_image_keypoints.txt'))
    np.testing.assert_allclose(kpts[:, :2], ref_kpts[:, :2], atol=1e-3)


def assert_arrays_are_equal(a, b, rtol=0.01, atol=0.1, verbose=True):
    """
    Test if two numpy arrays are equal up to tolerance.
    """
    # check that the shapes are the same
    np.testing.assert_equal(a.shape, b.shape, verbose=verbose)

    # check that the arrays elements are the same
    np.testing.assert_allclose(a, b, rtol=rtol, atol=atol, verbose=verbose)


def test_matching():
    """
    """
    computed = sift.keypoints_match(np.loadtxt(data_path('units/sift1.txt')),
                                    np.loadtxt(data_path('units/sift2.txt')))
    expected = np.loadtxt(data_path('expected_output/units/unit_keypoints_match.txt'))
    assert_arrays_are_equal(computed, expected)


def test_matches_on_rpc_roi():
    """
    Test SIFT matching of two image ROIs.
    """
    img1 = data_path('input_triplet/img_01.tif')
    img2 = data_path('input_triplet/img_02.tif')
    rpc1 = rpcm.rpc_from_geotiff(img1)
    rpc2 = rpcm.rpc_from_geotiff(img2)
    computed = s2p.sift.matches_on_rpc_roi(img1, img2, rpc1, rpc2, 100, 100,
                                           200, 200)
    expected = np.loadtxt(data_path('expected_output/units/matches_on_rpc_roi.txt'))
    assert_arrays_are_equal(computed, expected)
