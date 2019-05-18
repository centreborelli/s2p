# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os
import numpy as np

import s2p
from s2p import sift
from s2p import config
from tests_utils import data_path


def test_image_keypoints():
    """
    """
    s2p.common.mkdir_p(s2p.config.cfg['temporary_dir'])
    kpts = sift.image_keypoints(data_path('input_triplet/img_02.tif'), 100,
                                100, 200, 200)

    test_kpts = np.loadtxt(kpts)
    ref_kpts = np.loadtxt(data_path('expected_output/units/unit_image_keypoints.txt'))

    test_set = set(map(tuple, test_kpts[:, 0:2]))
    ref_set = set(map(tuple, ref_kpts[:, 0:2]))

    print(str(test_kpts.shape[0]-len(test_set)) + " spatially redundant kpts found in test")
    print(str(ref_kpts.shape[0]-len(ref_set)) + " spatially redundant kpts found in ref")

    common_set = test_set.intersection(ref_set)

    print(str(len(test_set)-len(common_set)) + " kpts found in test but not in ref")
    print(str(len(ref_set)-len(common_set)) + " kpts found in ref but not in test")

    dist_tol = 0.01

    nb_test_not_in_ref = 0
    for i in range(test_kpts.shape[0]):
        found = False
        for j in range(ref_kpts.shape[0]):
            dist = np.linalg.norm(test_kpts[i, 0:1]-ref_kpts[j, 0:1])
            if dist < dist_tol:
                found = True
        if not found:
            print("KeyPoint not found: "+str((test_kpts[i, 0:1])))
            nb_test_not_in_ref += 1

    print(str(nb_test_not_in_ref)+" test kpts have no spatially close match in ref")

    nb_ref_not_in_test = 0
    for i in range(test_kpts.shape[0]):
        found = False
        for j in range(ref_kpts.shape[0]):
            dist = np.linalg.norm(test_kpts[i, 0:1]-ref_kpts[j, 0:1])
            if dist < dist_tol:
                found = True
        if not found:
            print("KeyPoint not found: "+str((test_kpts[i, 0:1])))
            nb_ref_not_in_test += 1

    print(str(nb_ref_not_in_test)+" ref kpts have no spatially close match in test")

    np.testing.assert_equal(nb_ref_not_in_test, 0)
    np.testing.assert_equal(nb_test_not_in_ref, 0)


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
    s2p.common.mkdir_p(s2p.config.cfg['temporary_dir'])

    computed = sift.keypoints_match(data_path('units/sift1.txt'),
                                    data_path('units/sift2.txt'))
    expected = np.loadtxt(data_path('expected_output/units/unit_keypoints_match.txt'))
    assert_arrays_are_equal(computed, expected)


def test_matches_on_rpc_roi():
    """
    Test SIFT matching of two image ROIs.
    """
    img1 = data_path('input_triplet/img_01.tif')
    img2 = data_path('input_triplet/img_02.tif')
    rpc1 = s2p.rpc_utils.rpc_from_geotiff(data_path('input_triplet/img_01.tif'))
    rpc2 = s2p.rpc_utils.rpc_from_geotiff(data_path('input_triplet/img_02.tif'))
    computed = s2p.sift.matches_on_rpc_roi(img1, img2, rpc1, rpc2, 100, 100,
                                           200, 200)
    expected = np.loadtxt(data_path('expected_output/units/matches_on_rpc_roi.txt'))
    assert_arrays_are_equal(computed, expected)
