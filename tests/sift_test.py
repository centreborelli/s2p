# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

from tests_utils import TestWithDefaultConfig, data_path
from s2plib import sift
import numpy as np

class TestSifts(TestWithDefaultConfig):
    """
    Test for the sift module
    """
    def test_image_keypoints(self):
        #from s2plib import sift
        kpts = sift.image_keypoints(data_path('testdata/input_triplet/img_02.tif'),100,100,200,200)

        test_kpts = np.loadtxt(kpts)
        ref_kpts  = np.loadtxt(data_path('testdata/expected_output/units/unit_image_keypoints.txt'))

        test_set = set(map(tuple,test_kpts[:,0:2]))
        ref_set = set(map(tuple,ref_kpts[:,0:2]))

        print(str(test_kpts.shape[0]-len(test_set))+" spatially redundant kpts found in test")
        print(str(ref_kpts.shape[0]-len(ref_set))+" spatially redundant kpts found in ref")

        common_set = test_set.intersection(ref_set)

        print(str(len(test_set)-len(common_set))+" kpts found in test but not in ref")
        print(str(len(ref_set)-len(common_set))+" kpts found in ref but not in test")

        dist_tol = 0.01

        nb_test_not_in_ref = 0
        for i in range(test_kpts.shape[0]):
            found = False
            for j in range(ref_kpts.shape[0]):
                dist = np.linalg.norm(test_kpts[i,0:1]-ref_kpts[j,0:1])
                if dist<dist_tol:
                    found = True
            if not found:
                print("KeyPoint not found: "+str((test_kpts[i,0:1])))
                nb_test_not_in_ref+=1

        print(str(nb_test_not_in_ref)+" test kpts have no spatially close match in ref")

        nb_ref_not_in_test = 0
        for i in range(test_kpts.shape[0]):
            found = False
            for j in range(ref_kpts.shape[0]):
                dist = np.linalg.norm(test_kpts[i,0:1]-ref_kpts[j,0:1])
                if dist<dist_tol:
                    found = True
            if not found:
                print("KeyPoint not found: "+str((test_kpts[i,0:1])))
                nb_ref_not_in_test+=1

        print(str(nb_ref_not_in_test)+" ref kpts have no spatially close match in test")

        np.testing.assert_equal(nb_ref_not_in_test,0)
        np.testing.assert_equal(nb_test_not_in_ref,0)


    def test_matching(self):

        test_matches = sift.keypoints_match(data_path('testdata/units/sift1.txt'),data_path('testdata/units/sift2.txt'))
        expected_matches = np.loadtxt(data_path('testdata/expected_output/units/unit_keypoints_match.txt'))

        # Check that numbers of matches are the same
        np.testing.assert_equal(test_matches.shape[0],expected_matches.shape[0],verbose=True)

        # Check that all matches are the same
        np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)