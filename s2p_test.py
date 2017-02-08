#!/usr/bin/env python

# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>

from __future__ import print_function
import numpy as np
from osgeo import gdal
import argparse
import os

import s2p
import s2plib

def unit_image_keypoints():
    try:
        os.mkdir('s2p_tmp')
    except:
        pass
    kpts = s2plib.sift.image_keypoints('testdata/input_triplet/img_02.tif',100,100,200,200)

    test_kpts = np.loadtxt(kpts)
    ref_kpts  = np.loadtxt('testdata/expected_output/units/unit_image_keypoints.txt')

    # Check that the number of keypoints is the same
    np.testing.assert_equal(test_kpts.shape[0],ref_kpts.shape[0],verbose=True)
    
    # Check that all keypoints are the same
    np.testing.assert_allclose(test_kpts, ref_kpts, rtol=.01, atol=1,verbose=True)

def unit_matching():
    try:
        os.mkdir('s2p_tmp')
    except:
        pass
    test_matches = s2plib.sift.keypoints_match('testdata/units/sift1.txt','testdata/units/sift2.txt')
    expected_matches = np.loadtxt('testdata/expected_output/units/unit_keypoints_match.txt')

    # Check that numbers of matches are the same
    np.testing.assert_equal(test_matches.shape[0],expected_matches.shape[0],verbose=True)
    
    # Check that all matches are the same
    np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)
    
def unit_matches_from_rpc():
    try:
        os.mkdir('s2p_tmp')
    except:
        pass

    s2plib.config.cfg['disable_srtm'] = True
    rpc1 = s2plib.rpc_model.RPCModel('testdata/input_pair/rpc_01.xml')
    rpc2 = s2plib.rpc_model.RPCModel('testdata/input_pair/rpc_02.xml')

    test_matches = s2plib.rpc_utils.matches_from_rpc(rpc1,rpc2,100,100,200,200,5)
    expected_matches = np.loadtxt('testdata/expected_output/units/unit_matches_from_rpc.txt')
    s2plib.config.cfg['disable_srtm'] = False
    
    np.testing.assert_equal(test_matches.shape[0],125,verbose=True)
    np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)

def end2end_pair():
    s2p.main('testdata/input_pair/config.json')
    computed = gdal.Open('test_pair/dsm.tif').ReadAsArray()
    expected = gdal.Open('testdata/expected_output/pair/dsm.tif').ReadAsArray()
    
    # compare shapes
    np.testing.assert_equal(computed.shape, expected.shape,verbose=True)
    # compare number of valid pixels
    n_computed = np.count_nonzero(np.isfinite(computed))
    n_expected = np.count_nonzero(np.isfinite(expected))
    np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100,verbose=True)
    
    # check largest difference
    print('99th percentile abs difference', 
          np.nanpercentile(np.abs(computed - expected), 99))
    assert(np.nanpercentile(np.abs(computed - expected), 99) < 1)

def end2end_triplet():
    s2p.main('testdata/input_triplet/config.json')
    computed = gdal.Open('test_triplet/dsm.tif').ReadAsArray()
    expected = gdal.Open('testdata/expected_output/triplet/dsm.tif').ReadAsArray()
    
    # compare shapes
    np.testing.assert_equal(computed.shape, expected.shape,verbose=True)
    # compare number of valid pixels
    n_computed = np.count_nonzero(np.isfinite(computed))
    n_expected = np.count_nonzero(np.isfinite(expected))
    np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100,verbose=True)
    
    # check largest difference
    print('99th percentile abs difference', 
          np.nanpercentile(np.abs(computed - expected), 99))
    assert(np.nanpercentile(np.abs(computed - expected), 99) < 1)

# Register tests
registered_tests = { 'unit_image_keypoints' : (unit_image_keypoints,[]),
                     'unit_matching' : (unit_matching,[]),
                     'unit_matches_from_rpc' : (unit_matches_from_rpc,[]),
                     'end2end_pair' : (end2end_pair, []),
                     'end2end_triplet' : (end2end_triplet, [])}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=('S2P: Satellite Stereo '
                                                  'Pipeline tests suite'))
    
    parser.add_argument('--tests',nargs='+',metavar='test', help='(List of tests to run)')
    parser.add_argument('--all',help=('Run all tests'),action='store_true')

    args = parser.parse_args()

    # If nothing provided, display list of tests
    if not args.tests and not args.all:
        parser.print_help()
        
        print(os.linesep+'available tests:')
        for test,commands in registered_tests.iteritems():
            print('\t'+test)
        exit(1)

    tests_to_run = args.tests

    if args.all:
        tests_to_run = registered_tests.keys()

    print('The following tests will be run: '+str(tests_to_run))

    # First, export the default config to start each test from a clean config
    test_default_cfg = s2plib.config.cfg
    
    failed = []
    
    for test in tests_to_run:
        if test in registered_tests:
            print('Running test '+test+'...'+os.linesep)
            command,args = registered_tests[test]
            try:
                # Ensure each test starts from the default cfg
                s2plib.config.cfg = test_default_cfg
                command(*args)
                print('Success.'+os.linesep)
            except AssertionError as e:
                print(e)
                print('Failure.'+os.linesep)
                failed.append(test)
        else:
            print('Test '+test+' not found')
    
    if len(failed)==0:
        print('All tests passes')
        exit(0)
    else:
        print('The following tests failed:')
        for test in failed:
            print('\t'+test)
        exit(1)
