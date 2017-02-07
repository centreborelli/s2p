#!/usr/bin/env python

# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>

from __future__ import print_function
import numpy as np
from osgeo import gdal
import argparse
import os

import s2p

def end2end_pair():
    s2p.main('testdata/input_pair/config.json')
    computed = gdal.Open('test_pair/dsm.tif').ReadAsArray()
    expected = gdal.Open('testdata/expected_output/pair/dsm.tif').ReadAsArray()
    
    # compare shapes
    np.testing.assert_equal(computed.shape, expected.shape)
    # compare number of valid pixels
    n_computed = np.count_nonzero(np.isfinite(computed))
    n_expected = np.count_nonzero(np.isfinite(expected))
    np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100)
    
    # check largest difference
    print('99th percentile abs difference', 
          np.nanpercentile(np.abs(computed - expected), 99))
    assert(np.nanpercentile(np.abs(computed - expected), 99) < 1)

def end2end_triplet():
    s2p.main('testdata/input_triplet/config.json')
    computed = gdal.Open('test_triplet/dsm.tif').ReadAsArray()
    expected = gdal.Open('testdata/expected_output/triplet/dsm.tif').ReadAsArray()
    
    # compare shapes
    np.testing.assert_equal(computed.shape, expected.shape)
    # compare number of valid pixels
    n_computed = np.count_nonzero(np.isfinite(computed))
    n_expected = np.count_nonzero(np.isfinite(expected))
    np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100)
    
    # check largest difference
    print('99th percentile abs difference', 
          np.nanpercentile(np.abs(computed - expected), 99))
    assert(np.nanpercentile(np.abs(computed - expected), 99) < 1)

# Register tests
registered_tests = { 'end2end_pair' : (end2end_pair, []),
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
        
    for test in tests_to_run:
        if test in registered_tests:
            print('Running test '+test+'...'+os.linesep)
            command,args = registered_tests[test]
            command(*args)
            print('Success.'+os.linesep)
        else:
            print('Test '+test+' not found')
            exit(1)
    exit(0)
            
