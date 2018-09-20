#!/usr/bin/env python

# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>

from __future__ import print_function

import numpy as np
import argparse
import os
import json
import shutil
import multiprocessing
import collections
import subprocess
import glob

import s2p
from utils import s2p_mosaic
import s2plib


############### Tests functions  #######################

def unit_gdal_version():
    try:
        import gdal
        version_num = int(gdal.VersionInfo('VERSION_NUM'))
        if (version_num < 2010000):
            raise AssertionError(("The version of GDAL should be at least 2.1.\n",
                                  "Recommended fix for Ubuntu 16.04:\n",
                                  "add-apt-repository -y ppa:ubuntugis/ppa\n",
                                  "apt-get update\n",
                                  "apt-get install gdal-bin libgdal-dev\n"))
    except ImportError:
        raise AssertionError('GDAL does not seem to be installed.')


def unit_image_keypoints():

    kpts = s2plib.sift.image_keypoints('testdata/input_triplet/img_02.tif',100,100,200,200)

    test_kpts = np.loadtxt(kpts)
    ref_kpts  = np.loadtxt('testdata/expected_output/units/unit_image_keypoints.txt')

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


def unit_matching():

    test_matches = s2plib.sift.keypoints_match('testdata/units/sift1.txt','testdata/units/sift2.txt')
    expected_matches = np.loadtxt('testdata/expected_output/units/unit_keypoints_match.txt')

    # Check that numbers of matches are the same
    np.testing.assert_equal(test_matches.shape[0],expected_matches.shape[0],verbose=True)

    # Check that all matches are the same
    np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)


# test the plyflatten executable
def unit_plyflatten():
    f = "testdata/input_ply/cloud.ply"                       # input cloud
    e = "testdata/expected_output/plyflatten/dsm_40cm.tiff"  # expected output
    o = s2plib.common.tmpfile(".tiff")                       # actual output
    s2plib.common.run("echo %s | plyflatten 0.4 %s" % (f,o)) # compute dsm
    s = "\"%w %h %v %Y\n\"" # statistics to compare: width,height,avg,numnans
    X = s2plib.common.tmpfile(".txt")
    Y = s2plib.common.tmpfile(".txt")
    s2plib.common.run("imprintf %s %s > %s" % (s, o, X))     # actual stats
    s2plib.common.run("imprintf %s %s > %s" % (s, e, Y))     # expected stats
    s2plib.common.run("diff %s %s" % (X, Y)) # compare stats




def unit_matches_from_rpc():

    rpc1 = s2plib.rpc_model.RPCModel('testdata/input_pair/rpc_01.xml')
    rpc2 = s2plib.rpc_model.RPCModel('testdata/input_pair/rpc_02.xml')

    test_matches = s2plib.rpc_utils.matches_from_rpc(rpc1,rpc2,100,100,200,200,5)
    expected_matches = np.loadtxt('testdata/expected_output/units/unit_matches_from_rpc.txt')

    np.testing.assert_equal(test_matches.shape[0],125,verbose=True)
    np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)


def unit_distributed_plyflatten(config):
    print('Configuration file: ',config)

    print('Running end2end with distributed plyflatten dsm ...')

    test_cfg = s2p.read_config_file(config)
    test_cfg['skip_existing'] = True
    s2p.main(test_cfg)

    outdir = test_cfg['out_dir']
    computed = s2plib.common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))

    print('Running plyflatten dsm reference ...')

    clouds = '\n'.join(glob.glob(os.path.join(outdir, "tiles", "*", "*", "cloud.ply")))
    out_dsm = os.path.join(outdir, "dsm_ref.tif")
    cmd = ['plyflatten', str(test_cfg['dsm_resolution']), out_dsm]
    if 'utm_bbx' in test_cfg:
        bbx = test_cfg['utm_bbx']
        global_xoff = bbx[0]
        global_yoff = bbx[3]
        global_xsize = int(np.ceil((bbx[1]-bbx[0]) / test_cfg['dsm_resolution']))
        global_ysize = int(np.ceil((bbx[3]-bbx[2]) / test_cfg['dsm_resolution']))
        cmd += ['-srcwin', '"{} {} {} {}"'.format(global_xoff, global_yoff,
                                                  global_xsize, global_ysize)]

    run_cmd = "ls %s | %s" % (clouds.replace('\n', ' '), " ".join(cmd))
    s2plib.common.run(run_cmd)

    expected = s2plib.common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm_ref.tif'))

    end2end_compare_dsm(computed,expected,0,0)


def end2end_compare_dsm(computed,expected,absmean_tol,percentile_tol):
    # compare shapes
    np.testing.assert_equal(computed.shape, expected.shape,verbose=True)
    # compare number of valid pixels
    n_computed = np.count_nonzero(np.isfinite(computed))
    n_expected = np.count_nonzero(np.isfinite(expected))

    np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100,verbose=True)

    diff = computed-expected

    # Strip nan from diff
    diff = diff[np.where(np.isfinite(diff))]

    # check mean difference
    meandiff = np.mean(diff)
    print('mean-difference:',meandiff,'(abs. tolerance='+str(absmean_tol)+')')
    assert(np.abs(meandiff) <= absmean_tol)

    # check largest difference
    percentile = np.nanpercentile(np.abs(diff), 99)
    print('99th percentile abs difference',percentile,'(tolerance='+str(percentile_tol)+')')
    assert(percentile <= percentile_tol)


def end2end(config,ref_dsm,absmean_tol=0.025,percentile_tol=1.):

    print('Configuration file: ',config)
    print('Reference DSM:',ref_dsm,os.linesep)

    test_cfg = s2p.read_config_file(config)
    s2p.main(test_cfg)

    outdir = test_cfg['out_dir']

    computed = s2plib.common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))
    expected = s2plib.common.gdal_read_as_array_with_nans(ref_dsm)

    end2end_compare_dsm(computed,expected,absmean_tol,percentile_tol)

def end2end_cluster(config):
    print('Configuration file: ',config)

    print('Running end2end in sequential mode to get reference DSM ...')

    test_cfg = s2p.read_config_file(config)
    test_cfg['skip_existing'] = True
    s2p.main(test_cfg)

    outdir = test_cfg['out_dir']
    expected = s2plib.common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))
    print('Running end2end in cluster mode ...')
    test_cfg_cluster = dict()
    test_cfg_cluster.update(test_cfg)
    test_cfg_cluster['out_dir'] = test_cfg_cluster['out_dir'] + "_cluster"
    test_cfg_cluster['skip_existing'] = True

    print("Running initialisation step ...")
    s2p.main(test_cfg_cluster,["initialisation"])

    # Retrieve tiles list
    outdir = test_cfg_cluster['out_dir']
    tiles_file = os.path.join(outdir,'tiles.txt')

    tiles = s2p.read_tiles(tiles_file)

    print('Found '+str(len(tiles))+' tiles to process')

    for step in s2p.ALL_STEPS:
        if s2p.ALL_STEPS[step] is True:
            print('Running %s on each tile...' % step)
            for tile in tiles:
                print('tile : %s' % tile)
                tile_cfg_cluster = s2p.read_config_file(tile)
                s2p.main(tile_cfg_cluster, [step])
        else:
            print('Running %s...' % step)
            print('test_cfg_cluster : %s' % test_cfg_cluster)
            s2p.main(test_cfg_cluster, [step])

    computed = s2plib.common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))

    end2end_compare_dsm(computed,expected,0,0)

def end2end_mosaic(config,ref_height_map,absmean_tol=0.025,percentile_tol=1.):

    test_cfg = s2p.read_config_file(config)
    outdir = test_cfg['out_dir']
    test_cfg['skip_existing'] = True
    s2p.main(test_cfg)

    tiles_file = os.path.join(outdir,'tiles.txt')
    global_height_map = os.path.join(outdir,'height_map.tif')

    s2p_mosaic.main(tiles_file,global_height_map,'pair_1/height_map.tif')

    computed = s2plib.common.gdal_read_as_array_with_nans(global_height_map)
    expected = s2plib.common.gdal_read_as_array_with_nans(ref_height_map)

    end2end_compare_dsm(computed,expected,absmean_tol,percentile_tol)


############### Registered tests #######################

registered_tests = [('unit_gdal_version', (unit_gdal_version,[])),
                    ('unit_image_keypoints', (unit_image_keypoints,[])),
                    ('unit_matching', (unit_matching,[])),
                    ('unit_plyflatten', (unit_plyflatten,[])),
                    ('unit_matches_from_rpc', (unit_matches_from_rpc,[])),
                    ('end2end_pair', (end2end, ['testdata/input_pair/config.json','testdata/expected_output/pair/dsm.tif',0.025,1])),
                    ('end2end_triplet', (end2end, ['testdata/input_triplet/config.json','testdata/expected_output/triplet/dsm.tif',0.05,2])),
                    ('end2end_cluster', (end2end_cluster, ['testdata/input_triplet/config.json'])),
                    ('end2end_mosaic', (end2end_mosaic, ['testdata/input_triplet/config.json','testdata/expected_output/triplet/height_map.tif',0.05,2])),
                    ('end2end_geometric', (end2end, ['testdata/input_triplet/config_geo.json', 'testdata/expected_output/triplet/dsm_geo.tif',0.05,2])),
                    ('unit_distributed_plyflatten', (unit_distributed_plyflatten, ['testdata/input_triplet/config.json']))]

registered_tests = collections.OrderedDict(registered_tests)


############### Tests main  #######################

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
        for test in registered_tests:
            print('\t'+test)
        exit(1)

    tests_to_run = args.tests

    if args.all:
        tests_to_run = registered_tests.keys()

    print('The following tests will be run: '+str(tests_to_run))

    # First, export the default config to start each test from a clean config
    s2plib.config.cfg["temporary_dir"] = "/tmp"
    test_default_cfg = s2plib.config.cfg.copy()

    # Ensure default temporary dir exists
    if not os.path.isdir(test_default_cfg['temporary_dir']):
        os.mkdir(test_default_cfg['temporary_dir'])

    failed = []

    for test in tests_to_run:
        if test in registered_tests:
            print('Running test '+test+'...'+os.linesep)
            command,args = registered_tests[test]
            try:
                # Ensure each test starts from the default cfg
                s2plib.config.cfg.clear()
                s2plib.config.cfg.update(test_default_cfg)
                command(*args)
                print('Success.'+os.linesep)
            except AssertionError as e:
                print(e)
                print('Failure.'+os.linesep)
                failed.append(test)
        else:
            print('Test '+test+' not found')

    if len(failed)==0:
        print('All tests passed')
        exit(0)
    else:
        print('The following tests failed:')
        for test in failed:
            print('\t'+test)
        exit(1)
