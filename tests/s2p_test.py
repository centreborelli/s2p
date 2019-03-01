#!/usr/bin/env python

# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

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
from s2plib import sift, config, rpc_model, rpc_utils, common
import unittest


def data_path(data_path):
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(here,data_path)

class TestWithDefaultConfig(unittest.TestCase):
    def __init__(self, name):
        super(TestWithDefaultConfig,self).__init__(name)
        self.test_default_cfg = config.cfg.copy()

    def setUp(self):
        config.cfg.clear()
        config.cfg.update(self.test_default_cfg)
        # Ensure default temporary dir exists
        if not os.path.isdir(config.cfg['temporary_dir']):
            os.mkdir(config.cfg['temporary_dir'])

    def tearDown(self):
        config.cfg.clear()
        config.cfg.update(self.test_default_cfg)
        common.garbage_cleanup()
      
class TestGdal(unittest.TestCase):
    def test_gdal_version(self):
        try:
            from osgeo import gdal
            version_num = int(gdal.VersionInfo('VERSION_NUM'))
            if (version_num < 2010000):
                raise AssertionError(("The version of GDAL should be at least 2.1.\n",
                                      "Recommended fix for Ubuntu 16.04:\n",
                                      "add-apt-repository -y ppa:ubuntugis/ppa\n",
                                      "apt-get update\n",
                                      "apt-get install gdal-bin libgdal-dev\n"))
        except ImportError:
            raise AssertionError('GDAL does not seem to be installed.')

class TestSifts(TestWithDefaultConfig):
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


    # test the plyflatten executable
    def test_plyflatten(self):
        f = data_path("testdata/input_ply/cloud.ply")                       # input cloud
        e = data_path("testdata/expected_output/plyflatten/dsm_40cm.tiff")  # expected output
        o = common.tmpfile(".tiff")                       # actual output
        common.run("echo %s | plyflatten 0.4 %s" % (f,o)) # compute dsm
        s = "\"%w %h %v %Y\n\"" # statistics to compare: width,height,avg,numnans
        X = common.tmpfile(".txt")
        Y = common.tmpfile(".txt")
        common.run("imprintf %s %s > %s" % (s, o, X))     # actual stats
        common.run("imprintf %s %s > %s" % (s, e, Y))     # expected stats
        common.run("diff %s %s" % (X, Y)) # compare stats


    def test_matches_from_rpc(self):

        rpc1 = rpc_model.RPCModel(data_path('testdata/input_pair/rpc_01.xml'))
        rpc2 = rpc_model.RPCModel(data_path('testdata/input_pair/rpc_02.xml'))

        test_matches = rpc_utils.matches_from_rpc(rpc1,rpc2,100,100,200,200,5)
        expected_matches = np.loadtxt(data_path('testdata/expected_output/units/unit_matches_from_rpc.txt'))

        np.testing.assert_equal(test_matches.shape[0],125,verbose=True)
        np.testing.assert_allclose(test_matches,expected_matches,rtol=0.01,atol=0.1,verbose=True)

def compare_dsm(computed,expected,absmean_tol,percentile_tol):
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
    

class TestEnd2End(TestWithDefaultConfig):
    def test_end2end_pair(self):
        self.end2end(data_path('testdata/input_pair/config.json'),data_path('testdata/expected_output/pair/dsm.tif'),0.025,1)
    def test_end2end_triplet(self):
        self.end2end(data_path('testdata/input_triplet/config.json'),data_path('testdata/expected_output/triplet/dsm.tif'),0.05,2)
    def test_end2end_geo(self):
        self.end2end(data_path('testdata/input_triplet/config_geo.json'), data_path('testdata/expected_output/triplet/dsm_geo.tif'),0.05,2)
    def test_end2end_cluster(self):
        self.end2end_cluster(data_path('testdata/input_triplet/config.json'))
    def test_end2end_mosaic(self):
        self.end2end_mosaic(data_path('testdata/input_triplet/config.json'),data_path('testdata/expected_output/triplet/height_map.tif'),0.05,2)
    def test_distributed_plyflatten(self):
        self.distributed_plyflatten()

    def distributed_plyflatten(self):
        config_file = data_path('testdata/input_triplet/config.json')

        print('Running end2end with distributed plyflatten dsm ...')
        test_cfg = s2p.read_config_file(config_file)
        test_cfg['skip_existing'] = True
        s2p.main(test_cfg)

        outdir = test_cfg['out_dir']
        computed = common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))

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
        common.run(run_cmd)

        expected = common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm_ref.tif'))

        compare_dsm(computed,expected,0,0)

    def end2end(self,config_file,ref_dsm,absmean_tol=0.025,percentile_tol=1.):

        print('Configuration file: ',config_file)
        print('Reference DSM:',ref_dsm,os.linesep)
     
        test_cfg = s2p.read_config_file(config_file)
        s2p.main(test_cfg)

        outdir = test_cfg['out_dir']

        computed = common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))
        expected = common.gdal_read_as_array_with_nans(ref_dsm)

        compare_dsm(computed,expected,absmean_tol,percentile_tol)

    def end2end_cluster(self,config_file):
        print('Configuration file: ',config_file)

        print('Running end2end in sequential mode to get reference DSM ...')
        test_cfg = s2p.read_config_file(config_file)
        test_cfg['skip_existing'] = True

        s2p.main(test_cfg)
        outdir = test_cfg['out_dir']

        expected = common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))
        
        print('Running end2end in cluster mode ...')
        test_cfg_cluster = dict()
        test_cfg_cluster.update(test_cfg)
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

        computed = common.gdal_read_as_array_with_nans(os.path.join(outdir,'dsm.tif'))

        compare_dsm(computed,expected,0,0)

    def end2end_mosaic(self,config_file,ref_height_map,absmean_tol=0.025,percentile_tol=1.):
        test_cfg = s2p.read_config_file(config_file)
        outdir = test_cfg['out_dir']
        test_cfg['skip_existing'] = True
        s2p.main(test_cfg)

        tiles_file = os.path.join(outdir,'tiles.txt')
        global_height_map = os.path.join(outdir,'height_map.tif')

        s2p_mosaic.main(tiles_file,global_height_map,'pair_1/height_map.tif')

        computed = common.gdal_read_as_array_with_nans(global_height_map)
        expected = common.gdal_read_as_array_with_nans(ref_height_map)

        compare_dsm(computed,expected,absmean_tol,percentile_tol)