# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

from __future__ import print_function

import os
import glob
import numpy as np

import s2p
from s2plib import common
from utils import s2p_mosaic
from tests_utils import data_path


def compare_dsm(computed, expected, absmean_tol, percentile_tol):
    """
    Helper function for end2end tests that compares 2 raster dsm.
    asserts on absolute mean difference and 99% error

    Parameters:
        computed: DSM from test
        expected: Reference DSM
        absmean_tol: Tolerance on absolute mean difference
        percentile_tol: Tolerance on 99% error
    """
    # compare shapes
    np.testing.assert_equal(computed.shape, expected.shape, verbose=True)

    # compare number of valid pixels
    n_computed = np.count_nonzero(np.isfinite(computed))
    n_expected = np.count_nonzero(np.isfinite(expected))

    np.testing.assert_allclose(n_computed, n_expected, rtol=.01, atol=100,
                               verbose=True)
    diff = computed - expected

    # Strip nan from diff
    diff = diff[np.where(np.isfinite(diff))]

    # check mean difference
    meandiff = np.mean(diff)
    print('mean-difference:', meandiff, '(abs. tolerance='+str(absmean_tol)+')')
    assert(np.abs(meandiff) <= absmean_tol)

    # check largest difference
    percentile = np.nanpercentile(np.abs(diff), 99)
    print('99th percentile abs difference', percentile,
          '(tolerance='+str(percentile_tol)+')')
    assert(percentile <= percentile_tol)


def end2end(config_file, ref_dsm, absmean_tol=0.025, percentile_tol=1.):
    print('Configuration file: ', config_file)
    print('Reference DSM:', ref_dsm, os.linesep)

    test_cfg = s2p.read_config_file(config_file)
    s2p.main(test_cfg)

    outdir = test_cfg['out_dir']

    computed = common.gdal_read_as_array_with_nans(os.path.join(outdir, 'dsm.tif'))
    expected = common.gdal_read_as_array_with_nans(ref_dsm)

    compare_dsm(computed, expected, absmean_tol, percentile_tol)


def end2end_mosaic(config_file, ref_height_map, absmean_tol=0.025, percentile_tol=1.):
    test_cfg = s2p.read_config_file(config_file)
    outdir = test_cfg['out_dir']
    test_cfg['skip_existing'] = True
    s2p.main(test_cfg)

    tiles_file = os.path.join(outdir, 'tiles.txt')
    global_height_map = os.path.join(outdir, 'height_map.tif')

    s2p_mosaic.main(tiles_file, global_height_map, 'pair_1/height_map.tif')

    computed = common.gdal_read_as_array_with_nans(global_height_map)
    expected = common.gdal_read_as_array_with_nans(ref_height_map)

    compare_dsm(computed, expected, absmean_tol, percentile_tol)


def test_end2end_pair():
    end2end(data_path('input_pair/config.json'),
            data_path('expected_output/pair/dsm.tif'), 0.025, 1)


def test_end2end_triplet():
    end2end(data_path('input_triplet/config.json'),
            data_path('expected_output/triplet/dsm.tif'), 0.05, 2)


def test_end2end_geo():
    end2end(data_path('input_triplet/config_geo.json'),
            data_path('expected_output/triplet/dsm_geo.tif'), 0.05, 2)


def test_end2end_mosaic():
    end2end_mosaic(data_path('input_triplet/config.json'),
                   data_path('expected_output/triplet/height_map.tif'), 0.05, 2)


def test_distributed_plyflatten():
    config_file = data_path('input_triplet/config.json')

    print('Running end2end with distributed plyflatten dsm ...')
    test_cfg = s2p.read_config_file(config_file)
    test_cfg['skip_existing'] = True
    s2p.main(test_cfg)

    outdir = test_cfg['out_dir']
    computed = common.gdal_read_as_array_with_nans(os.path.join(outdir,
                                                                'dsm.tif'))

    print('Running plyflatten dsm reference ...')

    clouds = '\n'.join(glob.glob(os.path.join(outdir, "tiles", "*", "*",
                                              "cloud.ply")))
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

    expected = common.gdal_read_as_array_with_nans(os.path.join(outdir, 'dsm_ref.tif'))

    compare_dsm(computed, expected, 0, 0)
