# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2017, Carlo de Franchis <carlo.de-franchis@polytechnique.org>
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os
import sys
import glob
import numpy as np
from plyflatten import plyflatten_from_plyfiles_list

import s2p
from s2p import common
from tests_utils import data_path

here = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(os.path.dirname(here)))
from utils import s2p_mosaic


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
    test_cfg = s2p.read_config_file(config_file)
    s2p.main(test_cfg)

    outdir = test_cfg['out_dir']

    computed = common.rio_read_as_array_with_nans(os.path.join(outdir, 'dsm.tif'))
    expected = common.rio_read_as_array_with_nans(ref_dsm)

    compare_dsm(computed, expected, absmean_tol, percentile_tol)


def test_end2end_pair():
    end2end(data_path('input_pair/config.json'),
            data_path('expected_output/pair/dsm.tif'), 0.025, 1)


def test_end2end_triplet():
    end2end(data_path('input_triplet/config.json'),
            data_path('expected_output/triplet/dsm.tif'), 0.05, 2)


def test_mosaic(tmp_path):

    test_cfg = s2p.read_config_file(data_path('input_triplet/config.json'))

    # we assume that test_end2end_triplet has already run
    s2p_mosaic.main(os.path.join(test_cfg['out_dir'], 'tiles.txt'),
                    os.path.join(tmp_path, 'height_map.tif'),
                    'pair_1/height_map.tif')

    computed = common.rio_read_as_array_with_nans(os.path.join(tmp_path, 'height_map.tif'))
    expected = common.rio_read_as_array_with_nans(data_path('expected_output/triplet/height_map.tif'))

    compare_dsm(computed, expected, absmean_tol=0.05, percentile_tol=2)


def test_distributed_plyflatten():

    test_cfg = s2p.read_config_file(data_path('input_triplet/config.json'))

    # we assume that test_end2end_triplet has already run
    computed = common.rio_read_as_array_with_nans(os.path.join(test_cfg['out_dir'],
                                                               'dsm.tif'))

    clouds_list = glob.glob(os.path.join(test_cfg["out_dir"], "tiles", "*",
                                         "*", "cloud.ply"))
    raster, _ = plyflatten_from_plyfiles_list(clouds_list,
                                              resolution=test_cfg['dsm_resolution'],
                                              roi=None)
    expected = raster[:, :, 0]

    compare_dsm(computed, expected, 0, 0)
