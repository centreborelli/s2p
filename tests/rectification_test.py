import os

import numpy as np
import pytest
from rpcm import rpc_from_geotiff

import s2p
from tests_utils import data_path


@pytest.fixture(name='matches')
def fixture_matches():
    matches = np.loadtxt(
        data_path(os.path.join('expected_output', 'units', 'unit_matches_from_rpc.txt'))
    )
    return matches


@pytest.fixture(name='images')
def fixture_images():
    res = []
    for i in [1, 2]:
        im = data_path(os.path.join('input_pair', 'img_0{}.tif'.format(i)))
        rpc = rpc_from_geotiff(im)
        res.append(im)
        res.append(rpc)
    return res


def test_rectification_homographies(matches):
    """
    Test for rectification.rectification_homographies().
    """
    x, y, w, h = 100, 100, 200, 200
    H1, H2, F = s2p.rectification.rectification_homographies(matches, x, y, w, h)

    for variable, filename in zip([H1, H2, F], ['H1.txt', 'H2.txt', 'F.txt']):
        expected = np.loadtxt(data_path(os.path.join('expected_output', 'units', filename)))
        np.testing.assert_allclose(variable, expected, rtol=0.01, atol=1e-6, verbose=True)


def test_rectify_pair_with_matches(tmp_path, matches, images):
    """
    Test running rectification.rectify_pair() with some matches.
    """
    im1, rpc1, im2, rpc2 = images
    s2p.rectification.rectify_pair(
        im1=im1,
        im2=im2,
        rpc1=rpc1,
        rpc2=rpc2,
        x=100,
        y=100,
        w=200,
        h=200,
        out1=str(tmp_path / 'out1.tiff'),
        out2=str(tmp_path / 'out2.tiff'),
        sift_matches=matches,
        method='sift',
    )
