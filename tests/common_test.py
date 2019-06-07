# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os

import s2p
from tests_utils import data_path


def test_plyflatten():
    f = data_path("input_ply/cloud.ply")                       # input cloud
    e = data_path("expected_output/plyflatten/dsm_40cm.tiff")  # expected output
    o = s2p.common.tmpfile(".tiff")                       # actual output
    raster, profile = s2p.rasterization.plyflatten_from_plyfiles_list([f], resolution=0.4) # compute dsm
    s2p.common.rasterio_write(o, raster[:,:,0], profile=profile) # write dsm
    s = "\"%w %h %v %Y\n\"" # statistics to compare: width, height, avg, numnans
    X = s2p.common.tmpfile(".txt")
    Y = s2p.common.tmpfile(".txt")
    s2p.common.run("imprintf %s %s > %s" % (s, o, X))     # actual stats
    s2p.common.run("imprintf %s %s > %s" % (s, e, Y))     # expected stats
    s2p.common.run("diff %s %s" % (X, Y)) # compare stats
