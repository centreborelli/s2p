# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os

import s2p
from tests_utils import data_path


def test_plyflatten():
    f = data_path("input_ply/cloud.ply")                       # input cloud
    e = data_path("expected_output/plyflatten/dsm_40cm.tiff")  # expected output

    s2p.common.mkdir_p(s2p.config.cfg['temporary_dir'])
    o = s2p.common.tmpfile(".tiff")                       # actual output
    s2p.common.run("echo %s | plyflatten 0.4 %s" % (f, o)) # compute dsm
    s = "\"%w %h %v %Y\n\"" # statistics to compare: width, height, avg, numnans
    X = s2p.common.tmpfile(".txt")
    Y = s2p.common.tmpfile(".txt")
    s2p.common.run("imprintf %s %s > %s" % (s, o, X))     # actual stats
    s2p.common.run("imprintf %s %s > %s" % (s, e, Y))     # expected stats
    s2p.common.run("diff %s %s" % (X, Y)) # compare stats
