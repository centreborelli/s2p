# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

from tests_utils import TestWithDefaultConfig, data_path
from s2plib import common

class TestCommon(TestWithDefaultConfig):
    """
    Test for the common module
    """
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