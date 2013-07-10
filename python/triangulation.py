#!/usr/bin/env python

from python import common

def compute_point_cloud(rpc1, H1, rpc2, H2, disp, mask, im1, cloud):
    common.run("disp_to_h %s %s %s %s %s %s %s %s" % (rpc1, H1, rpc2, H2, disp,
                                                        mask, im1, cloud))


