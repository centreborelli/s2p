# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import numpy as np
from nose.tools import assert_equals
from nose.tools import assert_almost_equals
from s2p import evaluation

def test_finite_distances():
    x = np.array([1, 1, 1])
    l = np.array([0, 1, 0])
    assert_equals(evaluation.distance_point_to_line(x, l), 1)
    x = np.array([-1, -1, -1])
    l = np.array([0, 1, 0])
    assert_equals(evaluation.distance_point_to_line(x, l), 1)
    x = np.array([1, 1, 1])
    l = np.array([3, 2, -1])
    assert_equals(evaluation.distance_point_to_line(x, l), 4/np.sqrt(13))

def test_infinite_distances():
    x = np.array([1, 1, 0])
    l = np.array([0, 1, 0])
    assert_equals(evaluation.distance_point_to_line(x, l), np.finfo(float).max)
    x = np.array([1, 1, -7])
    l = np.array([0, 0, 2.3])
    assert_equals(evaluation.distance_point_to_line(x, l), np.finfo(float).max)
