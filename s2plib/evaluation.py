# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from __future__ import print_function
import numpy as np

def distance_point_to_line(x, l):
    """
    Computes the distance between a point and a line expressed in homogeneous
    coordinates.

    Args:
        x: 3-vector, containing the homogeneous coordinates of a point
        l: 3-vector, containing the homogeneous coordinates of a line

    Returns:
        the distance between x and l
        In the case where x is an ideal point or l is the line at infinity, the
        returned distance is infinity
    """
    if (np.abs(x[2]) < np.finfo(float).eps):
        # x is an ideal point, distance is +infty
        return np.finfo(float).max
    if (np.hypot(l[0], l[1]) < np.finfo(float).eps):
        # l is the line at infinity, distance is +infty
        return np.finfo(float).max

    num = np.abs(np.dot(x, l))
    den = np.hypot(l[0], l[1]) * np.abs(x[2])
    return num/den
