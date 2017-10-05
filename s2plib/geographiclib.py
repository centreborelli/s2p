# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import subprocess
import numpy as np
from s2plib import common

def geodetic_to_geocentric(lat, lon, alt):
    """
    Converts WGS84 ellipsoidal coordinates to geocentric cartesian coordinates.

    Args:
        lat: latitude, in degrees between -90 and 90
        lon: longitude, between -180 and 180
        alt: altitude, in meters, above the WGS84 reference ellipsoid

    Returns:
        x, y, z: the cartesian coordinates of the input point, expressed in the
            geocentric frame.

    The conversion is made by a commandline tool, CartConvert, which is part of
    the GeographicLib library:
    http://geographiclib.sourceforge.net/html/intro.html
    """
    pts = np.vstack([lat, lon, alt]).T
    out = common.run_binary_on_list_of_points(pts, 'CartConvert')
    return out[:, 0], out[:, 1], out[:, 2]


def geocentric_to_geodetic(x, y, z):
    """
    Converts geocentric cartesian coordinates to WGS84 ellipsoidal coordinates.

    Args:
        x, y, z: the cartesian coordinates of the input point, expressed in the
            geocentric frame.

    Returns:
        lat: latitude, in degrees between -90 and 90
        lon: longitude, between -180 and 180
        alt: altitude, in meters, above the WGS84 reference ellipsoid

    The conversion is made by a commandline tool, CartConvert, which is part of
    the GeographicLib library:
    http://geographiclib.sourceforge.net/html/intro.html
    """
    pts = np.vstack([x, y, z]).T
    out = common.run_binary_on_list_of_points(pts, 'CartConvert', '-r')
    return out[:, 0], out[:, 1], out[:, 2]


def geodetic_to_mercator(lat, lon, ref_lon=0):
    """
    Converts WGS84 ellipsoidal coordinates to mercator coordinates, using a
    reference longitude.

    Args:
        lat: latitude, in degrees between -90 and 90
        lon: longitude, between -180 and 180
        ref_lon (optional, default 0): reference longitude, in degrees

    Returns:
        x, y: the mercator coordinates of the input point
    """
    r = 6378.1 * 1000
    c = 2 * np.pi / 360.0
    x = r * (lon - ref_lon) * c
    y = r * np.log( (1 + np.sin(lat*c)) / np.cos(lat*c))
    return x, y


def geoid_above_ellipsoid(lat, lon):
    """
    Computes the height, in meters, of the EGM96 geoid above the WGS84 ellipsoid.

    Args:
        lat: latitude, in degrees between -90 and 90
        lon: longitude, between -180 and 180

    Returns:
        the height in meters. It should be between -106 and 85 meters.


    The conversion is made by a commandline tool, GeoidEval, which is part of
    the GeographicLib library:
    http://geographiclib.sourceforge.net/html/intro.html
    """
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    geoid_dir = os.path.join(parent_dir, 'c')
    geoid_name = 'egm96-15'
    p = subprocess.Popen(['echo', str(lat), str(lon)], stdout=subprocess.PIPE)
    q = subprocess.Popen(['GeoidEval',
                          '-d', geoid_dir,
                          '-n', geoid_name], stdin=p.stdout, stdout=subprocess.PIPE)
    height = float(q.stdout.readline().split()[0])
    return height
