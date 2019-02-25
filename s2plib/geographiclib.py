# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import subprocess


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
