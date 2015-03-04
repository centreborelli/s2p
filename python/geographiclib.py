# Copyright (C) 2013, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2013, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2013, Julien Michel <julien.michel@cnes.fr>

import subprocess
import numpy as np
import common

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


def geodetic_to_utm(lat, lon, zone=None):
    """
    Converts WGS84 ellipsoidal coordinates to UTM coordinates, using
    the most appropriate zone. Please note that forcing an UTM zone
    may raise an exception if the lat, lon falls too far from the
    requested UTM zone.

    Args:
        lat: latitude, in degrees between -90 and 90
        lon: longitude, between -180 and 180
        zone: None (automatic zone determination) or string ('31N' for instance)

    Returns:
        x, y, zone: the UTM coordinates of the input point, and the zone
    """

    command = ['GeoConvert']

    # if zone is None, we let GeoConvert guess the zone
    if zone is None:
        command.append("-u")
        command.append("-s")
    else:
        command.append("-u")
        command.append("-z")
        command.append(zone)

    p1 = subprocess.Popen(['echo', str(lat), str(lon)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(command, stdin=p1.stdout, stdout=subprocess.PIPE)
    line = p2.stdout.readline()

    splits = line.split()
    zone = splits[0]
    x = float(splits[1])
    y = float(splits[2])

    return x, y, zone


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
    p1 = subprocess.Popen(['echo', str(lat), str(lon)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['GeoidEval'], stdin=p1.stdout,
                                                        stdout=subprocess.PIPE)
    line = p2.stdout.readline()
    h = float(line.split()[0])
    return h
