# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import subprocess

import pyproj


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


def compute_utm_zone(lon, lat):
    """
    Compute the UTM zone containing the point with given longitude and latitude.

    Args:
        lon (float): longitude of the point
        lat (float): latitude of the point

    Returns:
        str: UTM zone number + hemisphere (eg: '30N')
    """
    # UTM zone number starts from 1 at longitude -180,
    # and increments by 1 every 6 degrees of longitude
    zone = int((lon + 180) // 6 + 1)

    hemisphere = "N" if lat >= 0 else "S"
    utm_zone = "{}{}".format(zone, hemisphere)
    return utm_zone


def epsg_code_from_utm_zone(utm_zone):
    """
    Compute the EPSG code of a given UTM zone.

    Args:
        utm_zone (str): UTM zone number + hemisphere (e.g. "30N" or "30S")

    Returns:
        epsg (int): EPSG code
    """
    zone_number = int(utm_zone[:-1])
    hemisphere = utm_zone[-1]

    if hemisphere not in ["N", "S"]:
        raise ValueError("unknown hemisphere {} in utm_zone {}".format(hemisphere,
                                                                       utm_zone))

    # EPSG = CONST + ZONE where CONST is
    # - 32600 for positive latitudes
    # - 32700 for negative latitudes
    const = 32600 if hemisphere == "N" else 32700
    return const + zone_number


def utm_proj(utm_zone):
    """
    Return a pyproj.Proj object that corresponds to the given utm_zone string.

    Args:
        utm_zone (str): UTM zone number + hemisphere (e.g. "30N" or "30S")

    Returns:
        pyproj.Proj: object that can be used to transform coordinates
    """
    zone_number = utm_zone[:-1]
    hemisphere = utm_zone[-1]
    return pyproj.Proj(
        proj='utm',
        zone=zone_number,
        ellps='WGS84',
        datum='WGS84',
        south=(hemisphere == 'S'),
    )


def pyproj_transform(x, y, in_epsg, out_epsg, z=None):
    """
    Wrapper around pyproj to convert coordinates from an EPSG system to another.

    Args:
        x (scalar or array): x coordinate(s), expressed in in_epsg
        y (scalar or array): y coordinate(s), expressed in in_epsg
        in_epsg (int): EPSG code of the input coordinate system
        out_epsg (int): EPSG code of the output coordinate system
        z (scalar or array): z coordinate(s), expressed in in_epsg

    Returns:
        scalar or array: x coordinate(s), expressed in out_epsg
        scalar or array: y coordinate(s), expressed in out_epsg
        scalar or array (optional if z): z coordinate(s), expressed in out_epsg
    """
    transformer = pyproj.Transformer.from_crs(in_epsg, out_epsg, always_xy=True)
    if z is None:
        return transformer.transform(x, y)
    else:
        return transformer.transform(x, y, z)


def lonlat_to_utm(lon, lat, utm_zone):
    """
    Compute UTM easting and northing of a given lon, lat point.

    Args:
        lon (float): longitude
        lat (float): latitude
        utm_zone (str): UTM zone, e.g. "14N" or "14S"

    Returns:
        easting, northing
    """
    e, n = pyproj_transform(lon, lat, 4326, epsg_code_from_utm_zone(utm_zone))
    return e, n


def lonlat_to_geocentric(lon, lat, alt):
    """
    Compute geocentric cartesian coordinates of a given lon, lat, alt point.

    Args:
        lon (float or list): longitude(s)
        lat (float or list): latitude(s)
        alt (float or list): latitude(s)

    Returns:
        three floats or three lists of floats
    """
    x, y, z = pyproj_transform(lon, lat, 4326, 4978, alt)
    return x, y, z
