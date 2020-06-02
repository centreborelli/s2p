# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import os
import subprocess

import pyproj
import numpy as np


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


def crs_proj(projparams):
    """
    Return a pyproj.Proj object that corresponds
    to the given parameters

    Args:
        projparams (int, str, dict, pyproj.CRS): PROJ parameters

    Returns:
        pyproj.Proj: object that can be used to transform coordinates
    """
    return pyproj.Proj(projparams)


def pyproj_crs(projparams):
    """
    Wrapper around pyproj to return a pyproj.crs.CRS object that corresponds
    to the given parameters

    Args:
        projparams (int, str, dict): CRS parameters

    Returns:
        pyproj.crs.CRS: object that defines a CRS
    """
    return pyproj.crs.CRS(projparams)


def pyproj_transform(x, y, in_crs, out_crs, z=None):
    """
    Wrapper around pyproj to convert coordinates from an EPSG system to another.

    Args:
        x (scalar or array): x coordinate(s), expressed in in_crs
        y (scalar or array): y coordinate(s), expressed in in_crs
        in_crs (pyproj.crs.CRS or int): input coordinate reference system or EPSG code
        out_crs (pyproj.crs.CRS or int): output coordinate reference system or EPSG code
        z (scalar or array): z coordinate(s), expressed in in_crs

    Returns:
        scalar or array: x coordinate(s), expressed in out_crs
        scalar or array: y coordinate(s), expressed in out_crs
        scalar or array (optional if z): z coordinate(s), expressed in out_crs
    """
    transformer = pyproj.Transformer.from_crs(in_crs, out_crs, always_xy=True)
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


def read_lon_lat_poly_from_geojson(geojson):
    """
    Read a (lon, lat) polygon from a geojson file or dict.

    Args:
        geojson: file path to a geojson file containing a single polygon,
            or content of the file as a dict.
            The geojson's top-level type should be either FeatureCollection,
            Feature, or Polygon.

    Returns:
        ll_poly: list of polygon vertices (lon, lat) coordinates
    """
    # extract lon lat from geojson file or dict
    if isinstance(geojson, str):
        with open(geojson, 'r') as f:
            a = json.load(f)
    else:
        a = geojson

    if a["type"] == "FeatureCollection":
        a = a["features"][0]

    if a["type"] == "Feature":
        a = a["geometry"]

    ll_poly = np.array(a["coordinates"][0])
    return ll_poly


def utm_bbx(ll_poly, utm_zone=None):
    """
    Compute the UTM bounding box of a given (lon, lat) polygon.

    Args:
        ll_poly ()
        utm_zone (): force the UTM zone number. If not specified, the default UTM
            zone for the given geography is used.

    Returns:
       4-tuple with easting min/max and northing min/max
    """
    if not utm_zone:
        utm_zone = compute_utm_zone(*ll_poly.mean(axis=0))

    # convert lon lat polygon to UTM
    utm_proj = utm_proj(utm_zone)
    easting, northing = pyproj.transform(pyproj.Proj(init="epsg:4326"),
                                         utm_proj, ll_poly[:, 0], ll_poly[:, 1])

    # return UTM bounding box
    east_min = min(easting)
    east_max = max(easting)
    nort_min = min(northing)
    nort_max = max(northing)
    return east_min, east_max, nort_min, nort_max
