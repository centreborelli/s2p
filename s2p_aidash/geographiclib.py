# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import geojson
from distutils.version import LooseVersion

import pyproj
import numpy as np
import rasterio
from rasterio.crs import CRS as RioCRS
from pyproj.enums import WktVersion


def geoid_to_ellipsoid(lat, lon, z):
    """
    Converts a height, in meters, from the EGM96 geoid datum
    to the WGS84 ellipsoid datum.

    Args:
        lat: latitude, in degrees between -90 and 90
        lon: longitude, between -180 and 180
        z: height, in meters w.r.t. the EGM96 geoid datum

    Returns:
        the height in meters w.r.t. the WGS84 ellipsoid datum

    The conversion is made by PROJ through its python wrapper pyproj
    """
    # WGS84 with ellipsoid height as vertical axis
    ellipsoid = pyproj.CRS.from_epsg(4979)
    # WGS84 with Gravity-related height (EGM96)
    geoid = pyproj.CRS("EPSG:4326+5773")
    transformer = pyproj.Transformer.from_crs(geoid, ellipsoid)
    height = transformer.transform(lat, lon, z)[-1]
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


def rasterio_crs(projparams):
    """
    Return a rasterio.crs.CRS object that corresponds to the given parameters.
    See: https://pyproj4.github.io/pyproj/stable/crs_compatibility.html#converting-from-pyproj-crs-crs-to-rasterio-crs-crs

    Args:
        projparams (int, str, dict, pyproj.CRS): PROJ parameters

    Returns:
        rasterio.crs.CRS: object that can be used with rasterio
    """
    proj_crs = pyproj_crs(projparams)
    if LooseVersion(rasterio.__gdal_version__) < LooseVersion("3.0.0"):
        rio_crs = RioCRS.from_wkt(proj_crs.to_wkt(WktVersion.WKT1_GDAL))
    else:
        rio_crs = RioCRS.from_wkt(proj_crs.to_wkt())
    return rio_crs


def pyproj_crs(projparams):
    """
    Wrapper around pyproj to return a pyproj.crs.CRS object that corresponds
    to the given parameters

    Args:
        projparams (int, str, dict): CRS parameters

    Returns:
        pyproj.crs.CRS: object that defines a CRS
    """
    if isinstance(projparams, str):
        try:
            projparams = int(projparams)
        except (ValueError, TypeError):
            pass
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


def read_lon_lat_poly_from_geojson(poly):
    """
    Read a (lon, lat) polygon from a geojson file or dict.

    Args:
        poly (str or dict): file path to a geojson file containing a single
            polygon, or content of the file as a dict. The geojson's top-level
            type should be either FeatureCollection, Feature, or Polygon.

    Returns:
        list: polygon vertices (lon, lat) coordinates
    """
    if isinstance(poly, str):
        with open(poly, "r") as f:
            a = geojson.load(f)
    else:
        a = poly

    if a["type"] == "FeatureCollection":
        a = a["features"][0]

    if a["type"] == "Feature":
        a = a["geometry"]

    return np.asarray(a["coordinates"][0])


def crs_bbx(ll_poly, crs=None):
    """
    Compute the UTM bounding box of a given (lon, lat) polygon.

    Args:
        ll_poly ()
        crs (pyproj.crs.CRS): pyproj CRS object. If not specified, the default CRS of the UTM
            zone for the given geography is used.

    Returns:
       4-tuple with easting min/max and northing min/max
    """
    if not crs:
        utm_zone = compute_utm_zone(*ll_poly.mean(axis=0))
        epsg = epsg_code_from_utm_zone(utm_zone)
        crs = pyproj_crs(epsg)

    # convert lon lat polygon to target CRS
    easting, northing = pyproj.transform(pyproj.Proj(init="epsg:4326"),
                                         crs, ll_poly[:, 0], ll_poly[:, 1])

    # return UTM bounding box
    east_min = min(easting)
    east_max = max(easting)
    nort_min = min(northing)
    nort_max = max(northing)
    return east_min, east_max, nort_min, nort_max
