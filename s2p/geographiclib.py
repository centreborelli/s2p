# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

from typing import Tuple, Union, List, Optional

import geojson
import numpy as np
import numpy.typing as npt
import pyproj


def geoid_to_ellipsoid(lat: float, lon: float, z: float) -> float:
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


def compute_utm_zone(lon: float, lat: float) -> str:
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


def epsg_code_from_utm_zone(utm_zone: str) -> int:
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


def pyproj_crs(projparams: Union[int, str, dict]) -> pyproj.CRS:
    """
    Wrapper around pyproj to return a pyproj.CRS object that corresponds
    to the given parameters

    Args:
        projparams (int, str, dict): CRS parameters

    Returns:
        pyproj.CRS: object that defines a CRS
    """
    if isinstance(projparams, str):
        try:
            projparams = int(projparams)
        except ValueError:
            pass
    return pyproj.CRS(projparams)


def pyproj_transform(x, y, in_crs, out_crs, z=None):
    """
    Wrapper around pyproj to convert coordinates from an EPSG system to another.

    Args:
        x (scalar or array): x coordinate(s), expressed in in_crs
        y (scalar or array): y coordinate(s), expressed in in_crs
        in_crs (pyproj.CRS or int): input coordinate reference system or EPSG code
        out_crs (pyproj.CRS or int): output coordinate reference system or EPSG code
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


def lonlat_to_utm(lon: float, lat: float, utm_zone: str) -> Tuple[float, float]:
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


def lonlat_to_geocentric(
        lon: Union[float, List[float]],
        lat: Union[float, List[float]],
        alt: Union[float, List[float]]
    ) -> Union[Tuple[float, float, float], Tuple[List[float], List[float], List[float]]]:
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


def read_lon_lat_poly_from_geojson(poly: Union[str, dict]) -> npt.NDArray[np.float64]:
    """
    Read a (lon, lat) polygon from a geojson file or dict.

    Args:
        poly (str or dict): file path to a geojson file containing a single
            polygon, or content of the file as a dict. The geojson's top-level
            type should be either FeatureCollection, Feature, or Polygon.

    Returns:
        numpy array: polygon vertices (lon, lat) coordinates
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


def crs_bbx(ll_poly: npt.NDArray[np.float64],
            crs: Optional[pyproj.CRS] = None,
            align: Optional[float] = None) -> Tuple[int, int, int, int]:
    """
    Compute the bounding box of a (lon, lat) polygon in a given CRS.

    Args:
        ll_poly (np.ndarray): array of shape (n, 2)
        crs (pyproj.CRS): pyproj CRS object. If not specified, the default
            CRS of the UTM zone for the given geography is used.
        align (float): adjust the bounds, by slightly expanding them, in order to
            fit a regular grid of resolution `align` (in crs units)

    Returns:
       4-tuple (left, bottom, right, top)
    """
    if not crs:
        utm_zone = compute_utm_zone(*ll_poly.mean(axis=0))
        epsg = epsg_code_from_utm_zone(utm_zone)
        crs = pyproj_crs(epsg)

    # convert lon lat polygon to target CRS
    easting, northing = pyproj.transform(pyproj.Proj(init="epsg:4326"),
                                         crs, ll_poly[:, 0], ll_poly[:, 1])

    # CRS bounding box
    left = min(easting)
    bottom = min(northing)
    right = max(easting)
    top = max(northing)

    if align:
        left = align * int(np.floor(left / align))
        bottom = align * int(np.floor(bottom / align))
        right = align * int(np.ceil(right / align))
        top = align * int(np.ceil(top / align))

    return left, bottom, right, top
