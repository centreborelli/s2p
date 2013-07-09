import subprocess
import numpy as np

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
    p1 = subprocess.Popen(['echo', str(lat), str(lon), str(alt)],
                                             stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['CartConvert'], stdin=p1.stdout,
                                             stdout=subprocess.PIPE)
    line = p2.stdout.readline()
    x = float(line.split()[0])
    y = float(line.split()[1])
    z = float(line.split()[2])
    return x, y, z

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
    p1 = subprocess.Popen(['echo', str(x), str(y), str(z)],
                                                 stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['CartConvert', '-r'], stdin=p1.stdout,
                                                 stdout=subprocess.PIPE)
    line = p2.stdout.readline()
    lat = float(line.split()[0])
    lon = float(line.split()[1])
    alt = float(line.split()[2])
    return lat, lon, alt


def geodetic_to_geocentric_array(lat, lon, alt):
    """
    converts lists of WGS84 (lat, lon, alt) coordinates to lists of geocentric
    cartesian (x, y, z), using the previous function.
    """
    x = np.zeros(len(lat))
    y = np.zeros(len(lat))
    z = np.zeros(len(lat))
    for i in range(len(lat)):
        x[i], y[i], z[i] = geodetic_to_geocentric(lat[i], lon[i], alt[i])

    return x, y, z


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


def srtm4(lon, lat):
    """
    Gives the SRTM height of a point. It is a wrapper to the srtm4 binary.

    Args:
        lon, lat: longitude and latitude

    Returns:
        the height, in meters above the WGS84 geoid (not ellipsoid)
    """
    p = subprocess.Popen(['srtm4', str(lon), str(lat)], stdout=subprocess.PIPE)
    line = p.stdout.readline()
    alt = float(line.split()[0])
    return alt
