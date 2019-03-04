"""
RPC model parsers, localization, and projection
Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
"""

import numpy as np
import utm


def apply_poly(poly, x, y, z):
    """
    Evaluates a 3-variables polynom of degree 3 on a triplet of numbers.

    Args:
        poly: list of the 20 coefficients of the 3-variate degree 3 polynom,
            ordered following the RPC convention.
        x, y, z: triplet of floats. They may be numpy arrays of same length.

    Returns:
        the value(s) of the polynom on the input point(s).
    """
    out = 0
    out += poly[0]
    out += poly[1]*y + poly[2]*x + poly[3]*z
    out += poly[4]*y*x + poly[5]*y*z +poly[6]*x*z
    out += poly[7]*y*y + poly[8]*x*x + poly[9]*z*z
    out += poly[10]*x*y*z
    out += poly[11]*y*y*y
    out += poly[12]*y*x*x + poly[13]*y*z*z + poly[14]*y*y*x
    out += poly[15]*x*x*x
    out += poly[16]*x*z*z + poly[17]*y*y*z + poly[18]*x*x*z
    out += poly[19]*z*z*z
    return out


def apply_rfm(num, den, x, y, z):
    """
    Evaluates a Rational Function Model (rfm), on a triplet of numbers.

    Args:
        num: list of the 20 coefficients of the numerator
        den: list of the 20 coefficients of the denominator
            All these coefficients are ordered following the RPC convention.
        x, y, z: triplet of floats. They may be numpy arrays of same length.

    Returns:
        the value(s) of the rfm on the input point(s).
    """
    return apply_poly(num, x, y, z) / apply_poly(den, x, y, z)


class RPCModel:
    def __init__(self, d):
        """
        Args:
            d (dict): dictionary read from a geotiff file with
                rasterio.open('/path/to/file.tiff', 'r').tags(ns='RPC')
        """
        self.row_offset = float(d['LINE_OFF'])
        self.col_offset = float(d['SAMP_OFF'])
        self.lat_offset = float(d['LAT_OFF'])
        self.lon_offset = float(d['LONG_OFF'])
        self.alt_offset = float(d['HEIGHT_OFF'])

        self.row_scale = float(d['LINE_SCALE'])
        self.col_scale = float(d['SAMP_SCALE'])
        self.lat_scale = float(d['LAT_SCALE'])
        self.lon_scale = float(d['LONG_SCALE'])
        self.alt_scale = float(d['HEIGHT_SCALE'])

        self.row_num = list(map(float, d['LINE_NUM_COEFF'].split()))
        self.row_den = list(map(float, d['LINE_DEN_COEFF'].split()))
        self.col_num = list(map(float, d['SAMP_NUM_COEFF'].split()))
        self.col_den = list(map(float, d['SAMP_DEN_COEFF'].split()))


    def projection(self, lon, lat, alt):
        nlon = (lon - self.lon_offset) / self.lon_scale
        nlat = (lat - self.lat_offset) / self.lat_scale
        nalt = (alt - self.alt_offset) / self.alt_scale
        col = apply_rfm(self.col_num, self.col_den, nlat, nlon, nalt)
        row = apply_rfm(self.row_num, self.row_den, nlat, nlon, nalt)
        col = col * self.col_scale + self.col_offset
        row = row * self.row_scale + self.row_offset
        return col, row


    def localization(self, col, row, alt, return_normalized=False):

        if not hasattr(self, 'lat_num'):
            return self.localization_iterative(col, row, alt, return_normalized)

        ncol = (col - self.col_offset) / self.col_scale
        nrow = (row - self.row_offset) / self.row_scale
        nalt = (alt - self.alt_offset) / self.alt_scale
        lon = apply_rfm(self.lon_num, self.lon_den, nrow, ncol, nalt)
        lat = apply_rfm(self.lat_num, self.lat_den, nrow, ncol, nalt)
        if not return_normalized:
            lon = lon * self.lon_scale + self.lon_offset
            lat = lat * self.lat_scale + self.lat_offset
        return lon, lat


    def localization_iterative(self, col, row, alt, return_normalized=False):
        """
        Iterative estimation of the localization function (image to ground),
        for a list of image points expressed in image coordinates.

        Args:
            col, row: image coordinates
            alt: altitude (in meters above the ellipsoid) of the corresponding
                3D point
            return_normalized: boolean flag. If true, then return normalized
                coordinates

        Returns:
            lon, lat, alt
        """
        # normalise input image coordinates
        ncol = (col - self.col_offset) / self.col_scale
        nrow = (row - self.row_offset) / self.row_scale
        nalt = (alt - self.alt_offset) / self.alt_scale

        # target point: Xf (f for final)
        Xf = np.vstack([ncol, nrow]).T

        # use 3 corners of the lon, lat domain and project them into the image
        # to get the first estimation of (lon, lat)
        # EPS is 2 for the first iteration, then 0.1.
        lon = -np.ones(len(Xf))
        lat = -np.ones(len(Xf))
        EPS = 2
        x0 = apply_rfm(self.col_num, self.col_den, lat, lon, nalt)
        y0 = apply_rfm(self.row_num, self.row_den, lat, lon, nalt)
        x1 = apply_rfm(self.col_num, self.col_den, lat, lon + EPS, nalt)
        y1 = apply_rfm(self.row_num, self.row_den, lat, lon + EPS, nalt)
        x2 = apply_rfm(self.col_num, self.col_den, lat + EPS, lon, nalt)
        y2 = apply_rfm(self.row_num, self.row_den, lat + EPS, lon, nalt)

        n = 0
        while not np.all((x0 - ncol) ** 2 + (y0 - nrow) ** 2 < 1e-18):
            X0 = np.vstack([x0, y0]).T
            X1 = np.vstack([x1, y1]).T
            X2 = np.vstack([x2, y2]).T
            e1 = X1 - X0
            e2 = X2 - X0
            u  = Xf - X0

            # project u on the base (e1, e2): u = a1*e1 + a2*e2
            # the exact computation is given by:
            #   M = np.vstack((e1, e2)).T
            #   a = np.dot(np.linalg.inv(M), u)
            # but I don't know how to vectorize this.
            # Assuming that e1 and e2 are orthogonal, a1 is given by
            # <u, e1> / <e1, e1>
            num = np.sum(np.multiply(u, e1), axis=1)
            den = np.sum(np.multiply(e1, e1), axis=1)
            a1 = np.divide(num, den)

            num = np.sum(np.multiply(u, e2), axis=1)
            den = np.sum(np.multiply(e2, e2), axis=1)
            a2 = np.divide(num, den)

            # use the coefficients a1, a2 to compute an approximation of the
            # point on the gound which in turn will give us the new X0
            lon += a1 * EPS
            lat += a2 * EPS

            # update X0, X1 and X2
            EPS = .1
            x0 = apply_rfm(self.col_num, self.col_den, lat, lon, nalt)
            y0 = apply_rfm(self.row_num, self.row_den, lat, lon, nalt)
            x1 = apply_rfm(self.col_num, self.col_den, lat, lon + EPS, nalt)
            y1 = apply_rfm(self.row_num, self.row_den, lat, lon + EPS, nalt)
            x2 = apply_rfm(self.col_num, self.col_den, lat + EPS, lon, nalt)
            y2 = apply_rfm(self.row_num, self.row_den, lat + EPS, lon, nalt)
            #n += 1

        #print('localization_iterative: %d iterations' % n)

        if not return_normalized:
            lon = lon * self.lon_scale + self.lon_offset
            lat = lat * self.lat_scale + self.lat_offset

        if np.size(lon) == 1 and np.size(lat) == 1:
            return lon[0], lat[0]
        else:
            return lon, lat


    def incidence_angles(self, lon, lat, z):
        """
        Compute the local incidence angles (zenith and azimuth).

        Args:
            self (rpc_model.RPCModel): camera model
            lon, lat, z (floats): longitude, latitude and altitude

        Return:
            zenith (float in [0, 90]): angle wrt the vertical, in degrees
            azimuth (float in [0, 360]): angle wrt to the north, clockwise, in degrees
        """
        # project the input 3D point in the image
        row, col = np.array(self.projection(lon, lat, z))

        # localize it with two different altitudes
        s = 100  # scale factor, in meters
        lon0, lat0 = self.localization(row, col, z + 0*s)
        lon1, lat1 = self.localization(row, col, z + 1*s)

        # convert to UTM
        zone_number = utm.conversion.latlon_to_zone_number(lat, lon)
        x0, y0 = utm.from_latlon(lat0, lon0, force_zone_number=zone_number)[:2]
        x1, y1 = utm.from_latlon(lat1, lon1, force_zone_number=zone_number)[:2]

        # compute local satellite incidence direction
        p0 = np.array([x0, y0, z + 0*s])
        p1 = np.array([x1, y1, z + 1*s])
        satellite_direction = (p1 - p0) / np.linalg.norm(p1 - p0)

        # return incidence angles
        zenith = np.degrees(np.arccos(satellite_direction @ [0, 0, 1]))
        azimut = np.degrees(np.mod(np.pi/2 - np.angle(np.complex(*satellite_direction[:2])), 2*np.pi))
        return zenith, azimut


    def __repr__(self):
        return """
    # Projection function coefficients
      col_num = {}
      col_den = {}
      row_num = {}
      row_den = {}

    # Offsets and Scales
      row_offset = {}
      col_offset = {}
      lat_offset = {}
      lon_offset = {}
      alt_offset = {}
      row_scale = {}
      col_scale = {}
      lat_scale = {}
      lon_scale = {}
      alt_scale = {}""".format(' '.join(['{: .4f}'.format(x) for x in self.col_num]),
                               ' '.join(['{: .4f}'.format(x) for x in self.col_den]),
                               ' '.join(['{: .4f}'.format(x) for x in self.row_num]),
                               ' '.join(['{: .4f}'.format(x) for x in self.row_den]),
                               self.row_offset,
                               self.col_offset,
                               self.lat_offset,
                               self.lon_offset,
                               self.alt_offset,
                               self.row_scale,
                               self.col_scale,
                               self.lat_scale,
                               self.lon_scale,
                               self.alt_scale)


    def write_to_file(self, path):
        """
        Write RPC coefficients to a txt file in IKONOS txt format.

        Args:
            path (str): path to the output txt file
        """
        with open(path, 'w') as f:

            # scale and offset
            f.write('LINE_OFF: {:.12f} pixels\n'.format(self.row_offset))
            f.write('SAMP_OFF: {:.12f} pixels\n'.format(self.col_offset))
            f.write('LAT_OFF: {:.12f} degrees\n'.format(self.lat_offset))
            f.write('LONG_OFF: {:.12f} degrees\n'.format(self.lon_offset))
            f.write('HEIGHT_OFF: {:.12f} meters\n'.format(self.alt_offset))
            f.write('LINE_SCALE: {:.12f} pixels\n'.format(self.row_scale))
            f.write('SAMP_SCALE: {:.12f} pixels\n'.format(self.col_scale))
            f.write('LAT_SCALE: {:.12f} degrees\n'.format(self.lat_scale))
            f.write('LONG_SCALE: {:.12f} degrees\n'.format(self.lon_scale))
            f.write('HEIGHT_SCALE: {:.12f} meters\n'.format(self.alt_scale))

            # projection function coefficients
            for i in range(20):
                f.write('LINE_NUM_COEFF_{:d}: {:.12f}\n'.format(i+1, self.row_num[i]))
            for i in range(20):
                f.write('LINE_DEN_COEFF_{:d}: {:.12f}\n'.format(i+1, self.row_den[i]))
            for i in range(20):
                f.write('SAMP_NUM_COEFF_{:d}: {:.12f}\n'.format(i+1, self.col_num[i]))
            for i in range(20):
                f.write('SAMP_DEN_COEFF_{:d}: {:.12f}\n'.format(i+1, self.col_den[i]))
