# Copyright (C) 2013, Carlo de Franchis <carlodef@gmail.com>
# Copyright (C) 2013, Gabriele Facciolo <gfacciol@gmail.com>

import sys
import numpy as np
from xml.etree.ElementTree import ElementTree


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


# this function was written to use numpy.polynomial.polynomial.polyval3d
# function, instead of our apply_poly function.
def reshape_coefficients_vector(c):
    """
    Transform a 1D array of coefficients of a 3D polynom into a 3D array.

    Args:
        c: 1D array of length 20, containing the coefficients of the
            3-variables polynom of degree 3, ordered with the RPC convention.
    Returns:
        a 4x4x4 ndarray, with at most 20 non-zero entries, containing the
        coefficients of input array.
    """
    out = np.zeros((4, 4, 4))
    out[0, 0, 0] = c[0]
    out[0, 1, 0] = c[1]
    out[1, 0, 0] = c[2]
    out[0, 0, 1] = c[3]
    out[1, 1, 0] = c[4]
    out[0, 1, 1] = c[5]
    out[1, 0, 1] = c[6]
    out[0, 2, 0] = c[7]
    out[2, 0, 0] = c[8]
    out[0, 0, 2] = c[9]
    out[1, 1, 1] = c[10]
    out[0, 3, 0] = c[11]
    out[2, 1, 0] = c[12]
    out[0, 1, 2] = c[13]
    out[1, 2, 0] = c[14]
    out[3, 0, 0] = c[15]
    out[1, 0, 2] = c[16]
    out[0, 2, 1] = c[17]
    out[2, 0, 1] = c[18]
    out[0, 0, 3] = c[19]
    return out

def apply_rfm_numpy(num, den, x, y, z):
    """
    Alternative implementation of apply_rfm, that uses numpy to evaluate
    polynoms.
    """
    c_num = reshape_coefficients_vector(num)
    c_den = reshape_coefficients_vector(den)
    a = np.polynomial.polynomial.polyval3d(x, y, z, c_num)
    b = np.polynomial.polynomial.polyval3d(x, y, z, c_den)
    return a/b

class RPCModel:
    def __init__(self, XMLfile):
        self.read_rpc(XMLfile)

    def read_rpc(self, XMLfile):
        self.filepath = XMLfile
        tree = ElementTree()
        tree.parse(XMLfile)

        # determine wether it's a pleiades or a worldview image
        a = tree.find('Metadata_Identification/METADATA_PROFILE') # PHR_SENSOR
        b = tree.find('IMD/IMAGE/SATID') # WV02
        if a is not None:
            if a.text == 'PHR_SENSOR':
                self.satellite = 'pleiades'
                self.read_rpc_pleiades(tree)
            else:
                print 'unknown sensor type'
        elif b is not None:
            if b.text == 'WV02':
                self.satellite = 'worldview2'
                self.read_rpc_worldview(tree)
            else:
                print 'unknown sensor type'

    def read_rpc_pleiades(self, tree):
        # direct model
        d = tree.find('Rational_Function_Model/Global_RFM/Direct_Model')
        direct = [float(child.text) for child in d]
        self.directLonNum = direct[:20]
        self.directLonDen = direct[20:40]
        self.directLatNum = direct[40:60]
        self.directLatDen = direct[60:80]
        self.directBias   = direct[80:]

        # inverse model
        i = tree.find('Rational_Function_Model/Global_RFM/Inverse_Model')
        inverse = [float(child.text) for child in i]
        self.inverseColNum = inverse[:20]
        self.inverseColDen = inverse[20:40]
        self.inverseLinNum = inverse[40:60]
        self.inverseLinDen = inverse[60:80]
        self.inverseBias   = inverse[80:]

        # validity domains
        v = tree.find('Rational_Function_Model/Global_RFM/RFM_Validity')
        vd = v.find('Direct_Model_Validity_Domain')
        self.firstRow = float(vd.find('FIRST_ROW').text)
        self.firstCol = float(vd.find('FIRST_COL').text)
        self.lastRow  = float(vd.find('LAST_ROW').text)
        self.lastCol  = float(vd.find('LAST_COL').text)

        vi = v.find('Inverse_Model_Validity_Domain')
        self.firstLon = float(vi.find('FIRST_LON').text)
        self.firstLat = float(vi.find('FIRST_LAT').text)
        self.lastLon  = float(vi.find('LAST_LON').text)
        self.lastLat  = float(vi.find('LAST_LAT').text)

        # scale and offset
        self.lonScale = float(v.find('LONG_SCALE').text)
        self.lonOff   = float(v.find('LONG_OFF').text)
        self.latScale = float(v.find('LAT_SCALE').text)
        self.latOff   = float(v.find('LAT_OFF').text)
        self.altScale = float(v.find('HEIGHT_SCALE').text)
        self.altOff   = float(v.find('HEIGHT_OFF').text)
        self.colScale = float(v.find('SAMP_SCALE').text)
        self.colOff   = float(v.find('SAMP_OFF').text)
        self.linScale = float(v.find('LINE_SCALE').text)
        self.linOff   = float(v.find('LINE_OFF').text)

    def read_rpc_worldview(self, tree):
        # inverse model
        im = tree.find('RPB/IMAGE')
        l = im.find('LINENUMCOEFList/LINENUMCOEF')
        self.inverseLinNum = [float(c) for c in l.text.split()]
        l = im.find('LINEDENCOEFList/LINEDENCOEF')
        self.inverseLinDen = [float(c) for c in l.text.split()]
        l = im.find('SAMPNUMCOEFList/SAMPNUMCOEF')
        self.inverseColNum = [float(c) for c in l.text.split()]
        l = im.find('SAMPDENCOEFList/SAMPDENCOEF')
        self.inverseColDen = [float(c) for c in l.text.split()]
        self.inverseBias = float(im.find('ERRBIAS').text)

        # scale and offset
        self.linOff   = float(im.find('LINEOFFSET').text)
        self.colOff   = float(im.find('SAMPOFFSET').text)
        self.latOff   = float(im.find('LATOFFSET').text)
        self.lonOff   = float(im.find('LONGOFFSET').text)
        self.altOff   = float(im.find('HEIGHTOFFSET').text)

        self.linScale = float(im.find('LINESCALE').text)
        self.colScale = float(im.find('SAMPSCALE').text)
        self.latScale = float(im.find('LATSCALE').text)
        self.lonScale = float(im.find('LONGSCALE').text)
        self.altScale = float(im.find('HEIGHTSCALE').text)

        # image dimensions
        self.lastRow = int(tree.find('IMD/NUMROWS').text)
        self.lastCol = int(tree.find('IMD/NUMCOLUMNS').text)


    def inverse_estimate(self, lon, lat, alt):
        cLon = (lon - self.lonOff) / self.lonScale
        cLat = (lat - self.latOff) / self.latScale
        cAlt = (alt - self.altOff) / self.altScale
        cCol = apply_rfm(self.inverseColNum, self.inverseColDen, cLat, cLon, cAlt)
        cLin = apply_rfm(self.inverseLinNum, self.inverseLinDen, cLat, cLon, cAlt)
        col = cCol*self.colScale + self.colOff
        lin = cLin*self.linScale + self.linOff
        return col, lin, alt


    def direct_estimate(self, col, lin, alt, return_normalized=False):

        if self.satellite == 'worldview2':
            return self.direct_estimate_iterative(col, lin, alt, return_normalized)

        cCol = (col - self.colOff) / self.colScale
        cLin = (lin - self.linOff) / self.linScale
        cAlt = (alt - self.altOff) / self.altScale
        cLon = apply_rfm(self.directLonNum, self.directLonDen, cLin, cCol, cAlt)
        cLat = apply_rfm(self.directLatNum, self.directLatDen, cLin, cCol, cAlt)
        lon = cLon*self.lonScale + self.lonOff
        lat = cLat*self.latScale + self.latOff
        if return_normalized:
           return cLon, cLat, cAlt
        return lon, lat, alt


    def direct_estimate_iterative_single_point(self, col, row, alt):
        """
        Iterative estimation of direct projection (image to ground), for a
        single point expressed in normalised coordinates.

        Args:
            col, row: normalised image coordinates (ie between -1 and 1)
            alt: altitude (in meters above the ellipsoid) of the corresponding
                3D point

        Returns:
            lon, lat, alt
        """
        # as a first approximation, lon is col and lat is row. The image
        # rows are supposed to be aligned with parallels (lat) and the image
        # columns are supposed to be aligned with the meridians (lon)
        lon = col
        lat = row
        x0 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon, alt)
        y0 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon, alt)
        # print 'initial clon, clat: ', lon, lat
        # print 'initial X0: ', x0, y0

        EPS = 0.1
        #n = 0
        while np.hypot(x0 - col, y0 - row) * self.colScale > 0.001:
            # debug print
           # sys.stdout.write('target: %02d, %02d\r' % (col, row))
           # sys.stdout.write('current position: %f, %f\r' % (x0, y0))
           # sys.stdout.flush()

            # build a base on the image
            x1 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon + EPS, alt)
            y1 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon + EPS, alt)
            x2 = apply_rfm(self.inverseColNum, self.inverseColDen, lat + EPS, lon, alt)
            y2 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat + EPS, lon, alt)

            X0 = np.array([x0, y0])
            X1 = np.array([x1, y1])
            X2 = np.array([x2, y2])
            Xf = np.array([col, row])

            # project Xf-X0 on the base X1-X0, X2-X0
            a1 = np.dot(Xf - X0, X1 - X0) / np.linalg.norm(X1 - X0)
            a2 = np.dot(Xf - X0, X2 - X0) / np.linalg.norm(X2 - X0)

            # use the coefficients a1, a2 to compute an approximation of the
            # point on the gound which in turn will give us the new X0
            lon += a1 * EPS
            lat += a2 * EPS
            x0 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon, alt)
            y0 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon, alt)
            #n += 1

        #print 'direct_estimate_iterative: %d iterations' % n
        return lon, lat, alt

    def direct_estimate_iterative(self, col, row, alt, return_normalized=False):
        """
        Iterative estimation of direct projection (image to ground), for a
        list of points expressed in real coordinates (not normalized).

        Args:
            col, row: tuples or numpy 1D arrays containing image coordinates
            alt: tuple or numpy 1D array containing altitudes (in meters above
                the ellipsoid) of the corresponding 3D points

        Returns:
            lon, lat, alt
        """
        # convert input args to numpy arrays, in case there is only one point
        if (type(col) is int) or (type(col) is float):
            col = np.array([col])
            row = np.array([row])
            alt = np.array([alt])

        # normalise input image coordinates
        cCol = (col - self.colOff) / self.colScale
        cRow = (row - self.linOff) / self.linScale
        cAlt = (alt - self.altOff) / self.altScale

        # iterate over the list of points
        n = len(col)
        cLon = np.zeros(n)
        cLat = np.zeros(n)
        for i in range(n):
            tmp = self.direct_estimate_iterative_single_point(cCol[i], cRow[i],
                    cAlt[i])
            cLon[i] = tmp[0]
            cLat[i] = tmp[1]

        # denormalize and return
        lon = cLon*self.lonScale + self.lonOff
        lat = cLat*self.latScale + self.latOff
        return lon, lat, alt




    def __repr__(self):
        return '''
    ### Inverse Model ###
        inverseColNum = {inverseColNum}
        inverseColDen = {inverseColDen}
        inverseLinNum = {inverseLinNum}
        inverseLinDen = {inverseLinDen}
        inverseBias   = {inverseBias}

    ### Scale and Offsets ###
        lonScale = {lonScale}
        lonOff   = {lonOff}
        latScale = {latScale}
        latOff   = {latOff}
        altScale = {altScale}
        altOff   = {altOff}
        colScale = {colScale}
        colOff   = {colOff}
        linScale = {linScale}
        linOff   = {linOff}'''.format(
        inverseColNum = self.inverseColNum,
        inverseColDen = self.inverseColDen,
        inverseLinNum = self.inverseLinNum,
        inverseLinDen = self.inverseLinDen,
        inverseBias   = self.inverseBias,
        lonScale      = self.lonScale,
        lonOff        = self.lonOff,
        latScale      = self.latScale,
        latOff        = self.latOff,
        altScale      = self.altScale,
        altOff        = self.altOff,
        colScale      = self.colScale,
        colOff        = self.colOff,
        linScale      = self.linScale,
        linOff        = self.linOff)


if __name__ == '__main__':
    # test on the first haiti image
    rpc = RPCModel('../pleiades_data/rpc/haiti/rpc01.xml')
    col, lin = 20000, 8000
    alt = 90
    print 'col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt)
    lon, lat, alt = rpc.direct_estimate(col, lin, alt)
    print 'lon={lon}, lat={lat}, alt={alt}'.format(lon=lon, lat=lat, alt=alt)
    col, lin, alt = rpc.inverse_estimate(lon, lat, alt)
    print 'col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt)
