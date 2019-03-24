# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>


from __future__ import print_function
import copy
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
    def __init__(self, rpc_file):
        self.nan_rpc()
        self.read_rpc(rpc_file)

    def nan_rpc(self):
        self.linOff = np.nan
        self.colOff = np.nan
        self.latOff = np.nan
        self.lonOff = np.nan
        self.altOff = np.nan
        self.linScale = np.nan
        self.colScale = np.nan
        self.latScale = np.nan
        self.lonScale = np.nan
        self.altScale = np.nan
        self.directLonNum = [np.nan] * 20
        self.directLonDen = [np.nan] * 20
        self.directLatNum = [np.nan] * 20
        self.directLatDen = [np.nan] * 20
        self.inverseLinNum = [np.nan] * 20
        self.inverseLinDen = [np.nan] * 20
        self.inverseColNum = [np.nan] * 20
        self.inverseColDen = [np.nan] * 20

    def read_rpc(self, rpc_file):
        self.filepath = rpc_file

        if rpc_file.lower().endswith('xml'):
            tree = ElementTree()
            tree.parse(rpc_file)
            self.tree = tree   # store the xml tree in the object
            self.read_rpc_xml(tree)
        else:
            # we assume that non xml rpc files follow the ikonos convention
            self.read_rpc_ikonos(rpc_file)

    def read_rpc_ikonos(self, rpc_file):
        lines = open(rpc_file).read().split('\n')
        for l in lines:
            ll = l.split()
            if len(ll) > 1: self.add_tag_rpc(ll[0], ll[1])

    def add_tag_rpc(self, tag, val):
        a = tag.split('_')
        if len(a) == 2:
            if a[1] == "OFF:":
                if   a[0] == "LINE":   self.linOff = float(val)
                elif a[0] == "SAMP":   self.colOff = float(val)
                elif a[0] == "LAT":    self.latOff = float(val)
                elif a[0] == "LONG":   self.lonOff = float(val)
                elif a[0] == "HEIGHT": self.altOff = float(val)
            elif a[1] == "SCALE:":
                if   a[0] == "LINE":   self.linScale = float(val)
                elif a[0] == "SAMP":   self.colScale = float(val)
                elif a[0] == "LAT":    self.latScale = float(val)
                elif a[0] == "LONG":   self.lonScale = float(val)
                elif a[0] == "HEIGHT": self.altScale = float(val)

        elif len(a) == 4 and a[2] == "COEFF":
            # remove ':', convert to int and decrease the coeff index
            a[3] = int(a[3][:-1]) - 1
            if a[0] == "LINE":
                if   a[1] == "NUM": self.inverseLinNum[a[3]] = float(val)
                elif a[1] == "DEN": self.inverseLinDen[a[3]] = float(val)
            elif a[0] == "SAMP":
                if   a[1] == "NUM": self.inverseColNum[a[3]] = float(val)
                elif a[1] == "DEN": self.inverseColDen[a[3]] = float(val)

    def read_rpc_xml(self, tree):
        # determine wether it's a pleiades, spot-6 or worldview image
        a = tree.find('Metadata_Identification/METADATA_PROFILE') # PHR_SENSOR
        b = tree.find('IMD/IMAGE/SATID') # WorldView
        if a is not None:
            if a.text in ['PHR_SENSOR', 'S6_SENSOR', 'S7_SENSOR']:
                self.read_rpc_pleiades(tree)
            else:
                print('unknown sensor type')
        elif b is not None:
            if b.text == 'WV02' or b.text == 'WV01' or b.text == 'WV03':
                self.read_rpc_worldview(tree)
            else:
                print('unknown sensor type')

    def parse_coeff(self, element, prefix, indices):
        tab = []
        for x in indices:
            tab.append(float(element.find("%s_%s" % (prefix, str(x))).text))
        return tab

    def read_rpc_pleiades(self, tree):
        # direct model
        d = tree.find('Rational_Function_Model/Global_RFM/Direct_Model')
        self.directLonNum = self.parse_coeff(d, "SAMP_NUM_COEFF", range(1, 21))
        self.directLonDen = self.parse_coeff(d, "SAMP_DEN_COEFF", range(1, 21))
        self.directLatNum = self.parse_coeff(d, "LINE_NUM_COEFF", range(1, 21))
        self.directLatDen = self.parse_coeff(d, "LINE_DEN_COEFF", range(1, 21))
        self.directBias = self.parse_coeff(d, "ERR_BIAS", ['X', 'Y'])
        
        # inverse model
        i = tree.find('Rational_Function_Model/Global_RFM/Inverse_Model')
        self.inverseColNum = self.parse_coeff(i, "SAMP_NUM_COEFF", range(1, 21))
        self.inverseColDen = self.parse_coeff(i, "SAMP_DEN_COEFF", range(1, 21))
        self.inverseLinNum = self.parse_coeff(i, "LINE_NUM_COEFF", range(1, 21))
        self.inverseLinDen = self.parse_coeff(i, "LINE_DEN_COEFF", range(1, 21))
        self.inverseBias = self.parse_coeff(i, "ERR_BIAS", ['ROW', 'COL'])
        
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
        # the -1 in line and column offsets is due to Pleiades RPC convention
        # that states that the top-left pixel of an image has coordinates
        # (1, 1)
        self.linOff   = float(v.find('LINE_OFF').text) - 1
        self.colOff   = float(v.find('SAMP_OFF').text) - 1
        self.latOff   = float(v.find('LAT_OFF').text)
        self.lonOff   = float(v.find('LONG_OFF').text)
        self.altOff   = float(v.find('HEIGHT_OFF').text)
        self.linScale = float(v.find('LINE_SCALE').text)
        self.colScale = float(v.find('SAMP_SCALE').text)
        self.latScale = float(v.find('LAT_SCALE').text)
        self.lonScale = float(v.find('LONG_SCALE').text)
        self.altScale = float(v.find('HEIGHT_SCALE').text)

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

        if np.isnan(self.directLatNum[0]):
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


    def direct_estimate_iterative(self, col, row, alt, return_normalized=False):
        """
        Iterative estimation of direct projection (image to ground), for a
        list (or array) of image points expressed in image coordinates.

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
        cCol = (col - self.colOff) / self.colScale
        cRow = (row - self.linOff) / self.linScale
        cAlt = (alt - self.altOff) / self.altScale

        # target point: Xf (f for final)
        Xf = np.vstack([cCol, cRow]).T

        # use 3 corners of the lon, lat domain and project them into the image
        # to get the first estimation of (lon, lat)
        # EPS is 2 for the first iteration, then 0.1.
        lon = -np.ones(len(Xf))
        lat = -np.ones(len(Xf))
        EPS = 2
        x0 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon, cAlt)
        y0 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon, cAlt)
        x1 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon + EPS, cAlt)
        y1 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon + EPS, cAlt)
        x2 = apply_rfm(self.inverseColNum, self.inverseColDen, lat + EPS, lon, cAlt)
        y2 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat + EPS, lon, cAlt)

        n = 0
        while not np.all((x0 - cCol) ** 2 + (y0 - cRow) ** 2 < 1e-18):
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
            x0 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon, cAlt)
            y0 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon, cAlt)
            x1 = apply_rfm(self.inverseColNum, self.inverseColDen, lat, lon + EPS, cAlt)
            y1 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat, lon + EPS, cAlt)
            x2 = apply_rfm(self.inverseColNum, self.inverseColDen, lat + EPS, lon, cAlt)
            y2 = apply_rfm(self.inverseLinNum, self.inverseLinDen, lat + EPS, lon, cAlt)
            #n += 1

        #print('direct_estimate_iterative: %d iterations' % n)

        if return_normalized:
           return lon, lat, cAlt

        # else denormalize and return
        lon = lon*self.lonScale + self.lonOff
        lat = lat*self.latScale + self.latOff
        return lon, lat, alt


    def __write_pleiades(self, filename):
        """
        Writes a new XML file with the rpc parameters
        If the read was performed on a pleiades RPC
        write can only be done using the pleiades format.
        """
        ## First transfer the coefficients back to the internal xml parsing tree
        tree = copy.deepcopy(self.tree)

        # list concatenation of direct model parameters
        direct = self.directLonNum + self.directLonDen + self.directLatNum \
            + self.directLatDen + self.directBias
        d = tree.find('Rational_Function_Model/Global_RFM/Direct_Model')
        for child,id in zip(d,range(82)):
           child.text = str(direct[id])

        # list concatenation of inverse model parameters
        inverse = self.inverseColNum + self.inverseColDen + self.inverseLinNum \
            + self.inverseLinDen + self.inverseBias
        i = tree.find('Rational_Function_Model/Global_RFM/Inverse_Model')
        for child,id in zip(i,range(82)):
           child.text = str(inverse[id])

        # validity domains
        v = tree.find('Rational_Function_Model/Global_RFM/RFM_Validity')
        vd = v.find('Direct_Model_Validity_Domain')
        vd.find('FIRST_ROW').text = str(self.firstRow)
        vd.find('FIRST_COL').text = str(self.firstCol)
        vd.find('LAST_ROW').text  = str(self.lastRow )
        vd.find('LAST_COL').text  = str(self.lastCol )

        vi = v.find('Inverse_Model_Validity_Domain')
        vi.find('FIRST_LON').text = str(self.firstLon)
        vi.find('FIRST_LAT').text = str(self.firstLat)
        vi.find('LAST_LON').text  = str(self.lastLon )
        vi.find('LAST_LAT').text  = str(self.lastLat )

        # scale and offset
        v.find('LINE_OFF').text     = str(self.linOff  )
        v.find('SAMP_OFF').text     = str(self.colOff  )
        v.find('LAT_OFF').text      = str(self.latOff  )
        v.find('LONG_OFF').text     = str(self.lonOff  )
        v.find('HEIGHT_OFF').text   = str(self.altOff  )
        v.find('LINE_SCALE').text   = str(self.linScale)
        v.find('SAMP_SCALE').text   = str(self.colScale)
        v.find('LAT_SCALE').text    = str(self.latScale)
        v.find('LONG_SCALE').text   = str(self.lonScale)
        v.find('HEIGHT_SCALE').text = str(self.altScale)

        ## Write the XML file!
        tree.write(filename)


    def __write_worldview(self, filename):
        """
        Writes a new XML file with the rpc parameters
        If the read was performed on a worldview RPC
        write can only be done using the worldview format.
        """
        ## First transfer the coefficients back to the internal xml parsing tree
        tree = copy.deepcopy(self.tree)
        v = tree.find('RPB/IMAGE')

        # inverse model parameters
        a = [str(x) for x in self.inverseLinNum]
        b = [str(x) for x in self.inverseLinDen]
        c = [str(x) for x in self.inverseColNum]
        d = [str(x) for x in self.inverseColDen]
        v.find('LINENUMCOEFList/LINENUMCOEF').text = ' '.join(a)
        v.find('LINEDENCOEFList/LINEDENCOEF').text = ' '.join(b)
        v.find('SAMPNUMCOEFList/SAMPNUMCOEF').text = ' '.join(c)
        v.find('SAMPDENCOEFList/SAMPDENCOEF').text = ' '.join(d)

        # scale and offset
        v.find('LINEOFFSET').text   = str(self.linOff)
        v.find('SAMPOFFSET').text   = str(self.colOff)
        v.find('LATOFFSET').text    = str(self.latOff)
        v.find('LONGOFFSET').text   = str(self.lonOff)
        v.find('HEIGHTOFFSET').text = str(self.altOff)
        v.find('LINESCALE').text    = str(self.linScale)
        v.find('SAMPSCALE').text    = str(self.colScale)
        v.find('LATSCALE').text     = str(self.latScale)
        v.find('LONGSCALE').text    = str(self.lonScale)
        v.find('HEIGHTSCALE').text  = str(self.altScale)

        # image dimensions
        tree.find('IMD/NUMROWS').text = str(self.lastRow)
        tree.find('IMD/NUMCOLUMNS').text = str(self.lastCol)

        ## Write the XML file!
        tree.write(filename)


    def __write_ikonos(self, filename):
        """
        Writes a text file with the rpc parameters in the Ikonos format.

        If the read was performed on an Ikonos RPC, write can only be done
        using the Ikonos format.
        """
        f = open(filename, "w")

        # scale and offset
        f.write('LINE_OFF: %.12f pixels\n'     % self.linOff  )
        f.write('SAMP_OFF: %.12f pixels\n'     % self.colOff  )
        f.write('LAT_OFF: %.12f degrees\n'     % self.latOff  )
        f.write('LONG_OFF: %.12f degrees\n'    % self.lonOff  )
        f.write('HEIGHT_OFF: %.12f meters\n'   % self.altOff  )
        f.write('LINE_SCALE: %.12f pixels\n'   % self.linScale)
        f.write('SAMP_SCALE: %.12f pixels\n'   % self.colScale)
        f.write('LAT_SCALE: %.12f degrees\n'   % self.latScale)
        f.write('LONG_SCALE: %.12f degrees\n'  % self.lonScale)
        f.write('HEIGHT_SCALE: %.12f meters\n' % self.altScale)

        # inverse model parameters
        for i in range(20):
            f.write('LINE_NUM_COEFF_%d: %.12e\n' % (i+1, self.inverseLinNum[i]))
        for i in range(20):
            f.write('LINE_DEN_COEFF_%d: %.12e\n' % (i+1, self.inverseLinDen[i]))
        for i in range(20):
            f.write('SAMP_NUM_COEFF_%d: %.12e\n' % (i+1, self.inverseColNum[i]))
        for i in range(20):
            f.write('SAMP_DEN_COEFF_%d: %.12e\n' % (i+1, self.inverseLinDen[i]))
        f.close()

    def write(self, filename):
        """
        Saves an rpc object to a file, choosing the Pleiades/Worldview/Ikonos
        format according to the type of the input rpc object

        Args:
            filename: path to the file
        """
        # distinguish 3 cases: pleiades, worldview or ikonos formats
        if hasattr(self, 'tree') and np.isfinite(self.directLatNum[0]):
            self.__write_pleiades(filename)
        elif hasattr(self, 'tree') and np.isnan(self.directLatNum[0]):
            self.__write_worldview(filename)
        else:
            self.__write_ikonos(filename)


    def __repr__(self):
        return '''
    ### Direct Model ###
        directLatNum = {directLatNum}
        directLatDen = {directLatDen}
        directLonNum = {directLonNum}
        directLonDen = {directLonDen}

    ### Inverse Model ###
        inverseColNum = {inverseColNum}
        inverseColDen = {inverseColDen}
        inverseLinNum = {inverseLinNum}
        inverseLinDen = {inverseLinDen}

    ### Scale and Offsets ###
        linOff   = {linOff}
        colOff   = {colOff}
        latOff   = {latOff}
        lonOff   = {lonOff}
        altOff   = {altOff}
        linScale = {linScale}
        colScale = {colScale}
        latScale = {latScale}
        lonScale = {lonScale}
        altScale = {altScale}'''.format(
        directLatNum  = self.directLatNum,
        directLatDen  = self.directLatDen,
        directLonNum  = self.directLonNum,
        directLonDen  = self.directLonDen,
        inverseColNum = self.inverseColNum,
        inverseColDen = self.inverseColDen,
        inverseLinNum = self.inverseLinNum,
        inverseLinDen = self.inverseLinDen,
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
    print('col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt))
    lon, lat, alt = rpc.direct_estimate(col, lin, alt)
    print('lon={lon}, lat={lat}, alt={alt}'.format(lon=lon, lat=lat, alt=alt))
    col, lin, alt = rpc.inverse_estimate(lon, lat, alt)
    print('col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt))
