from xml.etree.ElementTree import ElementTree

def apply_poly(poly, x, y, z):
    #TODO: use numpy poly3d function. Test if its faster
    """
    Evaluates a 3-variate polynom of degree 3, on a triplet of numbers.

    Args:
        poly: list of the 20 coefficients of the 3-variate degree 3 polynom
        x, y, z: triplet of floats

    Returns:
        the value of the polynom
    """
    return poly[0] +\
           poly[1]*y + poly[2]*x + poly[3]*z +\
           poly[4]*y*x + poly[5]*y*z +poly[6]*x*z +\
           poly[7]*y*y + poly[8]*x*x + poly[9]*z*z +\
           poly[10]*x*y*z +\
           poly[11]*y*y*y +\
           poly[12]*y*x*x + poly[13]*y*z*z + poly[14]*y*y*x +\
           poly[15]*x*x*x +\
           poly[16]*x*z*z + poly[17]*y*y*z + poly[18]*x*x*z +\
           poly[19]*z*z*z

def apply_rfm(num, den, x, y, z):
    return apply_poly(num, x, y, z) / apply_poly(den, x, y, z)

class RPCModel:
    def __init__(self, XMLfile):
        self.read_rpc(XMLfile)

    def read_rpc(self, XMLfile):
        self.filepath = XMLfile
        tree = ElementTree()
        tree.parse(XMLfile)

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

    def direct_estimate(self, col, lin, alt, return_normalized=False):
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

    def inverse_estimate(self, lon, lat, alt):
        cLon = (lon - self.lonOff) / self.lonScale
        cLat = (lat - self.latOff) / self.latScale
        cAlt = (alt - self.altOff) / self.altScale
        cCol = apply_rfm(self.inverseColNum, self.inverseColDen, cLat, cLon, cAlt)
        cLin = apply_rfm(self.inverseLinNum, self.inverseLinDen, cLat, cLon, cAlt)
        col = cCol*self.colScale+self.colOff
        lin = cLin*self.linScale+self.linOff
        return col, lin, alt

    def __repr__(self):
        return '''
    ### Direct Model ###
        directLonNum = {directLonNum}
        directLonDen = {directLonDen}
        directLatNum = {directLatNum}
        directLatDen = {directLatDen}
        directBias   = {directBias}

    ### Inverse Model ###
        inverseColNum = {inverseColNum}
        inverseColDen = {inverseColDen}
        inverseLinNum = {inverseLinNum}
        inverseLinDen = {inverseLinDen}
        inverseBias   = {inverseBias}

    ### Validity Domains ###
        firstLon = {firstLon}
        firstLat = {firstLat}
        lastLon  = {lastLon}
        lastLat  = {lastLat}

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
        directLonNum  = self.directLonNum,
        directLonDen  = self.directLonDen,
        directLatNum  = self.directLatNum,
        directLatDen  = self.directLatDen,
        directBias    = self.directBias,
        inverseColNum = self.inverseColNum,
        inverseColDen = self.inverseColDen,
        inverseLinNum = self.inverseLinNum,
        inverseLinDen = self.inverseLinDen,
        inverseBias   = self.inverseBias,
        firstLon      = self.firstLon,
        firstLat      = self.firstLat,
        lastLon       = self.lastLon,
        lastLat       = self.lastLat,
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
