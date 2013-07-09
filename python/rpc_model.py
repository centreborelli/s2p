from xml.etree.ElementTree import ElementTree

def apply_poly(poly, x, y, z):
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

def apply_rfm( num, den, x, y, z):
    return apply_poly(num, x, y, z) / apply_poly(den, x, y, z)

class RPCModel:
    def __init__(self, XMLfile):
        self.read_rpc(XMLfile)

    def read_rpc(self,XMLfile):
        self.filepath = XMLfile
        tree = ElementTree()
        tree.parse(XMLfile)
        # direct model
        d=tree.find('Rational_Function_Model/Global_RFM/Direct_Model')
        direct=[float(child.text) for child in d]
        self.directLonNum,self.directLonDen,self.directLatNum,self.directLatDen,self.directBias=direct[:20],direct[20:40],direct[40:60],direct[60:80],direct[80:]
        # inverse model
        i=tree.find('Rational_Function_Model/Global_RFM/Inverse_Model')
        inverse=[float(child.text) for child in i]
        self.inverseColNum,self.inverseColDen,self.inverseLinNum,self.inverseLinDen,self.inverseBias=inverse[:20],inverse[20:40],inverse[40:60],inverse[60:80],inverse[80:]
        # validity domains
        v = tree.find('Rational_Function_Model/Global_RFM/RFM_Validity')
        vd=v.find('Direct_Model_Validity_Domain')
        self.firstRow, self.firstCol, self.lastRow, self.lastCol=[float(vd.find(s).text) for s in 'FIRST_ROW', 'FIRST_COL', 'LAST_ROW', 'LAST_COL']
        vi=v.find('Inverse_Model_Validity_Domain')
        self.firstLon, self.firstLat, self.lastLon, self.lastLat= [float(vi.find(s).text) for s in 'FIRST_LON', 'FIRST_LAT', 'LAST_LON', 'LAST_LAT']
        #  scale and offset
        self.lonScale,self.lonOff,self.latScale,self.latOff,self.altScale,self.altOff,self.colScale,self.colOff,self.linScale,self.linOff = [float(v.find(s).text) for s in 'LONG_SCALE', 'LONG_OFF', 'LAT_SCALE', 'LAT_OFF', 'HEIGHT_SCALE', 'HEIGHT_OFF', 'SAMP_SCALE', 'SAMP_OFF', 'LINE_SCALE', 'LINE_OFF']

    def direct_estimate(self, col, lin, alt, return_normalized=False):
        cCol, cLin, cAlt=(col-self.colOff)/self.colScale, (lin-self.linOff)/self.linScale, (alt-self.altOff)/self.altScale
        cLon,cLat=apply_rfm(self.directLonNum, self.directLonDen, cLin, cCol, cAlt), apply_rfm(self.directLatNum, self.directLatDen, cLin, cCol, cAlt),
        lon,lat=cLon*self.lonScale+self.lonOff, cLat*self.latScale+self.latOff
        if return_normalized:
           return cLon,cLat,cAlt
        return lon,lat,alt

    def inverse_estimate(self, lon, lat, alt):
        cLon, cLat, cAlt=(lon-self.lonOff)/self.lonScale, (lat-self.latOff)/self.latScale, (alt-self.altOff)/self.altScale
        cCol,cLin=apply_rfm(self.inverseColNum, self.inverseColDen, cLat, cLon, cAlt), apply_rfm(self.inverseLinNum, self.inverseLinDen, cLat, cLon, cAlt),
        col,lin=cCol*self.colScale+self.colOff, cLin*self.linScale+self.linOff
        return col,lin,alt

    def __repr__(self):
        return '''
    ### Direct Model ###
        directLonNum={directLonNum}
        directLonDen={directLonDen}
        directLatNum={directLatNum}
        directLatDen={directLatDen}
        directBias={directBias}

    ### Inverse Model ###
        inverseColNum={inverseColNum}
        inverseColDen={inverseColDen}
        inverseLinNum={inverseLinNum}
        inverseLinDen={inverseLinDen}
        inverseBias={inverseBias}

    ### Validity Domains ###
        firstLon={firstLon}
        firstLat={firstLat}
        lastLon={lastLon}
        lastLat={lastLat}

    ### Scale and Offsets ###
        lonScale={lonScale}
        lonOff={lonOff}
        latScale={latScale}
        latOff={latOff}
        altScale={altScale}
        altOff={altOff}
        colScale={colScale}
        colOff={colOff}
        linScale={linScale}
        linOff={linOff}'''.format(directLonNum=self.directLonNum,directLonDen=self.directLonDen,directLatNum=self.directLatNum,directLatDen=self.directLatDen,directBias=self.directBias,
                                  inverseColNum=self.inverseColNum,inverseColDen=self.inverseColDen,inverseLinNum=self.inverseLinNum,inverseLinDen=self.inverseLinDen,inverseBias=self.inverseBias,
                                  firstLon=self.firstLon, firstLat=self.firstLat, lastLon=self.lastLon, lastLat=self.lastLat,
                                  lonScale=self.lonScale,lonOff=self.lonOff,latScale=self.latScale,latOff=self.latOff,altScale=self.altScale,altOff=self.altOff,colScale=self.colScale,colOff=self.colOff,linScale=self.linScale,linOff=self.linOff)


if __name__ == '__main__':
    # test on the first melbourne image
    rpc = RPCModel('../testdata/RPC_PHR1A_P_201202250024153_SEN_IPU_20120229_8029-001.XML')
    col, lin = 20000, 8000
    alt = 90 # I don't know what to put here
    print 'col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt)
    lon, lat, alt = rpc.direct_estimate(col, lin, alt)
    print 'lon={lon}, lat={lat}, alt={alt}'.format(lon=lon, lat=lat, alt=alt)
    col, lin, alt = rpc.inverse_estimate(lon, lat, alt)
    print 'col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt)
