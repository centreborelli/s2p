import numpy as np
import geographiclib


class ProjModel:
    def __init__(self, P):
        self.P = P

    def back_project(self, col, lin):
        """
        back-projects pixel (col, lin) to ray, according to Hartley and
        Zisserman "Multiple View Geometry" formula 6.14 (p.162). Warning: this
        formula is valid for finite cameras only. For general projective
        cameras, use formula 6.13.
        This function returns the 3D inhomogeneous coordinates of the direction
        of the ray, together with the camera center (see Zisserman's book)
        """
        M = self.P[:,:3]
        inv_M = np.linalg.inv(M)
        x = np.array([col,lin,1])
        p4 = self.P[:,3]
        return inv_M.dot(x), -inv_M.dot(p4)


    def direct_estimate(self, col, lin, alt):
        """
        computes the (lon,lat,alt) coordinates of the 3D point whose altitude
        above the WGS84 geoid is alt (in meters), and whose projection on the
        image is located in (col,lin)
        """
        # estimate cartesian coordinates of the 3D point
        # at the end of the for loop below, the coordinates are given by t*r+c
        t = 0
        step = 1000
        r,c = self.back_project(col,lin)
        for i in range(0,100):
            p0 = t*r + c
            p1 = p0 + step*r

            tmp_lat, tmp_lon, alt0 = geographiclib.geocentric_to_geodetic(
                                                      p0[0], p0[1], p0[2])
            tmp_lat, tmp_lon, alt1 = geographiclib.geocentric_to_geodetic(
                                                      p1[0], p1[1], p1[2])

            t_inc = (alt - alt0)/(alt1 - alt0)
            t += t_inc * step

            if np.max(np.fabs(t_inc)) < 0.001:
               break

        X = t*r+c
        x = X[0]
        y = X[1]
        z = X[2]

        # convert cartesian coordinates (x,y,z) to geodetic (lat,lon,alt)
        lat,lon,alt = geographiclib.geocentric_to_geodetic(x,y,z)
        return lon, lat, alt

    def inverse_estimate(self, lon, lat, alt):
        """
        computes the projection location (col, lin) of the 3D point
        (lon,lat,alt), on the image.  For symmetry, this function returns a
        triplet (col, lin, alt). The value of alt is just a copy of the one
        given in input.
        """
        # convert the geodetic coordinates (lat, lon, alt) to cartesian (x, y, z)
        # WARNING: the convention used by CartConvert is (lat, lon, alt), while
        # we use (lon, lat, alt)
        x, y, z = geographiclib.geodetic_to_geocentric(lat, lon, alt)

        # project the 3D point, represented as homogeneous cartesian coordinates
        p0 = self.P[0, :]
        p1 = self.P[1, :]
        p2 = self.P[2, :]
        col = p0[0] * x + p0[1] * y + p0[2] * z + p0[3];
        lin = p1[0] * x + p1[1] * y + p1[2] * z + p1[3];
        n   = p2[0] * x + p2[1] * y + p2[2] * z + p2[3];
        # homogeneous normalization
        return col/n, lin/n, alt


if __name__ == '__main__':
    camera_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    P = ProjModel(camera_matrix)
    col, lin = 1000, 1600
    alt = 0 # I don't know what to put here
    print 'col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt)
    lon, lat, alt = P.direct_estimate(col, lin, alt)
    print 'lon={lon}, lat={lat}, alt={alt}'.format(lon=lon, lat=lat, alt=alt)
    col, lin, alt = P.inverse_estimate(lon, lat, alt)
    print 'col={col}, lin={lin}, alt={alt}'.format(col=col, lin=lin, alt=alt)
