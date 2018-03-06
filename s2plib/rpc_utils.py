# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>


from __future__ import print_function
import bs4
import utm
import datetime
import numpy as np

from s2plib import estimation
from s2plib import geographiclib
from s2plib import common
from s2plib import rpc_model
from s2plib.config import cfg

def print_distance_between_vectors(u, v, msg):
    """
    print min, max and mean of the coordinates of two vectors difference
    """
    tmp = u - v
    print('distance on %s: '%(msg), np.min(tmp), np.max(tmp), np.mean(tmp))


def find_corresponding_point(model_a, model_b, x, y, z):
    """
    Finds corresponding points in the second image, given the heights.

    Arguments:
        model_a, model_b: two instances of the rpc_model.RPCModel class, or of
            the projective_model.ProjModel class
        x, y, z: three 1D numpy arrrays, of the same length. x, y are the
        coordinates of pixels in the image, and z contains the altitudes of the
        corresponding 3D point.

    Returns:
        xp, yp, z: three 1D numpy arrrays, of the same length as the input. xp,
            yp contains the coordinates of the projection of the 3D point in image
            b.
    """
    t1, t2, t3 = model_a.direct_estimate(x, y, z)
    xp, yp, zp = model_b.inverse_estimate(t1, t2, z)
    return (xp, yp, z)



def compute_height(model_a, model_b, x1, y1, x2, y2):
    """
    Computes the height of a point given its location inside two images.

    Arguments:
        model_a, model_b: two instances of the rpc_model.RPCModel class, or of
            the projective_model.ProjModel class
        x1, y1: two 1D numpy arrrays, of the same length, containing the
            coordinates of points in the first image.
        x2, y2: two 2D numpy arrrays, of the same length, containing the
            coordinates of points in the second image.

    Returns:
        a 1D numpy array containing the list of computed heights.
    """
    n = len(x1)
    h0 = np.zeros(n)
    h0_inc = h0
    p2 = np.vstack([x2, y2]).T
    HSTEP = 1
    err = np.zeros(n)

    for i in range(100):
        tx, ty, tz = find_corresponding_point(model_a, model_b, x1, y1, h0)
        r0 = np.vstack([tx,ty]).T
        tx, ty, tz = find_corresponding_point(model_a, model_b, x1, y1, h0+HSTEP)
        r1 = np.vstack([tx,ty]).T
        a = r1 - r0
        b = p2 - r0
        # implements: h0_inc = dot(a,b) / dot(a,a)
        # For some reason, the formulation below causes massive memory leaks on
        # some systems.
        # h0_inc = np.divide(np.diag(np.dot(a, b.T)), np.diag(np.dot(a, a.T)))
        # Replacing with the equivalent:
        diagabdot = np.multiply(a[:, 0], b[:, 0]) + np.multiply(a[:, 1], b[:, 1])
        diagaadot = np.multiply(a[:, 0], a[:, 0]) + np.multiply(a[:, 1], a[:, 1])
        h0_inc = np.divide(diagabdot, diagaadot)
#        if np.any(np.isnan(h0_inc)):
#            print(x1, y1, x2, y2)
#            print(a)
#            return h0, h0*0
        # implements:   q = r0 + h0_inc * a
        q = r0 + np.dot(np.diag(h0_inc), a)
        # implements: err = sqrt(dot(q-p2, q-p2))
        tmp = q-p2
        err =  np.sqrt(np.multiply(tmp[:, 0], tmp[:, 0]) + np.multiply(tmp[:, 1], tmp[:, 1]))
#       print(np.arctan2(tmp[:, 1], tmp[:, 0])) # for debug
#       print(err) # for debug
        h0 = np.add(h0, h0_inc*HSTEP)
        # implements: if fabs(h0_inc) < 0.0001:
        if np.max(np.fabs(h0_inc)) < 0.001:
            break

    return h0, err


def approximate_rpc_as_projective(rpc_model, col_range, lin_range, alt_range,
        verbose=False):
    """
    Returns a least-square approximation of the RPC functions as a projection
    matrix. The approximation is optimized on a sampling of the 3D region
    defined by the altitudes in alt_range and the image tile defined by
    col_range and lin_range.
    """
    ### step 1: generate cartesian coordinates of 3d points used to fit the
    ###         best projection matrix
    # get mesh points and convert them to geodetic then to geocentric
    # coordinates
    cols, lins, alts = generate_point_mesh(col_range, lin_range, alt_range)
    lons, lats, alts = rpc_model.direct_estimate(cols, lins, alts)
    x, y, z = geographiclib.geodetic_to_geocentric(lats, lons, alts)

    ### step 2: estimate the camera projection matrix from corresponding
    # 3-space and image entities
    world_points = np.vstack([x, y, z]).T
    image_points = np.vstack([cols, lins]).T
    P = estimation.camera_matrix(world_points, image_points)

    ### step 3: for debug, test the approximation error
    if verbose:
        # compute the projection error made by the computed matrix P, on the
        # used learning points
        colPROJ = np.zeros(len(x))
        linPROJ = np.zeros(len(x))
        for i in range(len(x)):
            v = np.dot(P, [[x[i]],[y[i]],[z[i]],[1]])
            colPROJ[i] = v[0]/v[2]
            linPROJ[i] = v[1]/v[2]

        print('approximate_rpc_as_projective: (min, max, mean)')
        print_distance_between_vectors(cols, colPROJ, 'cols')
        print_distance_between_vectors(lins, linPROJ, 'rows')

    return P


def geodesic_bounding_box(rpc, x, y, w, h):
    """
    Computes a bounding box on the WGS84 ellipsoid associated to a Pleiades
    image region of interest, through its rpc function.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are
            the dimensions of the rectangle.

    Returns:
        4 geodesic coordinates: the min and max longitudes, and the min and
        max latitudes.
    """
    # compute altitude coarse extrema from rpc data
    m = rpc.altOff - rpc.altScale
    M = rpc.altOff + rpc.altScale

    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    x = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    y = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    a = np.array([m, M,   m,   M,   m,   M,   m,   M])

    # compute geodetic coordinates of corresponding world points
    lon, lat, alt = rpc.direct_estimate(x, y, a)

    # extract extrema
    # TODO: handle the case where longitudes pass over -180 degrees
    # for latitudes it doesn't matter since for latitudes out of [-60, 60]
    # there is no SRTM data
    return np.min(lon), np.max(lon), np.min(lat), np.max(lat)


def sample_bounding_box(lon_m, lon_M, lat_m, lat_M):
    """
    Samples a geodetic "rectangular" region with regularly spaced points.
    The sampling distance is the srtm resolution, ie 3 arcseconds.

    Args:
        lon_m, lon_M: min and max longitudes, between -180 and 180
        lat_m, lat_M: min and max latitudes, between -60 and 60

    Returns:
        a numpy array, of size N x 2, containing the list of sample locations
        in geodetic coordinates.
    """
    # check parameters
    assert lon_m > -180
    assert lon_M < 180
    assert lon_m < lon_M
    assert lat_m > -60
    assert lat_M < 60
    assert lat_m < lat_M

    # width of srtm bin: 6000x6000 samples in a tile of 5x5 degrees, ie 3
    # arcseconds (in degrees)
    srtm_bin = 1.0/1200

    # round down lon_m, lat_m and round up lon_M, lat_M so they are integer
    # multiples of 3 arcseconds
    lon_m, lon_M = round_updown(lon_m, lon_M, srtm_bin)
    lat_m, lat_M = round_updown(lat_m, lat_M, srtm_bin)

    # compute the samples: one in the center of each srtm bin
    lons = np.arange(lon_m, lon_M, srtm_bin) + .5 * srtm_bin
    lats = np.arange(lat_m, lat_M, srtm_bin) + .5 * srtm_bin

    # put all the samples in an array. There should be a more pythonic way to
    # do this
    out = np.zeros((len(lons)*len(lats), 2))
    for i in range(len(lons)):
        for j in range(len(lats)):
            out[i*len(lats)+j, 0] = lons[i]
            out[i*len(lats)+j, 1] = lats[j]

    return out


def round_updown(a, b, q):
    """
    Rounds down (resp. up) a (resp. b) to the closest multiple of q.

    Args:
        a: float value to round down
        b: float value to round up
        q: float value defining the targeted multiples

    Returns:
        the two modified values
    """
    a = q*np.floor(a/q)
    b = q*np.ceil(b/q)
    return a, b


def altitude_range_coarse(rpc, scale_factor=1):
    """
    Computes a coarse altitude range using the RPC informations only.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        scale_factor: factor by which the scale offset is multiplied

    Returns:
        the altitude validity range of the RPC.
    """
    m = rpc.altOff - scale_factor * rpc.altScale
    M = rpc.altOff + scale_factor * rpc.altScale
    return m, M


def altitude_range(rpc, x, y, w, h, margin_top=0, margin_bottom=0):
    """
    Computes an altitude range using SRTM data.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.
        margin_top: margin (in meters) to add to the upper bound of the range
        margin_bottom: margin (usually negative) to add to the lower bound of
            the range

    Returns:
        lower and upper bounds on the altitude of the world points that are
        imaged by the RPC projection function in the provided ROI. To compute
        these bounds, we use SRTM data. The altitudes are computed with respect
        to the WGS84 reference ellipsoid.
    """
    # TODO: iterate the procedure used here to get a finer estimation of the
    # bounding box on the ellipsoid and thus of the altitude range. For flat
    # regions it will not improve much, but for mountainous regions there is a
    # lot to improve.

    # find bounding box on the ellipsoid (in geodesic coordinates)
    lon_m, lon_M, lat_m, lat_M = geodesic_bounding_box(rpc, x, y, w, h)

    # if bounding box is out of srtm domain, return coarse altitude estimation
    if (lat_m < -60 or lat_M > 60 or cfg['disable_srtm']):
        print("WARNING: returning coarse range from rpc")
        return altitude_range_coarse(rpc, cfg['rpc_alt_range_scale_factor'])

    # sample the bounding box with regular step of 3 arcseconds (srtm
    # resolution)
    ellipsoid_points = sample_bounding_box(lon_m, lon_M, lat_m, lat_M)

    # compute srtm height on all these points
    srtm = common.run_binary_on_list_of_points(ellipsoid_points, 'srtm4',
                                               env_var=('SRTM4_CACHE',
                                                        cfg['srtm_dir']))
    h = np.ravel(srtm)

    # srtm data may contain 'nan' values (meaning no data is available there).
    # These points are most likely water (sea) and thus their height with
    # respect to geoid is 0. Thus we replace the nans with 0.
    # TODO: this should not be zero, but the geoid/ellipsoid offset
    srtm[np.isnan(h)] = 0

    # extract extrema (and add a +-100m security margin)
    h_m = np.round(h.min()) + margin_bottom
    h_M = np.round(h.max()) + margin_top

    return h_m, h_M


def utm_zone(rpc, x, y, w, h):
    """
    Compute the UTM zone where the ROI probably falls (or close to its border).

    Args:
        rpc: instance of the rpc_model.RPCModel class, or path to the xml file
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle.

    Returns:
        a string of the form '18N' or '18S' where 18 is the utm zone
        identificator.
    """
    # read rpc file
    if not isinstance(rpc, rpc_model.RPCModel):
        rpc = rpc_model.RPCModel(rpc)

    # determine lat lon of the center of the roi, assuming median altitude
    lon, lat = rpc.direct_estimate(x + .5*w, y + .5*h, rpc.altOff)[:2]

    # compute the utm zone number and add the hemisphere letter
    zone = utm.conversion.latlon_to_zone_number(lat, lon)
    if lat < 0:
        return '%dS' % zone
    else:
        return '%dN' % zone


def utm_roi_to_img_roi(rpc, roi):
    """
    """
    # define utm rectangular box
    x, y, w, h = [roi[k] for k in ['x', 'y', 'w', 'h']]
    box = [(x, y), (x+w, y), (x+w, y+h), (x, y+h)]

    # convert utm to lon/lat
    utm_z = roi['utm_band']
    north = roi['hemisphere'] == 'N'
    box_latlon = [utm.to_latlon(p[0], p[1], utm_z, northern=north) for p in box]

    # project lon/lat vertices into the image
    if not isinstance(rpc, rpc_model.RPCModel):
        rpc = rpc_model.RPCModel(rpc)
    img_pts = [rpc.inverse_estimate(p[1], p[0], rpc.altOff)[:2] for p in
               box_latlon]

    # return image roi
    x, y, w, h = common.bounding_box2D(img_pts)
    return {'x': x, 'y': y, 'w': w, 'h': h}


def kml_roi_process(rpc, kml):
    """
    """
    # extract lon lat from kml
    with open(kml, 'r') as f:
        a = bs4.BeautifulSoup(f, "lxml").find_all('coordinates')[0].text.split()
    ll_bbx = np.array([list(map(float, x.split(','))) for x in a])[:4, :2]

    # save lon lat bounding box to cfg dictionary
    lon_min = min(ll_bbx[:, 0])
    lon_max = max(ll_bbx[:, 0])
    lat_min = min(ll_bbx[:, 1])
    lat_max = max(ll_bbx[:, 1])
    cfg['ll_bbx'] = (lon_min, lon_max, lat_min, lat_max)

    # convert lon lat bbox to utm
    z = utm.conversion.latlon_to_zone_number((lat_min + lat_max) * .5,
                                             (lon_min + lon_max) * .5)
    utm_bbx = np.array([utm.from_latlon(p[1], p[0], force_zone_number=z)[:2] for
                        p in ll_bbx])
    east_min = min(utm_bbx[:, 0])
    east_max = max(utm_bbx[:, 0])
    nort_min = min(utm_bbx[:, 1])
    nort_max = max(utm_bbx[:, 1])
    cfg['utm_bbx'] = (east_min, east_max, nort_min, nort_max)

    # project lon lat vertices into the image
    if not isinstance(rpc, rpc_model.RPCModel):
        rpc = rpc_model.RPCModel(rpc)
    img_pts = [rpc.inverse_estimate(p[0], p[1], rpc.altOff)[:2] for p in ll_bbx]

    # return image roi
    x, y, w, h = common.bounding_box2D(img_pts)
    return {'x': x, 'y': y, 'w': w, 'h': h}


def generate_point_mesh(col_range, row_range, alt_range):
    """
    Generates image coordinates (col, row, alt) of 3D points located on the grid
    defined by col_range and row_range, at uniformly sampled altitudes defined
    by alt_range.
    Args:
        col_range: triplet (col_min, col_max, n_col), where n_col is the
            desired number of samples
        row_range: triplet (row_min, row_max, n_row)
        alt_range: triplet (alt_min, alt_max, n_alt)

    Returns:
        3 lists, containing the col, row and alt coordinates.
    """
    # input points in col, row, alt space
    cols, rows, alts = [np.linspace(v[0], v[1], v[2]) for v in
            [col_range, row_range, alt_range]]

    # make it a kind of meshgrid (but with three components)
    # if cols, rows and alts are lists of length 5, then after this operation
    # they will be lists of length 5x5x5
    cols, rows, alts =\
            (  cols+0*rows[:,np.newaxis]+0*alts[:,np.newaxis,np.newaxis]).reshape(-1),\
            (0*cols+  rows[:,np.newaxis]+0*alts[:,np.newaxis,np.newaxis]).reshape(-1),\
            (0*cols+0*rows[:,np.newaxis]+  alts[:,np.newaxis,np.newaxis]).reshape(-1)

    return (cols, rows, alts)


def ground_control_points(rpc, x, y, w, h, m, M, n):
    """
    Computes a set of ground control points (GCP), corresponding to RPC data.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are
            the dimensions of the rectangle.
        m, M: minimal and maximal altitudes of the ground control points
        n: cube root of the desired number of ground control points.

    Returns:
        a list of world points, given by their geodetic (lon, lat, alt)
        coordinates.
    """
    # points will be sampled in [x, x+w] and [y, y+h]. To avoid always sampling
    # the same four corners with each value of n, we make these intervals a
    # little bit smaller, with a dependence on n.
    col_range = [x+(1.0/(2*n))*w, x+((2*n-1.0)/(2*n))*w, n]
    row_range = [y+(1.0/(2*n))*h, y+((2*n-1.0)/(2*n))*h, n]
    alt_range = [m, M, n]
    col, row, alt = generate_point_mesh(col_range, row_range, alt_range)
    return rpc.direct_estimate(col, row, alt)


def corresponding_roi(rpc1, rpc2, x, y, w, h):
    """
    Uses RPC functions to determine the region of im2 associated to the
    specified ROI of im1.

    Args:
        rpc1, rpc2: two instances of the rpc_model.RPCModel class, or paths to
            the xml files
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the first view. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle.

    Returns:
        four integers defining a ROI in the second view. This ROI is supposed
        to contain the projections of the 3D points that are visible in the
        input ROI.
    """
    # read rpc files
    if not isinstance(rpc1, rpc_model.RPCModel):
        rpc1 = rpc_model.RPCModel(rpc1)
    if not isinstance(rpc2, rpc_model.RPCModel):
        rpc2 = rpc_model.RPCModel(rpc2)
    m, M = altitude_range(rpc1, x, y, w, h, 0, 0)

    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    c = np.array([m, M,   m,   M,   m,   M,   m,   M])

    # corresponding points in im2
    xx, yy = find_corresponding_point(rpc1, rpc2, a, b, c)[0:2]

    # return coordinates of the bounding box in im2
    out = common.bounding_box2D(np.vstack([xx, yy]).T)
    return np.round(out)


def matches_from_rpc(rpc1, rpc2, x, y, w, h, n):
    """
    Uses RPC functions to generate matches between two Pleiades images.

    Args:
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the first view. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle. In the first view, the matches
            will be located in that ROI.
        n: cube root of the desired number of matches.

    Returns:
        an array of matches, one per line, expressed as x1, y1, x2, y2.
    """
    m, M = altitude_range(rpc1, x, y, w, h, 100, -100)
    lon, lat, alt = ground_control_points(rpc1, x, y, w, h, m, M, n)
    x1, y1, h1 = rpc1.inverse_estimate(lon, lat, alt)
    x2, y2, h2 = rpc2.inverse_estimate(lon, lat, alt)

    return np.vstack([x1, y1, x2, y2]).T


def world_to_image_correspondences_from_rpc(rpc, x, y, w, h, n):
    """
    Uses RPC functions to generate a set of world to image correspondences.

    Args:
        rpc: an instance of the rpc_model.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle. The correspondences
            will be located in that ROI.
        n: cube root of the desired number of correspondences.

    Returns:
        an array of correspondences, one per line, expressed as X, Y, Z, x, y.
    """
    m, M = altitude_range(rpc, x, y, w, h, 100, -100)
    lon, lat, alt = ground_control_points(rpc, x, y, w, h, m, M, n)
    x, y, h = rpc.inverse_estimate(lon, lat, alt)
    X, Y, Z = geographiclib.geodetic_to_geocentric(lat, lon, alt)

    return np.vstack([X, Y, Z]).T, np.vstack([x, y]).T


def alt_to_disp(rpc1, rpc2, x, y, alt, H1, H2, A=None):
    """
    Converts an altitude into a disparity.

    Args:
        rpc1: an instance of the rpc_model.RPCModel class for the reference
            image
        rpc2: an instance of the rpc_model.RPCModel class for the secondary
            image
        x, y: coordinates of the point in the reference image
        alt: altitude above the WGS84 ellipsoid (in meters) of the point
        H1, H2: rectifying homographies
        A (optional): pointing correction matrix

    Returns:
        the horizontal disparity of the (x, y) point of im1, assuming that the
        3-space point associated has altitude alt. The disparity is made
        horizontal thanks to the two rectifying homographies H1 and H2.
    """
    xx, yy = find_corresponding_point(rpc1, rpc2, x, y, alt)[0:2]
    p1 = np.vstack([x, y]).T
    p2 = np.vstack([xx, yy]).T

    if A is not None:
        print("rpc_utils.alt_to_disp: applying pointing error correction")
        # correct coordinates of points in im2, according to A
        p2 = common.points_apply_homography(np.linalg.inv(A), p2)

    p1 = common.points_apply_homography(H1, p1)
    p2 = common.points_apply_homography(H2, p2)
    # np.testing.assert_allclose(p1[:, 1], p2[:, 1], atol=0.1)
    disp = p2[:, 0] - p1[:, 0]
    return disp


def srtm_disp_range_estimation(rpc1, rpc2, x, y, w, h, H1, H2, A=None,
        margin_top=0, margin_bottom=0):
    """
    Args:
        rpc1: an instance of the rpc_model.RPCModel class for the reference
            image
        rpc2: an instance of the rpc_model.RPCModel class for the secondary
            image
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the reference image. (x, y) is the top-left corner, and
            (w, h) are the dimensions of the rectangle.
        H1, H2: rectifying homographies
        A (optional): pointing correction matrix
        margin_top: margin (in meters) to add to the upper bound of the range
        margin_bottom: margin (negative) to add to the lower bound of the range

    Returns:
        the min and max horizontal disparity observed on the 4 corners of the
        ROI with the min/max altitude assumptions given by the srtm. The
        disparity is made horizontal thanks to the two rectifying homographies
        H1 and H2.
    """
    m, M = altitude_range(rpc1, x, y, w, h, margin_top, margin_bottom)

    return altitude_range_to_disp_range(m, M, rpc1, rpc2, x, y, w, h, H1, H2,
                                        A, margin_top, margin_bottom)

def altitude_range_to_disp_range(m, M, rpc1, rpc2, x, y, w, h, H1, H2, A=None,
                                 margin_top=0, margin_bottom=0):
    """
    Args:
        m: min altitude over the tile
        M: max altitude over the tile
        rpc1: instance of the rpc_model.RPCModel class for the reference image
        rpc2: instance of the rpc_model.RPCModel class for the secondary image
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the reference image. (x, y) is the top-left corner, and
            (w, h) are the dimensions of the rectangle.
        H1, H2: rectifying homographies
        A (optional): pointing correction matrix
        margin_top: margin (in meters) to add to the upper bound of the range
        margin_bottom: margin (negative) to add to the lower bound of the range

    Returns:
        the min and max horizontal disparity observed on the 4 corners of the
        ROI with the min/max altitude assumptions given as parameters. The
        disparity is made horizontal thanks to the two rectifying homographies
        H1 and H2.
    """
    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    c = np.array([m, M,   m,   M,   m,   M,   m,   M])

    # compute the disparities of these 8 points
    d = alt_to_disp(rpc1, rpc2, a, b, c, H1, H2, A)

    # return min and max disparities
    return np.min(d), np.max(d)

def compute_ms_panchro_offset(dim_pan, dim_ms):
    """
    Computes the offset, in panchro pixels,  between a panchro image and the
    corresponding 4x zoomed ms image.

    Args:
        dim_pan: path to the xml file DIM_*.XML for the panchro image
        dim_ms:  path to the xml file DIM_*.XML for the ms (ie color) image

    Returns:
        (off_col, off_row): the offset to apply to the ms image before fusion
        with the panchro.
    """
    # column offset
    first_col_ms = int(common.grep_xml(dim_ms, "FIRST_COL"))
    first_col_pan = int(common.grep_xml(dim_pan, "FIRST_COL"))
    off_col = 4 * first_col_ms - first_col_pan

    # row offset
    t_e_pan = float(common.grep_xml(dim_pan, "LINE_PERIOD"))
    t_init_ms  = common.grep_xml(dim_ms, "START")
    t_init_pan = common.grep_xml(dim_pan, "START")
    t_ms =  datetime.datetime.strptime(t_init_ms[:26], "%Y-%m-%dT%H:%M:%S.%f")
    t_pan = datetime.datetime.strptime(t_init_pan[:26], "%Y-%m-%dT%H:%M:%S.%f")

    delta_t = 1000 * (t_ms - t_pan)
    off_row = int(total_seconds(delta_t) / t_e_pan)

    #print("t_e_pan: %f" % t_e_pan)
    #print("t_init_ms:  %s" % t_init_ms)
    #print("t_init_pan: %s" % t_init_pan)
    #print(off_col, off_row)

    return off_col, off_row


def total_seconds(td):
    """
    Return the total number of seconds contained in the duration.

    Args:
        td: datetime.timedelta object

    Returns:
        the equivalent time expressed in seconds

    This function implements the timedelta.total_seconds() method available in
    python 2.7, to make the compute_ms_panchro_offset usable even with python
    2.6
    """
    return float((td.microseconds + (td.seconds + td.days * 24 * 3600) *
        10**6)) / 10**6


def crop_corresponding_areas(out_dir, images, roi, zoom=1):
    """
    Crops areas corresponding to the reference ROI in the secondary images.

    Args:
        out_dir:
        images: sequence of dicts containing the paths to input data
        roi: dictionary containing the ROI definition
        zoom: integer zoom out factor
    """
    rpc_ref = images[0]['rpc']
    for i, image in enumerate(images[1:]):
        x, y, w, h = corresponding_roi(rpc_ref, image['rpc'], roi['x'],
                                       roi['y'], roi['w'], roi['h'])
        if zoom == 1:
            common.image_crop_gdal(image['img'], x, y, w, h, '%s/roi_sec_%d.tif' % (out_dir, i))
        else:
            # gdal is used for the zoom because it handles BigTIFF files, and
            # before the zoom out the image may be that big
            tmp = common.image_crop_gdal(image['img'], x, y, w, h)
            common.image_zoom_gdal(tmp, zoom, '%s/roi_sec_%d.tif' % (out_dir, i), w, h)


def rpc_from_geotiff(geotiff_path, outrpcfile='.rpc'):
    """
    extracts the rpc from a geotiff file (including vsicurl)
    """
    import os
    import subprocess

    env = os.environ.copy()
    if geotiff_path.startswith(('http://', 'https://')):
        env['CPL_VSIL_CURL_ALLOWED_EXTENSIONS'] = geotiff_path[-3:]
        path = '/vsicurl/{}'.format(geotiff_path)
    else:
        path = geotiff_path

    f = open(outrpcfile, 'wb')
    x = subprocess.Popen(["gdalinfo", path], stdout=subprocess.PIPE).communicate()[0]
    x = x.splitlines()
    for l in x:

        if (b'SAMP_' not in l) and (b'LINE_' not in l) and (b'HEIGHT_' not in l) and (b'LAT_' not in l) and (b'LONG_' not in l) and (b'MAX_' not in l) and (b'MIN_' not in l):
              continue
        y = l.strip().replace(b'=',b': ')
        if b'COEFF' in y:
              z = y.split(b' ')
              t=1
              for j in z[1:]:
                      f.write(b'%s_%d: %s\n'%(z[0][:-1],t,j))
                      t+=1
        else:
              f.write((y+b'\n'))

    f.close()
    return rpc_model.RPCModel(outrpcfile)
