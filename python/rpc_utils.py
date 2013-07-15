import subprocess
import numpy as np
import estimation
import geographiclib


def print_distance_between_vectors(u, v, msg):
    """
    print min, max and mean of the coordinates of two vectors difference
    """
    tmp = u - v
    print 'distance on %s: '%(msg), np.min(tmp), np.max(tmp), np.mean(tmp)


def find_corresponding_point(model_a, model_b, x, y, z):
    """
    find corresponding point using projection model, which can be a projection
    matrix or a rpc.
    x,y is the position of a pixel in image a, and z the altitude of the
    corresponding 3D point. This function returns (xp,yp,z), where xp,yp is the
    projection of the 3D point in image b
    """
    t1, t2, t3 = model_a.direct_estimate(x, y, z)
    xp, yp, zp = model_b.inverse_estimate(t1, t2, z)
    return (xp, yp, z)


def compute_height(model_a, model_b, x, y, xp, yp):
    """
    compute the height of a point given its location inside two images, using
    rpc functions or projection model
    """
    h0 = 0
    p1 = np.array([x, y])
    p2 = np.array([xp, yp])
    HSTEP = 1
    for i in range(100):
        tx, ty, tz = find_corresponding_point(model_a, model_b, p1[0], p1[1], h0)
        r0 = np.array([tx,ty])
        tx, ty, tz = find_corresponding_point(model_a, model_b, p1[0], p1[1], h0+HSTEP)
        r1 = np.array([tx,ty])

        a = r1 - r0
        b = p2 - r0
        # implements:   h0inc = dot(a,b) / dot(a,a)
        h0inc = (a[0]*b[0] + a[1]*b[1]) / (a[0]*a[0]+a[1]*a[1])
        # implements:   q = r0 + h0inc * a
        q = r0
        q[0] = q[0] + h0inc * a[0]
        q[1] = q[1] + h0inc * a[1]
        # implements: err = sqrt( dot(q-p2,q-p2) )
        tmp = q-p2
        err =  np.sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1])
        h0 += h0inc*HSTEP
        # implements: if fabs(h0inc) < 0.0001:
        if np.max(np.fabs(h0inc)) < 0.001:
            break

    return (h0,err)


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
    x, y, z = geographiclib.geodetic_to_geocentric_array(lats, lons, alts)

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
        for i in xrange(len(x)):
            v = np.dot(P, [[x[i]],[y[i]],[z[i]],[1]])
            colPROJ[i] = v[0]/v[2]
            linPROJ[i] = v[1]/v[2]

        print 'approximate_rpc_as_projective: (min, max, mean)'
        print_distance_between_vectors(cols, colPROJ, 'cols')
        print_distance_between_vectors(lins, linPROJ, 'rows')

    return P


def geodesic_bounding_box(rpc, x, y, w, h):
    """
    Computes a bounding box on the WGS84 ellipsoid associated to a Pleiades
    image region of interest, through its rpc function.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        x, y, w, h: four integers definig a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.

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

    # round down lon_m, lat_m and round up lon_M, lat_M so they are integer
    # multiples of 3 arcseconds
    lon_m, lon_M = round_updown(lon_m, lon_M, 1.0/1200.0)
    lat_m, lat_M = round_updown(lat_m, lat_M, 1.0/1200.0)

    # compute the samples
    srtm_step = 1.0/1200 # 6000x6000 samples in a tile of 5x5 degrees
    lons = np.arange(lon_m, lon_M, srtm_step)
    lats = np.arange(lat_m, lat_M, srtm_step)

    # put all the samples in an array. There should a more pythonic way to do
    # this
    out = np.zeros((len(lons)*len(lats), 2))
    for i in xrange(len(lons)):
        for j in xrange(len(lats)):
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

def run_binary_on_list_of_points(points, binary):
    """
    Runs a binary that reads its input on stdin.

    Args:
        points: numpy array containing all the input points, one per line
        binary: path to the binary. It is supposed to write one output value on
            stdout for each input point

    Returns:
        a numpy array containing all the output values.
    """
    # run the binary
    np.savetxt('/tmp/pts', points, '%.18f')
    p1 = subprocess.Popen(['cat', '/tmp/pts'], stdout = subprocess.PIPE)
    p2 = subprocess.Popen([binary], stdin = p1.stdout, stdout =
            subprocess.PIPE)

    # recover output values
    out = np.zeros(len(points))
    for i in range(len(points)):
        line = p2.stdout.readline()
        out[i] = float(line.split()[0])

    return out


def altitude_range_coarse(rpc):
    """
    Computes a coarse altitude range using the RPC informations only.

    Args:
        rpc: instance of the rpc_model.RPCModel class

    Returns:
        the altitude validity range of the RPC.
    """
    m = rpc.altOff - rpc.altScale
    M = rpc.altOff + rpc.altScale
    return m, M


def altitude_range(rpc, x, y, w, h):
    """
    Computes an altitude range using SRTM data.

    Args:
        rpc: instance of the rpc_model.RPCModel class
        x, y, w, h: four integers definig a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.

    Returns:
        lower and upper bounds on the altitude of the world points that are
        imaged by the RPC projection function in the provided ROI. To compute
        these bounds, we use SRTM data. The altitudes are computed with respect
        to the WGS84 reference ellipsoid.
    """

    # find bounding box on the ellipsoid (in geodesic coordinates)
    lon_m, lon_M, lat_m, lat_M = geodesic_bounding_box(rpc, x, y, w, h)

    # if bounding box is out of srtm domain, return coarse altitude estimation
    if (lat_m < -60 or lat_M > 60):
        return altitude_range_coarse(rpc)

    # sample the bounding box with regular step of 3 arcseconds (srtm
    # resolution)
    ellipsoid_points = sample_bounding_box(lon_m, lon_M, lat_m, lat_M)

    # compute srtm height on all these points
    srtm = run_binary_on_list_of_points(ellipsoid_points, 'srtm4')

    # srtm data may contain 'nan' values (meaning no data is available there).
    # These points are most likely water (sea) and thus their height with
    # respect to geoid is 0. But for safety we prefer to give up the precise
    # altitude estimation in these cases and use the coarse one.
    if np.isnan(np.sum(srtm)):
        return altitude_range_coarse(rpc)

    # offset srtm heights with the geoid - ellipsoid difference
    geoid = run_binary_on_list_of_points(ellipsoid_points, 'GeoidEval')
    h = geoid + srtm

    # extract extrema (and add a +-100m security margin)
    h_m = -100 + np.round(h.min())
    h_M =  100 + np.round(h.max())

    return h_m, h_M


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
        x, y, w, h: four integers definig a rectangular region of interest
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
    lon, lat, alt = rpc.direct_estimate(col, row, alt)
    return lon, lat, alt


def matches_from_rpc(rpc1, rpc2, x, y, w, h, n):
    """
    Uses RPC functions to generate matches between two Pleiades images.

    Args:
        rpc1, rpc2: two instances of the rpc_model.RPCModel class
        x, y, w, h: four integers definig a rectangular region of interest
            (ROI) in the first view. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle. In the first view, the matches
            will be located in that ROI.
        n: cube root of the desired number of matches.

    Returns:
        an array of matches, one per line, expressed as x1, y1, x2, y2.
    """
    m, M = altitude_range(rpc1, x, y, w, h)
    lon, lat, alt = ground_control_points(rpc1, x, y, w, h, m, M, n)
    x1, y1, h1 = rpc1.inverse_estimate(lon, lat, alt)
    x2, y2, h2 = rpc2.inverse_estimate(lon, lat, alt)

    return np.vstack([x1, y1, x2, y2]).T


def world_to_image_correspondences_from_rpc(rpc, x, y, w, h, n):
    """
    Uses RPC functions to generate a set of world to image correspondences.

    Args:
        rpc: an instance of the rpc_model.RPCModel class
        x, y, w, h: four integers definig a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle. The correspondences
            will be located in that ROI.
        n: cube root of the desired number of correspondences.

    Returns:
        an array of correspondences, one per line, expressed as X, Y, Z, x, y.
    """
    m, M = altitude_range(rpc, x, y, w, h)
    lon, lat, alt = ground_control_points(rpc, x, y, w, h, m, M, n)
    x, y, h = rpc.inverse_estimate(lon, lat, alt)
    X, Y, Z = geographiclib.geodetic_to_geocentric_array(lat, lon, alt)

    return np.vstack([X, Y, Z]).T, np.vstack([x, y]).T
