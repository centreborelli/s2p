# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>


import warnings
import rasterio
import numpy as np

import rpcm
import srtm4

from s2p import geographiclib
from s2p import common

warnings.filterwarnings("ignore", category=rasterio.errors.NotGeoreferencedWarning)


def find_corresponding_point(model_a, model_b, x, y, z):
    """
    Finds corresponding points in the second image, given the heights.

    Arguments:
        model_a, model_b: two instances of the rpcm.RPCModel class, or of

            the projective_model.ProjModel class
        x, y, z: three 1D numpy arrrays, of the same length. x, y are the
        coordinates of pixels in the image, and z contains the altitudes of the
        corresponding 3D point.

    Returns:
        xp, yp, z: three 1D numpy arrrays, of the same length as the input. xp,
            yp contains the coordinates of the projection of the 3D point in image
            b.
    """
    t1, t2 = model_a.localization(x, y, z)
    xp, yp = model_b.projection(t1, t2, z)
    return (xp, yp, z)


def geodesic_bounding_box(rpc, x, y, w, h):
    """
    Computes a bounding box on the WGS84 ellipsoid associated to a Pleiades
    image region of interest, through its rpc function.

    Args:
        rpc: instance of the rpcm.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are
            the dimensions of the rectangle.

    Returns:
        4 geodesic coordinates: the min and max longitudes, and the min and
        max latitudes.
    """
    # compute altitude coarse extrema from rpc data
    m = rpc.alt_offset - rpc.alt_scale
    M = rpc.alt_offset + rpc.alt_scale

    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    x = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    y = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    a = np.array([m, M,   m,   M,   m,   M,   m,   M])

    # compute geodetic coordinates of corresponding world points
    lon, lat = rpc.localization(x, y, a)

    # extract extrema
    # TODO: handle the case where longitudes pass over -180 degrees
    # for latitudes it doesn't matter since for latitudes out of [-60, 60]
    # there is no SRTM data
    return np.min(lon), np.max(lon), np.min(lat), np.max(lat)


def altitude_range_coarse(rpc, scale_factor=1):
    """
    Computes a coarse altitude range using the RPC informations only.

    Args:
        rpc: instance of the rpcm.RPCModel class
        scale_factor: factor by which the scale offset is multiplied

    Returns:
        the altitude validity range of the RPC.
    """
    m = rpc.alt_offset - scale_factor * rpc.alt_scale
    M = rpc.alt_offset + scale_factor * rpc.alt_scale
    return m, M


def min_max_heights_from_bbx(im, lon_m, lon_M, lat_m, lat_M, rpc,
                             exogenous_dem_geoid_mode, rpc_alt_range_scale_factor):
    """
    Compute min, max heights from bounding box

    Args:
        im: path to an image file
        lon_m, lon_M, lat_m, lat_M: bounding box

    Returns:
        hmin, hmax: min, max heights
    """
    # open image
    dataset = rasterio.open(im, 'r')

    # convert lon/lat to im projection
    x_im_proj, y_im_proj = geographiclib.pyproj_transform([lon_m, lon_M],
                                                          [lat_m, lat_M],
                                                          4326,
                                                          dataset.crs.to_epsg())

    # convert im projection to pixel
    pts = []
    pts.append(~dataset.transform * (x_im_proj[0], y_im_proj[0]))
    pts.append(~dataset.transform * (x_im_proj[1], y_im_proj[1]))
    px = [p[0] for p in pts]
    py = [p[1] for p in pts]

    # get footprint
    [px_min, px_max, py_min, py_max] = map(int, [np.amin(px),
                                                 np.amax(px)+1,
                                                 np.amin(py),
                                                 np.amax(py)+1])

    # limits of im extract
    x, y, w, h = px_min, py_min, px_max - px_min + 1, py_max - py_min + 1
    sizey, sizex = dataset.shape
    x0 = np.clip(x, 0, sizex-1)
    y0 = np.clip(y, 0, sizey-1)
    w -= (x0-x)
    h -= (y0-y)
    w = np.clip(w, 0, sizex - 1 - x0)
    h = np.clip(h, 0, sizey - 1 - y0)

    # get value for each pixel
    if (w != 0) and (h != 0):
        array = dataset.read(1, window=((y0, y0 + h), (x0, x0 + w))).astype(float)
        array[array == -32768] = np.nan
        hmin = np.nanmin(array)
        hmax = np.nanmax(array)

        if exogenous_dem_geoid_mode is True:
            offset = geographiclib.geoid_to_ellipsoid((lat_m + lat_M)/2, (lon_m + lon_M)/2, 0)
            hmin += offset
            hmax += offset
        return hmin, hmax
    else:
        print("WARNING: rpc_utils.min_max_heights_from_bbx: access window out of range")
        print("returning coarse range from rpc")
        return altitude_range_coarse(rpc, rpc_alt_range_scale_factor)


def altitude_range(cfg, rpc, x, y, w, h, margin_top=0, margin_bottom=0):
    """
    Computes an altitude range using the exogenous dem.

    Args:
        rpc: instance of the rpcm.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.
        margin_top: margin (in meters) to add to the upper bound of the range
        margin_bottom: margin (usually negative) to add to the lower bound of
            the range

    Returns:
        lower and upper bounds on the altitude of the world points that are
        imaged by the RPC projection function in the provided ROI. To compute
        these bounds, we use exogenous data. The altitudes are computed with respect
        to the WGS84 reference ellipsoid.
    """
    # TODO: iterate the procedure used here to get a finer estimation of the
    # bounding box on the ellipsoid and thus of the altitude range. For flat
    # regions it will not improve much, but for mountainous regions there is a
    # lot to improve.

    # find bounding box on the ellipsoid (in geodesic coordinates)
    lon_m, lon_M, lat_m, lat_M = geodesic_bounding_box(rpc, x, y, w, h)

    # compute heights on this bounding box
    if cfg['exogenous_dem'] is not None:
        rpc_alt_range_scale_factor = cfg['rpc_alt_range_scale_factor']
        exogenous_dem_geoid_mode = cfg['exogenous_dem_geoid_mode']
        h_m, h_M = min_max_heights_from_bbx(cfg['exogenous_dem'],
                                            lon_m, lon_M, lat_m, lat_M, rpc,
                                            exogenous_dem_geoid_mode,
                                            rpc_alt_range_scale_factor)
        h_m += margin_bottom
        h_M += margin_top
    elif cfg['use_srtm']:
        s = 0.001 / 12  # SRTM90 pixel spacing is 0.001 / 12 degrees
        points = [(lon, lat) for lon in np.arange(lon_m, lon_M, s)
                             for lat in np.arange(lat_m, lat_M, s)]
        lons, lats = np.asarray(points).T
        alts = srtm4.srtm4(lons, lats)  # TODO use srtm4 nn interpolation option
        h_m = min(alts) + margin_bottom
        h_M = max(alts) + margin_top
    else:
        h_m, h_M = altitude_range_coarse(rpc, cfg['rpc_alt_range_scale_factor'])

    return h_m, h_M


def utm_zone(rpc, x, y, w, h):
    """
    Compute the UTM zone where the ROI probably falls (or close to its border).

    Args:
        rpc: instance of the rpcm.RPCModel class, or path to a GeoTIFF file
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle.

    Returns:
        a string of the form '18N' or '18S' where 18 is the utm zone
        identificator.
    """
    # read rpc file
    if not isinstance(rpc, rpcm.RPCModel):
        rpc = rpcm.rpc_from_geotiff(rpc)

    # determine lat lon of the center of the roi, assuming median altitude
    lon, lat = rpc.localization(x + .5*w, y + .5*h, rpc.alt_offset)[:2]

    return geographiclib.compute_utm_zone(lon, lat)


def roi_process(rpc, ll_poly, use_srtm=False, exogenous_dem=None,
                exogenous_dem_geoid_mode=True):
    """
    Convert a (lon, lat) polygon into a rectangular bounding box in image space.

    Args:
        rpc (rpcm.RPCModel): camera model
        ll_poly (array): 2D array of shape (n, 2) containing the vertices
            (longitude, latitude) of the polygon
        use_srtm (bool): whether or not to use SRTM DEM to estimate the
            average ground altitude of the ROI.

    Returns:
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the image. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle.
    """
    if use_srtm and (exogenous_dem is not None):
        raise ValueError("use_srtm and exogenous_dem are mutually exclusive")

    # project lon lat vertices into the image
    lon, lat = np.mean(ll_poly, axis=0)
    if exogenous_dem is not None:
        with rasterio.open(exogenous_dem) as src:
            x, y = geographiclib.pyproj_transform(lon, lat, 4326, src.crs.to_epsg())
            z = list(src.sample([(x, y)]))[0][0]
            if exogenous_dem_geoid_mode is True:
                z = geographiclib.geoid_to_ellipsoid(lat, lon, z)
    elif use_srtm:
        z = srtm4.srtm4(lon, lat)
    else:
        z = rpc.alt_offset
    img_pts = rpc.projection(ll_poly[:, 0], ll_poly[:, 1], z)

    # return image roi
    x, y, w, h = common.bounding_box2D(list(zip(*img_pts)))
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

    return cols, rows, alts


def ground_control_points(rpc, x, y, w, h, m, M, n):
    """
    Computes a set of ground control points (GCP), corresponding to RPC data.

    Args:
        rpc: instance of the rpcm.RPCModel class
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
    lon, lat = rpc.localization(col, row, alt)
    return lon, lat, alt


def corresponding_roi(cfg, rpc1, rpc2, x, y, w, h):
    """
    Uses RPC functions to determine the region of im2 associated to the
    specified ROI of im1.

    Args:
        rpc1, rpc2: two instances of the rpcm.RPCModel class, or paths to
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
    if not isinstance(rpc1, rpcm.RPCModel):
        rpc1 = rpcm.RPCModel(rpc1)
    if not isinstance(rpc2, rpcm.RPCModel):
        rpc2 = rpcm.RPCModel(rpc2)
    m, M = altitude_range(cfg, rpc1, x, y, w, h, 0, 0)

    # build an array with vertices of the 3D ROI, obtained as {2D ROI} x [m, M]
    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
    c = np.array([m, M,   m,   M,   m,   M,   m,   M])

    # corresponding points in im2
    xx, yy = find_corresponding_point(rpc1, rpc2, a, b, c)[0:2]

    # return coordinates of the bounding box in im2
    out = common.bounding_box2D(np.vstack([xx, yy]).T)
    return np.round(out)


def matches_from_rpc(cfg, rpc1, rpc2, x, y, w, h, n):
    """
    Uses RPC functions to generate matches between two Pleiades images.

    Args:
        rpc1, rpc2: two instances of the rpcm.RPCModel class
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the first view. (x, y) is the top-left corner, and (w, h)
            are the dimensions of the rectangle. In the first view, the matches
            will be located in that ROI.
        n: cube root of the desired number of matches.

    Returns:
        an array of matches, one per line, expressed as x1, y1, x2, y2.
    """
    m, M = altitude_range(cfg, rpc1, x, y, w, h, 100, -100)
    lon, lat, alt = ground_control_points(rpc1, x, y, w, h, m, M, n)
    x1, y1 = rpc1.projection(lon, lat, alt)
    x2, y2 = rpc2.projection(lon, lat, alt)

    return np.vstack([x1, y1, x2, y2]).T


def alt_to_disp(rpc1, rpc2, x, y, alt, H1, H2, A=None):
    """
    Converts an altitude into a disparity.

    Args:
        rpc1: an instance of the rpcm.RPCModel class for the reference
            image
        rpc2: an instance of the rpcm.RPCModel class for the secondary
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


def exogenous_disp_range_estimation(cfg, rpc1, rpc2, x, y, w, h, H1, H2,
                                    A=None, margin_top=0, margin_bottom=0):
    """
    Args:
        rpc1: an instance of the rpcm.RPCModel class for the reference
            image
        rpc2: an instance of the rpcm.RPCModel class for the secondary
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
        ROI with the min/max altitude assumptions given by the exogenous dem. The
        disparity is made horizontal thanks to the two rectifying homographies
        H1 and H2.
    """
    m, M = altitude_range(cfg, rpc1, x, y, w, h, margin_top, margin_bottom)
    return altitude_range_to_disp_range(m, M, rpc1, rpc2, x, y, w, h, H1, H2,
                                        A, margin_top, margin_bottom)


def altitude_range_to_disp_range(m, M, rpc1, rpc2, x, y, w, h, H1, H2, A=None,
                                 margin_top=0, margin_bottom=0):
    """
    Args:
        m: min altitude over the tile
        M: max altitude over the tile
        rpc1: instance of the rpcm.RPCModel class for the reference image
        rpc2: instance of the rpcm.RPCModel class for the secondary image
        x, y, w, h: four integers defining a rectangular region of interest
            (ROI) in the reference image. (x, y) is the top-left corner, and
            (w, h) are the dimensions of the rectangle.
        H1, H2: rectifying homographies
        A (optional): pointing correction matrix

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


def gsd_from_rpc(rpc, z=0):
    """
    Compute an image ground sampling distance from its RPC camera model.

    Args:
        rpc (rpcm.RPCModel): camera model
        z (float, optional): ground elevation

    Returns:
        float (meters per pixel)
    """
    # we assume that RPC col/row offset match the coordinates of image center
    c = rpc.col_offset
    r = rpc.row_offset

    a = geographiclib.lonlat_to_geocentric(*rpc.localization(c+0, r, z), alt=z)
    b = geographiclib.lonlat_to_geocentric(*rpc.localization(c+1, r, z), alt=z)
    return np.linalg.norm(np.asarray(b) - np.asarray(a))
