from python import common
from python import rpc_model
from python import geographiclib
from python import pointing_accuracy
from python import rectification
from python import block_matching
from python import triangulation
from python import tile_composer
from python import fusion
from python import global_params
import multiprocessing
import numpy as np

N = multiprocessing.cpu_count()

def process_pair_single_tile(out_dir, img_name, ref_img_id=1, sec_img_id=2,
        x=None, y=None, w=None, h=None, A_global=None):
    """
    Computes a height map from a Pair of Pleiades images, without tiling

    Args:
        out_dir: path to the output directory
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        ref_img_id: index of the image used as the reference image of the pair
        sec_img_id: id of the image used as the secondary image of the pair
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle.
        A_global (optional): global pointing correction matrix, used for
            triangulation (but not for stereo-rectification)

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # input files
    im1 = 'pleiades_data/images/%s/im%02d.tif' % (img_name, ref_img_id)
    im2 = 'pleiades_data/images/%s/im%02d.tif' % (img_name, sec_img_id)
    rpc1 = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, ref_img_id)
    rpc2 = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, sec_img_id)
    prev1 = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, ref_img_id)

    # output files
    rect1 = '%s/rectified_ref.tif' % (out_dir)
    rect2 = '%s/rectified_sec.tif' % (out_dir)
    disp    = '%s/rectified_disp.pgm'   % (out_dir)
    mask    = '%s/rectified_mask.png'   % (out_dir)
    height  = '%s/rectified_height.tif' % (out_dir)
    rpc_err = '%s/rpc_err.tif'% (out_dir)
    height_unrect  = '%s/height.tif' % (out_dir)
    subsampling = '%s/subsampling.txt' % (out_dir)
    pointing = '%s/pointing_%02d_%02d.txt' % (out_dir, ref_img_id, sec_img_id)

    ## select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        x, y, w, h = common.get_roi_coordinates(rpc1, prev1)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor
    z = global_params.subsampling_factor
    assert(z > 0 and z == np.floor(z))
    if (z != 1):
        x = z * np.floor(x / z)
        y = z * np.floor(y / z)
        w = z * np.ceil(w / z)
        h = z * np.ceil(h / z)

    ## correct pointing error
    A = pointing_accuracy.compute_correction(img_name, x, y, w, h, ref_img_id,
        sec_img_id)

    ## save the subsampling factor and
    # the pointing correction matrix
    np.savetxt(pointing, A)
    np.savetxt(subsampling, np.array([z]))

    # ATTENTION if subsampling_factor is set the rectified images will be
    # smaller, and the homography matrices and disparity range will reflect
    # this fact

    ## rectification
    H1, H2, disp_min, disp_max = rectification.rectify_pair(im1, im2, rpc1,
        rpc2, x, y, w, h, rect1, rect2, A)

    ## block-matching
    block_matching.compute_disparity_map(rect1, rect2, disp, mask,
        global_params.matching_algorithm, disp_min, disp_max)

    ## triangulation
    if A_global is not None:
        A = A_global
    triangulation.compute_height_map(rpc1, rpc2, H1, H2, disp, mask, height,
        rpc_err, A)
    triangulation.transfer_map(height, H1, x, y, w, h, z, height_unrect)

    return height_unrect



def process_pair(out_dir, img_name, ref_img_id=1, sec_img_id=2, x=None, y=None,
        w=None, h=None, tw=None, th=None, ov=None):
    """
    Computes a height map from a Pair of Pleiades images, using tiles.

    Args:
        out_dir: path to the output directory
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        ref_img_id: id (1, 2 or 3) of the image used as the reference image of
            the pair
        sec_img_id: id of the image used as the secondary image of the
            pair
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        tw, th: dimensions of the tiles
        ov: width of overlapping bands between tiles

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problems.
    z = global_params.subsampling_factor
    assert(z > 0 and z == np.floor(z))
    if (z != 1):
        x = z * np.floor(x / z)
        y = z * np.floor(y / z)
        w = z * np.ceil(w / z)
        h = z * np.ceil(h / z)

    # compute global pointing correction (on the whole ROI)
    A = pointing_accuracy.compute_correction(img_name, x, y, w, h,
        ref_img_id, sec_img_id)

    # automatically compute optimal size for tiles
    # TODO: impose the constraint that ntx*nty is inferior to or equal to a
    # multiple of the number of cores
    if tw is None and th is None and ov is None:
        ov = z * np.ceil(50 / z)
        tw = 300
        th = 300
        while (np.ceil((w - ov) / (tw - ov)) - .2 > (w - ov) / (tw - ov)):
            tw += 1
        while (np.ceil((h - ov) / (th - ov)) - .2 > (h - ov) / (th - ov)):
            th += 1
    ntx = np.ceil((w - ov) / (tw - ov))
    nty = np.ceil((h - ov) / (th - ov))
    # ensure that the coordinates of each tile are multiples of the zoom factor
    if (z != 1):
        ov = z * np.floor(ov / z)
        tw = z * np.floor(tw / z)
        th = z * np.floor(th / z)
    print 'tiles size is tw, th = (%d, %d)' % (tw, th)
    print 'number of tiles is %d, %d' % (ntx, nty)
    print 'total number of tiles is %d' % (ntx * nty)

    # process the tiles
    processes = []
    tiles = []
    for j in np.arange(y, y + h - ov, th - ov):
        for i in np.arange(x, x + w - ov, tw - ov):
            common.wait_processes(processes, N-1)
            tile_dir = '%s/tile_%d_%d_%d_%d' % (out_dir, i, j, tw, th)
            tiles.append('%s/height.tif' % tile_dir)
            p = multiprocessing.Process(target=process_pair_single_tile,
                args=(tile_dir, img_name, ref_img_id, sec_img_id, i, j, tw,
                th, A))
            p.start()
            processes.append(p)
#            process_pair_single_tile(tile_dir, img_name, ref_img_id,
#                    sec_img_id, i, j, tw, th, A)

    # wait for all the processes to terminate
    common.wait_processes(processes, 0)

    # tiles composition
    out = '%s/height.tif' % out_dir
    tile_composer.mosaic(out, w, h, ov, tiles)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    return out


def process_triplet(out_dir, img_name, ref_img_id=2, left_img_id=1,
        right_img_id=3, x=None, y=None, w=None, h=None, thresh=3, tile_w=None,
        tile_h=None, overlap=None):
    """
    Computes a height map from three Pleiades images.

    Args:
        out_dir: path to the output directory
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        ref_img_id: id (1, 2 or 3) of the image used as the reference image of
            the triplet
        left_img_id: id of the image used as the secondary image of the first
            pair.
        right_img_id: id of the image used as the secondary image of the second
            pair.
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        thresh: threshold used for the fusion algorithm, in meters.
        tile_w, tile_h: dimensions of the tiles
        overlap: width of overlapping bands between tiles

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # create a directory for the experiment
    common.run('mkdir -p %s' % out_dir)

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        rpc = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
        prev = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, reference_image_id)
        x, y, w, h = common.get_roi_coordinates(rpc, prev)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # process the two pairs
    out_dir_left = '%s/left' % out_dir
    h_left = process_pair(out_dir_left, img_name, ref_img_id, left_img_id, x,
            y, w, h, tile_w, tile_h, overlap)

    out_dir_right = '%s/right' % out_dir
    h_right = process_pair(out_dir_right, img_name, ref_img_id, right_img_id, x,
            y, w, h, tile_w, tile_h, overlap)

    # merge the two height maps
    h = '%s/height_fusion.tif' % out_dir
    fusion.merge(h_left, h_right, thresh, h)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    return h


def generate_cloud(out_dir, img_name, ref_img_id, x, y, w, h, height_map,
        merc_x=None, merc_y=None):
    """
    Args:
        out_dir: output directory. The file cloud.ply will be written there
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        ref_img_id: id (1, 2 or 3) of the image used as the reference
            image. The height map has been resampled on its grid.
        x, y, w, h: four integers defining the rectangular ROI in the original
            panchro image. (x, y) is the top-left corner, and (w, h) are the
            dimensions of the rectangle.
        height_map: path to the height_map, produced by the process_pair or
            process_triplet function
        merc_{x,y}: mercator coordinates of the point we want to use as
            origin in the local coordinate system of the computed cloud
    """

    # needed input files
    rpc = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, ref_img_id)
    im = 'pleiades_data/images/%s/im%02d.tif' % (img_name, ref_img_id)
    im_color = 'pleiades_data/images/%s/im%02d_color.tif' % (img_name, ref_img_id)

    # output files
    common.run('mkdir -p %s' % out_dir)
    crop   = '%s/roi_ref%02d.tif' % (out_dir, ref_img_id)
    crop_color = '%s/roi_color_ref%02d.tif' % (out_dir, ref_img_id)
    cloud   = '%s/cloud.ply'  % (out_dir)

    # read the zoom value
    zoom = global_params.subsampling_factor

    # build the matrix of the zoom + translation transformation
    A = common.matrix_translation(-x, -y)
    f = 1.0/zoom
    Z = np.diag([f, f, 1])
    A = np.dot(Z, A)
    trans = common.tmpfile('.txt')
    np.savetxt(trans, A)

    # compute mercator offset
    if merc_x is None:
        lat = rpc_model.RPCModel(rpc).firstLat
        lon = rpc_model.RPCModel(rpc).firstLon
        merc_x, merc_y = geographiclib.geodetic_to_mercator(lat, lon)

    # crop the ROI and zoom
    if zoom == 1:
        common.image_crop_TIFF(im, x, y, w, h, crop)
    else:
        tmp_crop = common.image_crop_TIFF(im, x, y, w, h)
        common.image_safe_zoom_fft(tmp_crop, zoom, crop)

    # colorize, then generate point cloud
    try:
        with open(im_color):
            triangulation.colorize(crop, im_color, x, y, zoom, crop_color)
    except IOError:
        print 'no color image available for this dataset.'
        crop_color = common.image_qauto(crop)

    triangulation.compute_point_cloud(crop_color, height_map, rpc, trans,
        cloud, merc_x, merc_y)

    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())


if __name__ == '__main__':

    if len(sys.argv) == 2:
      config  = sys.argv[1]
    else:
      print """
      Incorrect syntax, use:
        > %s config.json

        Launches the s2p pipeline. All the parameters, paths to input and
        output files, are defined in the json configuration file.
      """ % sys.argv[0]
      sys.exit(1)

    # parse the json configuration file
    import json
    f = open(config)
    cfg = json.load(f)
    f.close()

    # launch the two functions
    if len(cfg['data']) == 3:
        dem = process_triplet(out_dir, img_name, ref_img_id=2, left_img_id=1,
            right_img_id=3, x, y, w, h, thresh=3)
    generate_cloud(out_dir, img_name, ref_img_id, x, y, w, h, dem)
