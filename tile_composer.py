from python import common
from python import fusion
from python import rpc_model
from python import geographiclib
from python import triangulation
from python import pointing_accuracy
from python import global_params
import main_script
import numpy as np
from multiprocessing import Process
from multiprocessing import cpu_count

img_name = 'uy1'
exp_name = 'campo1'

canvas_x = 5500
canvas_y = 13000
canvas_w = 7000
canvas_h = 10000

canvas_x = 5500
canvas_y = 13000
canvas_w = 3500
canvas_h = 3500

tile_w=2000
tile_h=2000
overlap=100

N = (3 * cpu_count()) / 4

def wait_processes(processes, n):
    """
    Wait until processes terminate.

    More precisely, wait until the number of running processes of the input
    list becomes less than the specified number.

    Args:
        processes: list of Process objects
        n: max number of running processes we want
    """
    while len(processes) > n:
        for p in processes:
            if not p.is_alive():
                processes.remove(p)
    return

def process_pair(img_name, exp_name, x, y, w, h, tile_w=1000, tile_h=1000,
        overlap=100, reference_image_id=1, secondary_image_id=2):
    """
    Computes a height map from a Pair of Pleiades images, using tiles.

    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        tile_w, tile_h: dimensions of the tiles
        overlap: width of overlapping bands between tiles
        reference_image_id: id (1, 2 or 3 for the tristereo datasets, and 1 or
            2 for the bistereo datasets) of the image used as the reference
            image of the pair
        secondary_image_id: id of the image used as the secondary image of the
            pair

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """
    # create a directory for the experiment
    #exp_dir = '/tmp/%s' % exp_name
    exp_dir = '/tmp/'

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
        tile_w = z * np.floor(tile_w / z)
        tile_h = z * np.floor(tile_h / z)
        overlap = z * np.floor(overlap / z)

    # compute global pointing correction (on the whole ROI)
    A = pointing_accuracy.compute_correction(img_name, x, y, w, h,
        reference_image_id, secondary_image_id)

    # generate the tiles - parallelized
    processes = []
    print np.arange(x, x + w, tile_w - overlap)
    print np.arange(y, y + h, tile_h - overlap)
    for i in np.arange(x, x + w, tile_w - overlap):
        for j in np.arange(y, y + h, tile_h - overlap):
#            wait_processes(processes, N-1)
            tile_exp = '%s_%d_%d' % (exp_name, i, j)
#            p = Process(target=main_script.process_pair, args=(img_name,
#                tile_exp, i, j, tile_w, tile_h, reference_image_id,
#                secondary_image_id, exp_dir, A))
#            p.start()
#            processes.append(p)
            main_script.process_pair(img_name, tile_exp, i, j, tile_w, tile_h,
                    reference_image_id, secondary_image_id, exp_dir, A)

    # wait for all the processes to terminate
#    wait_processes(processes, 0)

    # tiles composition
    out = '%s/%s_height_full.tif' % (exp_dir, exp_name)
    common.run('plambda zero:%dx%d "nan" > %s' % (w, h, out))
    for i in np.arange(x, x + w, tile_w - overlap):
        for j in np.arange(y, y + h, tile_h - overlap):
            tile_exp = '%s_%d_%d' % (exp_name, i, j)
            height = '%s/%s/height_unrect.tif' % (exp_dir, tile_exp)
            # weird usage of 'crop' with negative coordinates, to do
            # 'anti-crop' (ie pasting a small image onto a bigger image
            # of nans)
            common.run('crop %s %s/%s/height_to_compose.tif %d %d %d %d' % (height,
                exp_dir, tile_exp, (x-i)/z, (y-j)/z, w/z, h/z))
            # the above divisions should give integer results, as x, y, w, h,
            # tile_w, tile_h and overlap have been modified in order to be
            # multiples of the zoom factor.
            # paste the tile onto the full output image
            common.run('veco med %s %s/%s/height_to_compose.tif | iion - %s' % (out,
                exp_dir, tile_exp, out))


    # cleanup
    while common.garbage:
        common.run('rm ' + common.garbage.pop())

    return out


def process_triplet(img_name, exp_name, x=None, y=None, w=None, h=None,
        tile_w=1000, tile_h=1000, overlap=100, reference_image_id=2,
        left_image_id=1, right_image_id=3):
    """
    Computes a height map from three Pleiades images.

    Args:
        img_name: name of the dataset, located in the 'pleiades_data/images'
            directory
        exp_name: string used to identify the experiment
        x, y, w, h: four integers defining the rectangular ROI in the reference
            image. (x, y) is the top-left corner, and (w, h) are the dimensions
            of the rectangle. The ROI may be as big as you want, as it will be
            cutted into small tiles for processing.
        tile_w, tile_h: dimensions of the tiles
        overlap: width of overlapping bands between tiles
        reference_image_id: id (1, 2 or 3) of the image used as the reference
            image of the triplet
        left_image_id: id of the image used as the secondary image of the first
            pair.
        right_image_id: id of the image used as the secondary image of the
            second pair.

    Returns:
        path to the height map, resampled on the grid of the reference image.
    """

    # select ROI
    try:
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)
    except TypeError:
        rpc = 'pleiades_data/rpc/%s/rpc%02d.xml' % (img_name, reference_image_id)
        prev = 'pleiades_data/images/%s/prev%02d.jpg' % (img_name, reference_image_id)
        x, y, w, h = common.get_roi_coordinates(rpc, prev)
        print "ROI x, y, w, h = %d, %d, %d, %d" % (x, y, w, h)

    # process the two pairs
    exp_name_left = '%s_left' % exp_name
    h_left = process_pair(img_name, exp_name_left, x, y, w, h, tile_w, tile_h,
        overlap, reference_image_id, left_image_id)

    exp_name_right = '%s_right' % exp_name
    h_right = process_pair(img_name, exp_name_right, x, y, w, h, tile_w,
        tile_h, overlap, reference_image_id, right_image_id)

    # merge the two height maps
    h = '/tmp/%s_merged_height.tif' % exp_name
    fusion.merge(h_left, h_right, 3, h)

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

    # colorize, then generate point cloud
    tmp_crop = common.image_crop_TIFF(im, x, y, w, h)
    common.image_safe_zoom_fft(tmp_crop, zoom, crop)
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
