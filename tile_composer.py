from python import common
from python import fusion
from python import pointing_accuracy
from python import global_params
import main_script
import numpy as np
from multiprocessing import Process

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

exp_dir = '/tmp'

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
    # if subsampling_factor is > 1, (ie 2, 3, 4... it has to be int) then
    # ensure that the coordinates of the ROI are multiples of the zoom factor,
    # to avoid bad registration of tiles due to rounding problemes.
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
    A = pointing_accuracy.compute_correction(img_name, exp_name, x, y, w, h,
        reference_image_id, secondary_image_id)

    # generate the tiles - parallelized
    processes = []
    for i in np.arange(x, x + w, tile_w - overlap):
        for j in np.arange(y, y + h, tile_h - overlap):
            print i, j
            tile_exp = '%s_%d_%d' % (exp_name, i, j)
            p = Process(target=main_script.process_pair, args=(img_name,
                tile_exp, i, j, tile_w, tile_h, reference_image_id,
                secondary_image_id, A))
            p.start()
            processes.append(p)

    # prepare the tiles for composition
    # reverse the list of launched processes to access the processes with pop
    # method
    processes.reverse()
    for i in np.arange(x, x + w, tile_w - overlap):
        for j in np.arange(y, y + h, tile_h - overlap):
            p = processes.pop()
            p.join()
            tile_exp = '%s_%d_%d' % (exp_name, i, j)
            height = '%s/%s_height_unrect.tif' % (exp_dir, tile_exp)
            # weird usage of 'crop' with negative coordinates, to do
            # 'anti-crop' (ie pasting a small image onto a bigger image
            # of nans)
            common.run('crop %s %s/%s_height_to_compose.tif %d %d %d %d' % (height,
                exp_dir, tile_exp, (x-i)/z, (y-j)/z, w/z, h/z))
            # the above divisions should give integer results, as x, y, w, h,
            # tile_w, tile_h and overlap have been modified in order to be
            # multiples of the zoom factor.

    # compose the tiles
    out = '%s/%s_height_full.tif' % (exp_dir, exp_name)
    common.run('veco med %s/%s_*_height_to_compose.tif | iion - %s' % (exp_dir,
        exp_name, out))

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
    h = '%s/%s_height_merged.tif' % (exp_dir, exp_name)
    fusion.merge(h_left, h_right, 3, h)

    return h
