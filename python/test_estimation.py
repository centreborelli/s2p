import numpy as np
from numpy.testing import assert_array_almost_equal 
import estimation 

def test_normalize_2d_points():
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    new_pts, T = estimation.normalize_2d_points(pts)
    assert_array_almost_equal(new_pts, [[-1, -1], [1, -1], [1, 1], [-1, 1]])
    assert_array_almost_equal(T, [[2, 0, -1], [0, 2, -1], [0, 0, 1]])

def test_normalize_3d_points():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    new_pts, U = estimation.normalize_3d_points(pts)
    assert_array_almost_equal(U, [[2, 0, 0, -1], [0, 2, 0, -1], [0, 0, 2, -1], [0, 0, 0, 1]])
    assert_array_almost_equal(new_pts, [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                                        [-1, -1,  1], [1, -1,  1], [1, 1,  1], [-1, 1,  1]])

def test_camera_matrix():
def test_fundamental_matrix():
def test_fundamental_matrix_ransac():
def test_loop_zhang():
