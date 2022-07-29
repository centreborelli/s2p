# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>
# Copyright (C) 2015, Julien Michel <julien.michel@cnes.fr>

import numpy as np

import s2p


def rotation_matrix(t):
    R = np.eye(3)
    R[0, 0] =  np.cos(t)
    R[0, 1] =  np.sin(t)
    R[1, 0] = -np.sin(t)
    R[1, 1] =  np.cos(t)
    return R


def similarity_matrix(t, s):
    R = np.eye(3)
    R[0, 0] =  s * np.cos(t)
    R[0, 1] =  s * np.sin(t)
    R[1, 0] = -s * np.sin(t)
    R[1, 1] =  s * np.cos(t)
    return R


def test_affine_transformation():
    """
    """
    x =  np.array([[0, 0], [0, 1], [1, 0], [1, 1]])

    # list of transformations to be tested
    T = np.eye(3)
    I = np.eye(3)
    S = np.eye(3)
    A = np.eye(3)
    translations = []
    isometries = []
    similarities = []
    affinities = []

    for i in range(100):
        translations.append(T)
        isometries.append(I)
        similarities.append(S)
        affinities.append(A)
        T[:2, 2] = np.random.random(2)
        I = rotation_matrix(2*np.pi * np.random.random_sample())
        I[:2, 2] = np.random.random(2)
        S = similarity_matrix(2*np.pi * np.random.random_sample(),
                              np.random.random_sample())
        S[:2, 2] = 100 * np.random.random(2)
        A[:2, :] = np.random.random((2, 3))

    for B in translations + isometries + similarities + affinities:
        xx = s2p.common.points_apply_homography(B, x)
        E = s2p.estimation.affine_transformation(x, xx)
        np.testing.assert_array_almost_equal(E, B)
