# s2p (Satellite Stereo Pipeline) testing module
# Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

import os

import s2p


def data_path(p):
    """
    Build an absolute data path from an input test datafile.

    Args:
        p (str): path to the input test data

    Returns:
        str: absolute path to that data file
    """
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(here, "data", p)
