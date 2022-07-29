import subprocess
import time

import pytest

from s2p import parallel


def raise_exception(t, e):
    """
    Wait for t seconds then raise exception e.
    """
    time.sleep(t)
    raise e


def test_launch_calls_error():
    """
    Run several calls to an erroring function through a multiprocessing.Pool.
    """
    with pytest.raises(subprocess.CalledProcessError):
        parallel.launch_calls(
            raise_exception,
            [1, 1, 1, 1],
            2,
            subprocess.CalledProcessError(1, "failcmd"),
            tilewise=False,
        )
