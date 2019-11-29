import subprocess

import pytest

from s2p import common


def test_run_success():
    """
    Test s2p.common.run() success with Unix "true" utility command.
    """
    common.run("true")


def test_run_error():
    """
    Test s2p.common.run() error run with Unix "false" utility command.
    """
    with pytest.raises(subprocess.CalledProcessError):
        common.run("false")


def test_run_timeout():
    """
    Test s2p.common.run() timeout with Unix "sleep" utility command.
    """
    with pytest.raises(subprocess.TimeoutExpired):
        common.run("sleep 10", timeout=1)
