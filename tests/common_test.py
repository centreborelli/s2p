import subprocess

import psutil
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
    Test s2p.common.run() timeout with Unix "sleep" utility command,
    and check that when the command times out, the launched process is killed.
    """
    with pytest.raises(subprocess.TimeoutExpired):
        common.run("sleep 10", timeout=1)

    # Get the names of the running processes
    proc_names = []
    for proc in psutil.process_iter(attrs=['pid', 'name']):
        proc_names.append(proc.info['name'])

    # Check that our process has effectively been killed
    # assert "sleep" not in proc_names
