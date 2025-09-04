#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Toplevel pytest configuration for Pestifer
"""

from __future__ import annotations

import pytest
import os
import logging

logging.getLogger("pidibble").setLevel(logging.WARNING)
logging.getLogger("ycleptic").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)
try:
    logging.getLogger("numexp").setLevel(logging.WARNING)
except:
    pass

from logging import FileHandler, Formatter, StreamHandler
from pathlib import Path

# Keep a reference so we can clean up
_FILE_HANDLER: FileHandler | None = None
_CONSOLE_HANDLER: StreamHandler | None = None
_LOGFILE: Path | None = None

## pytest_plugins = ["pestifer.util.pytest_plugin"]


@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    """Causes each test to run in the directory in which the module is found **or** a subdirectory with the same base name as the module

    For example, say you have <packagename>/tests/unit/test_foo.py.  If the directory
    <packagename>/tests/unit/test_foo/ exists, the test will run in that directory; if not, it will run in <packagename>/tests/unit/

    :param request: pytest request
    :type request: pytest request
    :param monkeypatch: pytest monkeypatch
    :type monkeypatch: pytest.Monkeypatch
    """
    module_bn=request.fspath.basename
    module_name,_=os.path.splitext(module_bn)
    test_dir=request.fspath.dirname
    subdir=os.path.join(test_dir,module_name)
    if os.path.isdir(subdir):
        monkeypatch.chdir(subdir)
    else:
        monkeypatch.chdir(test_dir)

def pytest_addoption(parser):
    """
    Adds the ``--runslow`` option to the ``pytest`` command to allow for 
    tests marked as "slow" to run
    """
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )

def pytest_configure(config):
    """
    Adds a file handler and console handler to pytest's logging
    """
    global _FILE_HANDLER, _LOGFILE, _CONSOLE_HANDLER

    # repo root; adjust if your conftest.py lives elsewhere
    root = Path(config.rootpath)
    logdir = root / "tests" / "logs"
    logdir.mkdir(parents=True, exist_ok=True)

    # xdist-friendly: one file per worker
    worker = os.getpid()
    _LOGFILE = logdir / f"tests-{worker}.log"

    fh = FileHandler(_LOGFILE, mode="w", encoding="utf-8", delay=False)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s",
                              "%Y-%m-%d %H:%M:%S"))
    fh._pytest_added = True  # tag for cleanup
    _FILE_HANDLER = fh

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)   # adjust to taste
    root_logger.addHandler(fh)

    ch = StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s",
                              "%Y-%m-%d %H:%M:%S"))
    ch._pytest_added = True  # tag for cleanup
    _CONSOLE_HANDLER = ch
    root_logger.addHandler(ch)

    # Optional: announce path at start
    config._logfile_path = str(_LOGFILE)

    config.addinivalue_line("markers", "slow: mark test as slow to run")

def pytest_unconfigure(config):
    """Remove and close our handler cleanly."""
    global _FILE_HANDLER, _CONSOLE_HANDLER
    if _FILE_HANDLER:
        root_logger = logging.getLogger()
        try:
            root_logger.removeHandler(_FILE_HANDLER)
        finally:
            try:
                _FILE_HANDLER.close()
            except Exception:
                pass
        _FILE_HANDLER = None
    if _CONSOLE_HANDLER:
        root_logger = logging.getLogger()
        root_logger.removeHandler(_CONSOLE_HANDLER)
        _CONSOLE_HANDLER = None

def pytest_report_header(config):
    p = getattr(config, "_logfile_path", None)
    return f"Logging to: {p}" if p else None

def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)

