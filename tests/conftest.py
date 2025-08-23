#Author: Cameron F. Abrams, <cfa22@drexel.edu>
from __future__ import annotations

import pytest
import os
# conftest.py
import logging

from logging import FileHandler, Formatter
from pathlib import Path

# Keep a reference so we can clean up
_FILE_HANDLER: FileHandler | None = None
_LOGFILE: Path | None = None

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
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )

def pytest_configure(config):
    """Runs very early; attach our FileHandler to the root logger."""
    global _FILE_HANDLER, _LOGFILE

    # repo root; adjust if your conftest.py lives elsewhere
    root = Path(config.rootpath)
    logdir = root / "tests" / "logs"
    logdir.mkdir(parents=True, exist_ok=True)

    # xdist-friendly: one file per worker
    worker = os.getpid()
    _LOGFILE = logdir / f"tests-{worker}.log"

    fh = FileHandler(_LOGFILE, mode="w", encoding="utf-8", delay=False)
    fh.setFormatter(Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s",
                              "%Y-%m-%d %H:%M:%S"))
    fh._pytest_added = True  # tag for cleanup
    _FILE_HANDLER = fh

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)   # adjust to taste
    root_logger.addHandler(fh)

    # Optional: announce path at start
    config._logfile_path = str(_LOGFILE)

    config.addinivalue_line("markers", "slow: mark test as slow to run")

def pytest_unconfigure(config):
    """Remove and close our handler cleanly."""
    global _FILE_HANDLER
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

