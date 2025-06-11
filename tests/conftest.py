#Author: Cameron F. Abrams, <cfa22@drexel.edu>

import pytest
import os

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
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)