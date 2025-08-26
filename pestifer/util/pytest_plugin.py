# Author: ChatGPT 5
import os

def pytest_addoption(parser):
    group = parser.getgroup("integration-examples")
    group.addoption(
        "--generate-gold",
        action="store_true",
        help="Write golden files instead of comparing to them."
    )

def pytest_configure(config):
    from ._goldenmode import set_generate_gold
    flag = config.getoption("--generate-gold")
    set_generate_gold(flag)

    if flag:
        os.environ["PESTIFER_PYTEST_GENERATE_GOLD"] = "1"
    else:
        os.environ.pop("PESTIFER_PYTEST_GENERATE_GOLD", None)
