"""
Pytest configuration for integration tests
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import argparse
import importlib.resources as ir
import os
import pytest

from pestifer.core.examplemanager import ExampleManager
from pestifer.core.example import Example

@dataclass(frozen=True)
class Case:
    example_id: int

PKG_NAME = "pestifer"                   # your package import name
CASES_DIR = ("resources", "examples")  # default discovery root: pkg/testdata/cases/

# pytest_plugins = ["pestifer.util.pytest_plugin"]

def _default_cases_root() -> Path:
    base = ir.files(PKG_NAME)
    for part in CASES_DIR:
        base = base / part
    return Path(str(base)).resolve()

def pytest_addoption(parser: pytest.Parser) -> None:
    """
    Adds the ``--examples`` and ``--cases-root`` options to pytest CLI.
    """
    g = parser.getgroup("integration-examples")
    g.addoption(
        "--examples",
        action="store",
        default=None,
        help="Comma/range list of example indices to run, e.g. '1,2,5-7'. "
             "Overrides discovery."
    )
    g.addoption(
        "--cases-root",
        action="store",
        default=None,
        help="Root directory to discover examples (default: pestifer/resources/examples)"
    )

def _expand_examples(spec: str) -> list[int]:
    # e.g., "1,2,5-7"
    out: list[int] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if "-" in chunk:
            a, b = chunk.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        elif chunk:
            out.append(int(chunk))
    return sorted(set(out))

def _cases_root(config: pytest.Config) -> Path:
    cr = config.getoption("--cases-root")
    return Path(cr).resolve() if cr else _default_cases_root()

def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """
    Generates the test cases based on the available examples.
    """
    if "case" not in metafunc.fixturenames:
        return
    _example_manager = ExampleManager(examples_path=_cases_root(metafunc.config))
    explicit = metafunc.config.getoption("--examples")
    if explicit:
        indices = _expand_examples(explicit)
    else:
        indices = [ex.example_id for ex in _example_manager.examples.data]
    cases = [Case(example_id=i) for i in indices]
    metafunc.parametrize("case", cases)

@pytest.fixture
def per_case_dir(case: Case, request: pytest.FixtureRequest) -> Path:
    """
    Makes a subdirectory under the test module's directory (works with change_test_dir in the top-level conf.py).
    
    Example: pkg/tests/integration/test_cli_cases/ex01/
    """
    test_module_dir = Path(request.node.fspath).parent.resolve()
    d = test_module_dir / Path('__test_builds') / Example.folder_name_format.format(example_id=case.example_id)
    os.environ["PESTIFER_PYTEST_GENERATE_GOLD_EXAMPLE_ID"] = str(case.example_id)
    d.mkdir(parents=True, exist_ok=True)
    return d

@pytest.fixture
def make_namespace():
    """
    Build the argparse.Namespace expected by run_example.
    """
    def _make(example_id: int, outdir: Path) -> argparse.Namespace:
        return argparse.Namespace(
            example_id=example_id,
            config=None,
            output_dir=str(outdir),
            log_file=f"diagnostics.log",
            gpu=False,
            complete_config=False,
            log_level="debug",
        )
    return _make
