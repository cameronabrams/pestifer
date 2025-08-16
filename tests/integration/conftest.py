from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import argparse
import importlib.resources as ir
import re
import pytest

from pestifer.core.examplemanager import ExampleManager
from pestifer.core.example import Example

PKG_NAME = "pestifer"                   # your package import name
CASES_DIR = ("resources", "examples")  # default discovery root: pkg/testdata/cases/


@dataclass(frozen=True)
class Case:
    index: int
    folder: str

def _example_foldername(example_id: int) -> Path:
    """
    Convert an example index to a folder name like 'ex01', 'ex02', etc.
    """
    return Path(Example.folder_name_format.format(example_id=example_id))


def _default_cases_root() -> Path:
    base = ir.files(PKG_NAME)
    for part in CASES_DIR:
        base = base / part
    return Path(str(base)).resolve()

example_manager = ExampleManager(examples_path=_default_cases_root())
def _example_foldername_validate(foldername: Path | str) -> bool:
    """
    Validate an example folder name like 'ex01', 'ex02', etc.
    """
    testname = foldername if isinstance(foldername, str) else foldername.name
    return example_manager.root_folder_format_validator.fullmatch(testname)

def _example_foldername_extract(foldername: Path | str) -> dict:
    testname = foldername if isinstance(foldername, str) else foldername.name
    return example_manager.root_folder_format_validator.extract(testname)

def _discover_example_folders(root: Path) -> list[Path]:
    if not root.exists():
        return []
    subfolders = []
    # look for subdirectories named with an index pattern
    for p in root.iterdir():
        if p.is_dir():
            if _example_foldername_validate(p):
                subfolders.append(p)
    return sorted(subfolders)

def pytest_addoption(parser: pytest.Parser) -> None:
    g = parser.getgroup("integration-examples")
    g.addoption(
        "--examples",
        action="store",
        default=None,
        help="Comma/ran­ge list of example indices to run, e.g. '1,2,5-7'. "
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
    if "case" not in metafunc.fixturenames:
        return
    root = _cases_root(metafunc.config)
    explicit = metafunc.config.getoption("--examples")
    if explicit:
        indices = _expand_examples(explicit)
        example_subfolders = [_example_foldername(i) for i in indices if _example_foldername(i).exists()]
    else:
        example_subfolders = _discover_example_folders(root)
        indices = [int(_example_foldername_extract(p.name).get('example_id')) for p in example_subfolders]
    cases = [Case(index=i, folder=_example_foldername(i)) for i in indices]
    metafunc.parametrize("case", cases, folders=[c.folder for c in cases])

@pytest.fixture
def per_case_dir(case: Case, request: pytest.FixtureRequest) -> Path:
    """
    Make a subdirectory under the test module's directory (works with your change_test_dir).
    Example: pkg/tests/integration/test_cli_cases/ex01/
    """
    test_module_dir = Path(request.node.fspath).parent.resolve()
    d = test_module_dir / f"ex{case.index:02d}"
    d.mkdir(parents=True, exist_ok=True)
    return d

@pytest.fixture
def make_namespace():
    """
    Build the argparse.Namespace expected by run_example, mirroring your current runner.
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
