from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import argparse
import importlib.resources as ir
import re
import pytest

PKG_NAME = "pestifer"                   # your package import name
CASES_DIR = ("pestifer", "resources", "examples")  # default discovery root: pkg/testdata/cases/

@dataclass(frozen=True)
class Case:
    index: int
    id: str

def _default_cases_root() -> Path:
    base = ir.files(PKG_NAME)
    for part in CASES_DIR:
        base = base / part
    return Path(str(base)).resolve()

def _parse_index(name: str) -> int | None:
    """
    Accept names like 'ex01', 'ex1', 'case_01', 'case-12', '01', '1'
    """
    m = re.search(r'(?:ex|case[_-]?)?(\d+)$', name)
    return int(m.group(1)) if m else None

def _discover_indices(root: Path) -> list[int]:
    if not root.exists():
        return []
    indices: set[int] = set()
    # 1) look for subdirectories named with an index pattern
    for p in root.iterdir():
        if p.is_dir():
            idx = _parse_index(p.name)
            if idx is not None:
                indices.add(idx)
    # 2) also allow files like ex01.in or case_02.yaml as hints
    for p in root.rglob("*"):
        if p.is_file():
            idx = _parse_index(p.stem)
            if idx is not None:
                indices.add(idx)
    return sorted(indices)

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
        help="Root directory to discover examples (default: pestifer/pestifer/resources/examples)"
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
    else:
        indices = _discover_indices(root)
    cases = [Case(index=i, id=f"ex{ i:02d }") for i in indices]
    metafunc.parametrize("case", cases, ids=[c.id for c in cases])

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
    def _make(exnumber: int, outdir: Path) -> argparse.Namespace:
        nstr = f"{exnumber:02d}"
        return argparse.Namespace(
            index=exnumber,
            config=None,
            output_dir=str(outdir),
            log_file=f"diagnostics-{nstr}.log",
            gpu=False,
            complete_config=False,
            log_level="debug",
        )
    return _make
