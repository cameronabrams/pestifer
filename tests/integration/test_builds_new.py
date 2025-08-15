import os
from pathlib import Path
import logging
import pytest

from pestifer.core.pestifer import run_example  # adjust if your callable lives elsewhere

@pytest.mark.slow
def test_example_case(case, per_case_dir, make_namespace, caplog):
    """
    Mirrors your do_it(exnumber) logic but uses per-case subdirs under this test module.
    """
    caplog.set_level(logging.DEBUG)

    # Clean the per-case dir like your runner does
    for p in per_case_dir.iterdir():
        if p.is_dir():
            # don't nuke the folder itself; just clear contents
            for q in p.rglob("*"):
                try:
                    q.unlink()
                except IsADirectoryError:
                    pass
        else:
            p.unlink()

    # Chdir into the per-case directory
    os.chdir(per_case_dir)

    # Build args Namespace and run
    args = make_namespace(case.index, per_case_dir)
    controller = run_example(args)

    # Validate produced tarball
    prod_basename = controller.tasks[-1].specs["package"]["basename"]
    prod_tgz = Path(f"{prod_basename}.tar.gz")
    assert prod_tgz.exists(), f"Missing expected archive: {prod_tgz}"
