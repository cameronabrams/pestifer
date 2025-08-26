# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
This module defines the ``test_example_case`` function for testing example cases and generating their gold-standard outputs for the tests.
"""

import os
import logging
import pytest
import shutil

from pestifer.core.artifacts import DataArtifact
from pestifer.subcommands import RunExampleSubcommand
from pestifer.util._goldenmode import generate_gold, set_example_id

from ..conftest import _LOGFILE

@pytest.mark.slow
def test_example_case(case, per_case_dir, make_namespace, caplog):
    """
    Runs an example in a per-case subdir.
    """
    caplog.set_level(logging.DEBUG)
    global _LOGFILE
    # Clean the per-case dir
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
    args = make_namespace(case.example_id, per_case_dir)
    task = RunExampleSubcommand()
    controller = task.func(args)
    shutil.copy(_LOGFILE, per_case_dir / "pestifer_diagnostics.log")
    if generate_gold():
        set_example_id(case.example_id)
    else:
        results_artifact: DataArtifact = controller.pipeline.get_current_artifact('test_results')
        assert results_artifact is not None, "Expected test results artifact not found"
        passlist = []
        faillist = []
        for file_name, result in results_artifact.data.items():
            if result == "pass":
                passlist.append(f"Test passed for {file_name}: {result}")
            else:
                faillist.append(f"Test failed for {file_name}: {result}")
        assert len(faillist) == 0, f"Some tests failed:\n" + "\n".join(faillist)