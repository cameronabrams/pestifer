# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of a custom TerminateTask that identifies and exposes outputs for testing purposes.
"""

import logging

from ..tasks.basetask import BaseTask
from ..core.artifacts import *
from ..util.util import running_under_pytest

from ..util._goldenmode import generate_gold, report_example_id

logger = logging.getLogger(__name__)

class PytestTerminateTask(BaseTask):

    _yaml_header = 'pytest_terminate'

    def do(self) -> int:
        assert running_under_pytest()
        # Determine all testable file artifacts
        task_results = self.pipeline.get_artifact_collection_as_lists()
        testable_file_artifacts = FileArtifactList(list(filter(lambda a: a.pytestable, task_results['files']))).unique_paths()
        # delete any non-testable file artifacts
        for file_artifact in task_results['files']:
            if not file_artifact in testable_file_artifacts:
                file_artifact.path.unlink(missing_ok=True)
        for file_artifact_list in task_results['filelists']:
            for file_artifact in file_artifact_list:
                if not file_artifact in testable_file_artifacts:
                    if file_artifact.pytestable:
                        testable_file_artifacts.append(file_artifact)
                    else:
                        file_artifact.path.unlink(missing_ok=True)

        # if we are generating new standards, simply copy the remaining
        # file artifacts to the standards directory for this test.  It
        # needs to know the example_id for the test, so that will be in
        # the task specs.
        example_id = report_example_id()
        example_manager = self.resource_manager.example_manager
        logger.debug(f'PytestTerminateTask using example_id={example_id}')
        if generate_gold():
            logger.debug('PytestTerminateTask generating gold')
            if example_id is None:
                raise ValueError("example_id must be specified in the specs for PytestTerminateTask")
            example_manager.update_example(
                example_id, 
                outputs=[x.name for x in testable_file_artifacts],
                skip_sphinx=True
            )
        else:
            logger.debug('PytestTerminateTask comparing to gold')
            example = example_manager.examples.get_example_by_example_id(example_id)
            golden_folder = example_manager.path / example.outputspath
            results = {}
            for file_artifact in testable_file_artifacts:
                results[file_artifact.name] = "pass" if file_artifact.compare(golden_folder / file_artifact.name) else "fail"
            self.register(DataArtifact(key='test_results', data=results))
            logger.debug(f'Registered all test results at "test_results"')
        return 0