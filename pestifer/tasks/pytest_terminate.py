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
        """
        This method will gather all file artifacts labeled as "pytestable", delete other file artifacts,
        and then either 
        1. save the pytestables as gold
        2. compare the pytestables to the gold standard
        """
        assert running_under_pytest()
        # Determine all testable file artifacts
        artifact_collection = self.pipeline.get_all_file_artifacts()
        testable_file_artifacts = FileArtifactList(list(filter(lambda a: a.pytestable, artifact_collection)))
        # delete any non-testable file artifacts
        for file_artifact in artifact_collection:
            if not file_artifact in testable_file_artifacts:
                logger.debug(f'PytestTerminateTask deleting non-testable file artifact {file_artifact.name}')
                file_artifact.path.unlink(missing_ok=True)

        all_files = os.listdir('.')
        logger.debug('All files:')
        for f in all_files:
            if not f in artifact_collection:
                logger.debug(f'  Not an artifact -> {f}')

        logger.debug('Testable artifacts:')
        for f in testable_file_artifacts:
            logger.debug(f'  {f.name}')

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
            self.register(results, key='test_results')
            logger.debug(f'Registered all test results at "test_results"')
        return 0