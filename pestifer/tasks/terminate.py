# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`TerminateTask` class for terminating a pestifer build.
This task is a descendant of the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to prepare the system for termination.
It handles the copying of state files, writing chain maps, and packaging the system for NAMD runs.
The task also manages the state of the simulation, including the base molecule and various file extensions such as PSF, PDB, COOR, XSC, and VEL.
The task is designed to be used in a workflow where the simulation needs to be gracefully terminated and packaged for further analysis or continuation.
It ensures that all necessary files are collected and organized, making it easy to resume or analyze the simulation later.
"""
import logging
import shutil
import yaml
import os

from .mdtask import MDTask
from ..core.artifacts import *
from ..molecule.molecule import Molecule
from ..util.colors import PestiferColors
from ..util.stringthings import my_logger

logger = logging.getLogger(__name__)

class TerminateTask(MDTask):
    """
    TerminateTask class for terminating a pestifer build.
    This class inherits from the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to prepare the system for termination.
    It handles the copying of state files, writing chain maps, and packaging the system for NAMD runs.
    """
    _yaml_header = 'terminate'
    """
    YAML header for the TerminateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare TerminateTask objects in YAML task lists.
    It is typically used at the end of a build workflow to finalize the state of the simulation and prepare for termination.
    """

    def do(self) -> int:
        self.next_basename()
        if 'chainmapfile' in self.specs:
            self.write_chainmaps()
        self.result = self.test_standard()
        self.result += self.make_package() 
        self.result += self.cleanup()
        return self.result

    def write_chainmaps(self):
        """
        Write the chain maps to a YAML file.
        This method retrieves the base molecule from the state variables, gets the chain maps, and writes them to a specified YAML file.
        The chain maps are used to map chains in the molecular structure, which is useful for understanding the topology of the system.
        """
        bm: Molecule = self.get_current_artifact_data('base_molecule')
        if bm:
            maps = bm.get_chainmaps()
            with open(self.specs['chainmapfile'], 'w') as f:
                yaml.dump(maps, f)
            self.register(self.specs['chainmapfile'].replace('.yaml', ''), key='chainmapfile', artifact_type=YAMLFileArtifact)

    def make_package(self):
        """
        Create a package for a production NAMD run starting from the end of the build.
        """
        package_specs = self.specs.get('package', {})
        if not package_specs:
            return 0
        md_specs = package_specs.get('md', {})
        state_dir = package_specs.get('state_dir', '.')
        if not package_specs:
            logger.debug('No package specifications provided; packaging will not be performed.')
            return 0
        TarballContents = FileArtifactList()
        self.basename = self.specs.get('basename', 'my_system')
        state: StateArtifacts = self.get_current_artifact('state')
        new_state = {}
        for ext in ['psf', 'pdb', 'coor', 'xsc', 'vel']:
            fa: FileArtifact = getattr(state, ext, None)
            if fa and fa.exists():
                logger.debug(f'Including {fa.name} in package as {self.basename}.{ext}')
                shutil.copy(fa.name, self.basename + '.' + ext)
                new_fa = fa.copy(data=self.basename + '.' + ext)
                TarballContents.append(new_fa)
                new_state[ext] = new_fa
        self.register(new_state, key='state', artifact_type=StateArtifacts)
        result = 0
        if md_specs:
            logger.debug(f'Packaging for namd using basename {self.basename}')
            save_specs = self.specs
            self.specs = md_specs
            self.specs['basename'] = package_specs.get('basename', self.basename)
            result = self.namdrun(script_only=True)
            self.specs = save_specs
            TarballContents.append(self.get_current_artifact('namd'))
            constraints = self.specs.get('constraints', {})
            if constraints:
                self.make_constraint_pdb(constraints, statekey='consref')
                TarballContents.append(self.get_current_artifact('consref'))
        else:
            logger.debug(f'No NAMD configuration is included in the package.')
        TarballContents.extend(self.get_current_artifact_data('charmmff_parfiles'))
        TarballContents.extend(self.get_current_artifact_data('charmmff_streamfiles'))
        TarballContents.make_tarball(self.basename, arcname_prefix=state_dir, unique=True, remove=True)
        return result

    def cleanup(self):

        if not self.specs.get('cleanup', True):
            logger.debug('Cleanup disabled; skipping cleanup step.')
            return 0
        archive_dir = self.specs.get('archive_dir', 'archive')  

        file_artifacts: FileArtifactList = self.pipeline.get_all_file_artifacts()
        file_artifacts.sort(key=lambda x: x.name)
        logger.debug(f'{len(file_artifacts)} file artifacts to be included in archive:')
        my_logger([fa.name for fa in file_artifacts.data], logger.debug, depth=1)
        file_artifact_names = [fa.name for fa in file_artifacts.data]
        non_artifact_files = []
        cwd_files = os.listdir('.')
        for f in cwd_files:
            if f not in file_artifact_names:
                non_artifact_files.append(f)
        if len(non_artifact_files) > 0:
            logger.debug(f'Non-artifact files in current working directory:')
            my_logger(non_artifact_files, logger.debug, depth=1)

        file_artifacts.make_tarball('artifacts', remove=True, arcname_prefix=archive_dir, unique=True)
        return 0
    
    def test_standard(self):
        standard_specs = self.specs.get('test_standards', {})
        if not standard_specs:
            logger.debug('No test_standards specifications provided; skipping test_standards step.')
            return 0
        artifact_file_collection = self.pipeline.get_all_file_artifacts()
        testable_file_artifacts = FileArtifactList(list(filter(lambda a: a.pytestable, artifact_file_collection)))
        if len(testable_file_artifacts) == 0:
            logger.debug('No pytestable artifacts found; skipping test_standard step.')
            return 0
        standard_path = standard_specs.get('standard_path', None)
        if standard_path is None:
            logger.debug('No standard_path provided; skipping test_standard step.')
            return 0
        standard_path = Path(standard_path).resolve()
        if not standard_path.exists():
            # create it
            standard_path.mkdir(parents=True, exist_ok=True)
        elif not standard_path.is_dir():
            logger.debug(f'standard_path {standard_path} is not a directory; skipping test_standard step.')
            return 0
        # standards path cannot be CWD
        if standard_path.absolute() == Path.cwd().absolute():
            logger.debug('standard_path cannot be the current working directory; skipping test_standard step.')
            return 0
        logger.debug(f'Using standard_path: {standard_path.name}')

        all_files = [x for x in standard_path.iterdir() if x.is_file()]
        if len(all_files) == 0 and standard_path.is_dir():
            logger.debug(f'standard_path {standard_path} is empty; populating it with current testable artifacts.')
            for f in testable_file_artifacts:
                logger.debug(f'  Populating standard with {f.name}')
                shutil.copy(f.name, standard_path / f.name)
            return 0
        elif not standard_path.is_dir():
            logger.debug(f'standard_path {standard_path} does not exist; creating it and populating with current testable artifacts.')
            standard_path.mkdir(parents=True, exist_ok=True)
            for f in testable_file_artifacts.data:
                logger.debug(f'  Populating standard with {f.name}')
                shutil.copy(f.name, standard_path / f.name)
            return 0
        else:
            logger.debug('Testable artifacts:')
            results = {}
            for f in testable_file_artifacts.data:
                logger.debug(f'  {f.name}')
                results[f.name] = "pass" if f.compare(standard_path / f.name) else "fail"
            self.register(results, key='test_results')
            logger.debug(f'Registered all test results at "test_results"')
            my_logger('Standard Test Results:', logger.info)
            maxnamelen = max([len(n) for n in results.keys()]) if len(results) > 0 else 0
            name_format = f'- {{:>{maxnamelen}s}}'
            num_items = len(results)
            colors = dict(_pass=PestiferColors['emerald'], _fail=PestiferColors['alizarin'])
            for item, (name, result) in enumerate(results.items()):
                color = colors.get(f'_{result}', PestiferColors['alizarin'])
                logger.info(f'  {name_format.format(name)}: {color}{result}{color.OFF} ({item+1}/{num_items})')
            logger.info('*'*70)
        return 0