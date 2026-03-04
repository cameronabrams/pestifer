# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`ResourceManager` class for managing access to :mod:`pestifer.resources`.
"""

import logging
import os

from pathlib import Path
from typing import Callable

from .examplemanager import ExampleManager
from .labels import Labels

from .. import resources

from ..charmmff.charmmffcontent import CHARMMFFContent, charmmff_version_key
from ..util.gitutil import get_git_origin_url

logger = logging.getLogger(__name__)

class ResourceManager:
    """
    A class for managing pestifer's built-in resources.
    This class provides access to various resources such as CHARMM force fields, example input files, and TcL scripts.
    The resources are organized into directories, and the class provides methods to access these resources.
    """

    _base_resources = ('charmmff', 'examples', 'tcl')
    """
    A list of base resources that are managed by the ResourceManager.
    These resources include CHARMM force fields, example input files, and TcL scripts.
    """

    def __init__(self, charmmff_config: dict = {}):
        self.resources_path = os.path.dirname(resources.__file__)
        self.package_path = Path(self.resources_path).parent.parent

        resources_subfolders = os.listdir(self.resources_path)
        self.resource_dirs = [x for x in resources_subfolders if os.path.isdir(os.path.join(self.resources_path, x)) and x in ResourceManager._base_resources]
        assert all([x in [os.path.basename(_) for _ in self.resource_dirs] for x in ResourceManager._base_resources]), f'some resources seem to be missing; found {self.resource_dirs}, expected {ResourceManager._base_resources}'

        self.resource_path = {}
        for r in self._base_resources:
            self.resource_path[r] = os.path.join(self.resources_path, r)
            if not os.path.isdir(self.resource_path[r]):
                raise FileNotFoundError(f'Resource {r} not found at {self.resource_path[r]} -- your installation is likely incomplete')

        self._charmmff_config = charmmff_config
        self._charmmff_content = None
        self._example_manager = None
        self.labels = Labels

    def charmmff_version_dirs(self) -> list[Path]:
        """Return all available charmmff version directories, sorted by release date (oldest first)."""
        import re
        import datetime
        _month_num = {
            'jan': 1, 'feb': 2, 'mar': 3, 'apr': 4, 'may': 5, 'jun': 6,
            'jul': 7, 'aug': 8, 'sep': 9, 'oct': 10, 'nov': 11, 'dec': 12,
        }
        def _key_to_date(d: Path) -> datetime.date:
            mon = _month_num.get(d.name[:3], 0)
            yr = 2000 + int(d.name[3:])
            return datetime.date(yr, mon, 1)
        charmmff_dir = Path(self.resource_path['charmmff'])
        return sorted(
            [d for d in charmmff_dir.iterdir() if d.is_dir() and re.match(r'^[a-z]{3}\d{2}$', d.name)],
            key=_key_to_date
        )

    def charmmff_version_path(self, version_str: str = '') -> Path:
        """Return the version-specific charmmff directory path.

        If *version_str* is given (e.g. ``'July2024'``), the corresponding
        subdirectory is returned.  Otherwise the most-recently-modified
        version directory is used as the default.
        """
        if version_str:
            key = charmmff_version_key(version_str)
            path = Path(self.resource_path['charmmff']) / key
            if not path.is_dir():
                raise FileNotFoundError(f'CHARMMFF version directory for {version_str!r} not found at {path}')
            return path
        version_dirs = self.charmmff_version_dirs()
        if not version_dirs:
            raise FileNotFoundError(f'No CHARMMFF version directories found in {self.resource_path["charmmff"]}')
        return version_dirs[-1]  # newest by mtime

    @property
    def charmmff_content(self):
        if self._charmmff_content is None:
            charmmff_path = self.charmmff_version_path(self._charmmff_config.get('release', ''))
            logger.info(f'Using CHARMMFF version: {charmmff_path.name}')
            self._charmmff_content = CHARMMFFContent(
                charmmff_path,
                user_custom_directories=self._charmmff_config.get('user_custom', {}).get('searchpath', []),
                user_pdbrepository_paths=self._charmmff_config.get('pdbrepository', []),
            )
        return self._charmmff_content

    @charmmff_content.setter
    def charmmff_content(self, value):
        self._charmmff_content = value

    @property
    def example_manager(self):
        if self._example_manager is None:
            is_source_package_with_git = os.path.isdir(os.path.join(self.package_path, '.git'))
            if is_source_package_with_git:
                remote = get_git_origin_url()
                if remote:
                    logger.debug(f'Installed as source from {remote} into {self.package_path}')
                else:
                    logger.debug(f'Installed as source into {self.package_path} but the remote origin URL could not be determined.')
            docs_source_path = None
            if is_source_package_with_git:
                docs_source_path = os.path.join(self.package_path, 'docs', 'source')
                if os.path.isdir(docs_source_path):
                    logger.debug(f'Docs path {docs_source_path} exists; using it for example documentation')
                else:
                    logger.debug(f'Error: This is a source package but docs path {docs_source_path} does not exist')
                    docs_source_path = None
            self._example_manager = ExampleManager(examples_path=self.resource_path['examples'], docs_source_path=docs_source_path)
        return self._example_manager

    def __str__(self):
        cp = os.path.commonpath(list(self.resource_path.values()))
        retstr = f'Resources are found under\n    {cp}\n'
        for r, p in self.resource_path.items():
            retstr += f'        {p.replace(cp + os.sep, "") + os.sep}\n'
        return retstr

    def show(self, out_stream: Callable = print, components: dict = {}, fullnames=False, missing_fullnames={}):
        """
        Show a summary of the specified components.
        
        Parameters
        ----------
        out_stream : callable, optional
            A callable that takes a string and outputs it. Default is `print`.
        components : dict, optional
            A dictionary specifying the components to show. The keys are the component names
            and the values are the specifications for each component.
        fullnames : bool, optional
            If True, show full names for the components. Default is False.
        missing_fullnames : dict, optional
            A dictionary of missing full names for the components. Default is an empty dictionary.
        """
        for c, spec in components.items():
            if not c in self._base_resources:
                logger.warning(f'{c} is not a base resource; expected one of {", ".join(self._base_resources)}')
            path = self.get_resource_path(c)
            if c == 'examples':
                out_stream(f'\nExamples:\n\n{self.example_manager.report_examples(header=True)}')
            elif c == 'charmmff':
                if 'tarball' in spec:
                    out_stream(f'{self.charmmff_content.tarfilename}')
                if 'pdb' in spec:
                    self.charmmff_content.provision_pdbrepository()
                    self.charmmff_content.pdbrepository.show(out_stream, fullnames=fullnames, missing_fullnames=missing_fullnames)
                if 'custom' in spec:
                    path = self.get_charmmff_customdir()
                    with open(os.path.join(path, '00PESTIFER-README.txt'), 'r') as f:
                        msg = f.read()
                    out_stream(msg)
            elif c == 'tcl':
                path = self.get_tcldir()
                with open(os.path.join(path, '00PESTIFER-README.txt'), 'r') as f:
                    msg = f.read()
                out_stream(msg)

    def get_resource_path(self,r):
        """
        Get the path to a specific resource.

        Parameters
        ----------
        r : str
            The name of the resource to get the path for. This should be one of the base resources defined in ``ResourceManager._base_resources``.
        """
        return self.resource_path.get(r,None)
            
    def get_charmmff_customdir(self):
        """
        Get the path to the custom CHARMM force field directory.
        This directory is used for storing custom CHARMM force field files that are not part of the standard distribution.

        Returns
        -------
        str
            The path to the custom CHARMM force field directory.
        """
        return str(self.charmmff_content.charmmff_path / 'custom')

    def get_tcldir(self):
        """
        Get the path to the TcL scripts and packages directory.
        
        Returns
        -------
        str
            The path to the TcL scripts and packages directory.
        """
        return self.resource_path['tcl']
    
    def update_atomselect_macros(self):
        """
        Update the atomselect macros in the macros.tcl file based on the content of labels.py.
        This method updates atomselect macros for residues that aren't already part of
        atomselect macros in a base VMD installation.
        """
        p=os.path.join(self.resource_path['tcl'],'macros.tcl')
        with open(p,'w') as f:
            f.write('## Updates atomselect macros to account for residues not in the macros in a base VMD installation.\n')
            f.write('## This file is automatically generated by ``pestifer.core.resourcemanager.update_atomselect_macros().``\n')
            f.write('## Do not edit this file directly; changes will be overwritten.\n')
            f.write('## Instead, edit :mod:`pestifer.core.labels` and run ``pestifer modify-package --update-atomselect-macros.``\n')
            f.write('## The Tcl proc ``update_atomselect_macro`` is defined in the PestiferUtil Tcl package.\n\n')
            self.labels.update_atomselect_macros(f)

    def get_tcl_pkgdir(self):
        """
        Get the path to the TcL package directory.
        This directory contains TcL packages that are used by pestifer for various tasks, such as loading additional functionality or libraries.

        Returns
        -------
        str
            The path to the TcL package directory.
        """
        return os.path.join(self.resource_path['tcl'], 'pkg')
    
    def get_tcl_scriptsdir(self):
        """
        Get the path to the TcL scripts directory.
        This directory contains TcL scripts that are used by pestifer for various tasks, such as setting up simulations or processing data.

        Returns
        -------
        str
            The path to the TcL scripts directory.
        """
        return os.path.join(self.resource_path['tcl'], 'scripts')
    

    def append_example(self, scriptname: str, example_id: int = 0, author_name='', author_email='', title='', db_id='', auxiliary_inputs=[], outputs=[]):
        """
        Add an example to the pestifer package.  Minimally, the name of the example's YAML input file must be specified.

        Parameters
        ----------
        scriptname : str
            The name of the example's YAML input file (with or without the .yaml extension).
        example_id : int
            The unique identifier for the example. If not provided, a new ID will be assigned.
        author_name : str
            The name of the author.
        author_email : str
            The email address of the author.
        title : str
            The title of the example.
        db_id : str
            The database ID of the example.
        auxiliary_inputs : list
            A list of auxiliary input files for the example.
        outputs : list
            A list of output files for the example.
        """
        example_id = len(self.example_manager.examples) + 1 if example_id == 0 else example_id
        self.example_manager.append_example(example_id, scriptname, title=title, author_name=author_name, author_email=author_email, auxiliary_inputs=auxiliary_inputs, outputs=outputs, db_id=db_id)

    def update_example(self, example_id: int, shortname: str = '', title: str = '', db_id: str = '', author_name: str = '', author_email: str = '', auxiliary_inputs: list = [], outputs: list = []):
        """
        Update an existing example in the pestifer examples directory.

        Parameters
        ----------
        example_id : int
            The unique identifier for the example to update.
        shortname : str
            The short name of the example.
        title : str
            The title of the example.
        db_id : str
            The database ID of the example.
        author_name : str
            The name of the author.
        author_email : str
            The email address of the author.
        auxiliary_inputs : list
            A list of auxiliary input files for the example.
        outputs : list
            A list of output files for the example.
        """
        self.example_manager.update_example(example_id, shortname=shortname, title=title, db_id=db_id, author_name=author_name, author_email=author_email, auxiliary_inputs=auxiliary_inputs, outputs=outputs)

    def delete_example(self, example_id: int):
        """
        Delete an example from the pestifer examples directory.

        Parameters
        ----------
        example_index : int
            The index of the example to delete. This should be a positive integer indicating the position in the examples list (1-based).
        """
        self.example_manager.delete_example(example_id)

    def rename_example(self, example_id: int, new_name: str):
        """
        Rename an existing example in the pestifer examples directory.

        Parameters
        ----------
        example_id : int
            The unique identifier for the example to rename.
        new_name : str
            The new name for the example. This should be a valid filename without any path components.
        """
        self.example_manager.rename_example(example_id, new_name)

    def set_example_author(self, example_id: int, author_name: str, author_email: str):
        """
        Set the author information for an example in the pestifer examples directory.

        Parameters
        ----------
        example_id : int
            The unique identifier for the example to set the author for.
        author_name : str
            The name of the author.
        author_email : str
            The email address of the author.
        """
        self.example_manager.set_example_author(example_id, author_name, author_email)