# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This module defines classes that facilitate the handling of CHARMM force field content.
"""
import logging
import os
import re

from pathlib import Path

from .charmmfftop import CharmmMassDict, CharmmMassList, CharmmResiDict, CharmmResi
from .pdbrepository import PDBRepository, PDBInput
from ..core.labels import Labels
from ..util.cacheable_object import CacheableObject, TarBytesFS
from ..util.patch import apply_unified_diff
from ..util.util import countTime
from ..util.spinner_wrapper import with_spinner
from ..util.stringthings import my_logger

logger = logging.getLogger(__name__)

def parse_conditional_script(script_text):
    """
    Parse a conditional script text and return a dictionary with parsed lines and variables.

    This function processes a script that contains conditional statements and variable assignments.
    It supports 'set' commands to define variables and 'if' statements to conditionally execute blocks of code.

    This function was written by ChatGPT 4o.

    Parameters
    ----------
    script_text : str
        The script text to parse, which may contain 'set', 'if', and 'endif' statements.
        
    Returns
    -------
    dict
        A dictionary with two keys:

            - ``parsed``: a string containing the processed script text with comments and conditionals resolved.
            - ``vars``: a dictionary of variables defined in the script.
    """

    lines = script_text.strip().splitlines()
    vars = {}
    output = []
    i = 0
    stop_processing = False

    while i < len(lines):
        if stop_processing:
            break

        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue

        tokens = line.split()

        if tokens[0] == "set" and len(tokens) >= 3:
            var, value = tokens[1], tokens[2]
            if value.isdigit():
                vars[var] = int(value)
            else:
                vars[var] = value.strip('"')
            output.append(line)

        elif tokens[0] == "if" and tokens[2] == 'eq':
            var, op, value = tokens[1], tokens[2], tokens[3]
            if value.isdigit():
                value = int(value)
            else:
                value = value.strip('"')
            condition_result = False
            if op == "eq":
                condition_result = vars.get(var) == value
            block_lines = []
            contains_return = False
            i += 1
            while i < len(lines):
                inner_line = lines[i].strip()
                if inner_line == "endif":
                    break
                block_lines.append(inner_line)
                if inner_line == "return":
                    contains_return = True
                i += 1

            if condition_result:
                output.extend(block_lines)
                if contains_return:
                    stop_processing = True
            # move past the endif
        else:
            output.append(line)

        i += 1
    return dict(parsed="\n".join(output), vars=vars)

def extract_resi_pres_blocks(text: str, keywords: tuple[str, ...] = ('RESI', 'PRES')):
    """ 
    Extract blocks of text starting with RESI or PRES and ending before the next RESI, PRES, ATOMS, end, or EOF.

    This function uses a regular expression to find blocks of text that start with the specified keywords
    and continue until the next occurrence of one of the keywords or the end of the file.
    
    Parameters
    ----------
    text : str
        The input text from which to extract the blocks.
    keywords : tuple of str, optional
        The keywords that indicate the start of a block. 
        Defaults to (``RESI``, ``PRES``).

    Returns
    -------
    list of str
        A list of strings, each containing a block of text that starts with one 
        of the specified keywords.
    """
    
    keyword_pattern = '|'.join(re.escape(k) for k in keywords)
    # Include ATOMS and BONDS as additional block-end sentinels
    terminator_pattern = rf'^(?:{keyword_pattern}|ATOMS|BONDS|ANGLES|DIHEDRALS|end|END|read|READ)\b'
    # Match starting with RESI or PRES, ending before the next RESI, PRES, ATOMS, or EOF
    pattern = rf'(?i)^(({keyword_pattern})\s.*?)(?={terminator_pattern}|\Z)'
    matches = re.finditer(pattern, text, flags=re.DOTALL | re.MULTILINE)
    return [m.group(1).strip() for m in matches]

def extract_mass_lines(file_contents):
    """
    Extract lines containing mass information from the CHARMM force field files.
    This function scans the contents of a CHARMM force field file and returns lines that start with ``MASS``.

    Parameters
    ----------
    file_contents : str
        The contents of the CHARMM force field file to scan for mass information.
    """
    return [line for line in file_contents.splitlines() if line.strip().upper().startswith("MASS")]

class CHARMMFFResiTopCollection(CacheableObject):
    """
    A collection of CHARMM residue topology data.
    """

    @countTime
    def __init__(self, *args, **kwargs):
        is_custom = 'residues' in kwargs and len(kwargs['residues']) > 0
        if is_custom:
            self.build_custom(*args, **kwargs)
        else:
            super().__init__(*args, **kwargs)

    @with_spinner('Building residue topology collection from package resources...')
    def _build_from_resources(self, charmmff_path: str = '', **kwargs):
        """
        Build the collection from the specified resources.

        Parameters
        ----------
        path_or_tarball : str
            The path to the directory or tarball containing the CHARMM residue topology files.
        streamID_override : str
            An optional stream ID to override the default.

        Returns
        -------
        CHARMMFFResiTopCollection
            The constructed CHARMMFFResiTopCollection object.
        """
        local_charmmffcontent = CHARMMFFContent(charmmff_path)
        # this should read from cache

        logger.debug(f'Parsing all RESI and PRES objects...')
        local_charmmffcontent.find_resis_and_patches()
        self.residues = local_charmmffcontent.residues
        self.patches = local_charmmffcontent.patches
        logger.debug(f'CHARMMFFContentData initialized with {len(self.residues)} residues and {len(self.patches)} patches.')

class CHARMMFFContent(CacheableObject):
    """
    Holds all CHARMM force field content parsed for use within Pestifer.

    The CHARMM force field is downloadable directly from the
    MacKerell lab at the University of Maryland:
    https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
    Pestifer uses its own local copy of this tarball.    

    Parameters
    ----------
    charmmff_path : str, optional
        Path to the directory containing the CHARMM force field files. Default is current directory.
    tarfilename : str, optional
        Name of the tarball file containing the CHARMM force field files. Default is 'toppar_c36_jul24.tgz'.

    Attributes
    ----------
    charmmff_path : str
        Path to the directory containing the CHARMM force field files.
    tarfilename : str
        Name of the tarball file containing the CHARMM force field files.
    filenamemap : dict
        Maps file basenames to their full paths in the CHARMM force field content.
    """

    @countTime
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.deprovision()

    @with_spinner('No cache yet -- building all CHARMMFF content from package resources...')
    def _build_from_resources(self, charmmff_path: Path, **kwargs):
        """ Method to build the CHARMMFFContent object from resources, if the cache is stale. """
        self.filenamemap = {}
        if not charmmff_path.is_dir():
            raise NotADirectoryError(f'Expected a directory at {charmmff_path}, but it is not a directory')
        self.charmmff_path = charmmff_path.absolute()
        self.basename = self.charmmff_path.name
        self.parent_path = self.charmmff_path.parent
        self.charmm_elements = [x.name for x in list(self.charmmff_path.glob('*'))]
        tarfilename = kwargs.get('tarfilename', 'toppar_c36_jul24.tgz')
        skip_streams = kwargs.get('skip_streams', ['misc', 'cphmd'])
        self.file_patches: dict[str, str] = {}
        self._load_charmmff(tarfilename=tarfilename, skip_streams=skip_streams)
        # self._report()
        self._initialize_resi_to_topfile_map()
        self.provisioned = False
        """ Items below are created by provisioning at run-time """
        self.residues = CharmmResiDict({})
        self.patches = CharmmResiDict({})
        self.pdbrepository = None

    def _load_charmmff(self, tarfilename='toppar_c36_jul24.tgz', skip_streams=['misc', 'cphmd']):
        """ 
        Load the CHARMM force field tarball from the specified path.

        Parameters
        ----------
        tarfilename : str, optional
            The name of the tarball file containing the CHARMM force field files. Default is 'toppar_c36_jul24.tgz'.
        skip_streams : list of str, optional
            A list of stream names to skip when loading the CHARMM force field content. Default is ['misc', 'cphmd'].

        Raises
        -------
        FileNotFoundError
            If the specified tarball file does not exist in the CHARMM force field path.
        """
        def okfilename(name):
            """ Check if a filename is ok to use in the tarfile """
            extension_whitelist = ['.rtf', '.prm', '.str']
            name_blacklist = ['history', 'all22', 'ljpme']
            return not any(blacklist in name for blacklist in name_blacklist) and any(name.endswith(ext) for ext in extension_whitelist)
        tar_path = self.charmmff_path / tarfilename
        if not tar_path.exists():
            raise FileNotFoundError(f'CHARMM force field tarball {self.tarfilename} not found in {self.charmmff_path.name}')
        self.tarfilename = tarfilename  
        logger.debug(f'Loading CHARMM force field tarball {self.tarfilename} from {self.charmmff_path.name}...')
        self.toppar_fs = TarBytesFS.from_file(tar_path, compression='gzip')
        root_listing = [x['name'] for x in self.toppar_fs.ls('toppar') if okfilename(x['name'])]
        # self.contents = {}
        par    = {os.path.basename(x): x for x in root_listing if okfilename(x) and CHARMMFFContent.charmmff_filetype(x) == 'par'}
        top    = {os.path.basename(x): x for x in root_listing if okfilename(x) and CHARMMFFContent.charmmff_filetype(x) == 'top'}
        toppar = {os.path.basename(x): x for x in root_listing if okfilename(x) and CHARMMFFContent.charmmff_filetype(x) == 'toppar'}
        self.fs_resolver = {x: os.path.join('toppar', x) for x in par.keys()}
        self.fs_resolver.update({x: os.path.join('toppar', x) for x in top.keys()})
        self.fs_resolver.update({x: os.path.join('toppar', x) for x in toppar.keys()})

        self.massdict = CharmmMassDict({})

        stream_listing = [x['name'] for x in self.toppar_fs.ls('toppar/stream')]
        self.streams = [os.path.basename(x) for x in stream_listing if os.path.basename(x) not in skip_streams]
        self.streamfiles = {}
        for stream in self.streams:
            streamdir_listing = [x['name'] for x in self.toppar_fs.ls(f'toppar/stream/{stream}')]
            streamfiles = {os.path.basename(x): x for x in streamdir_listing if okfilename(x) and CHARMMFFContent.charmmff_filetype(x) == 'toppar'}
            self.fs_resolver.update({x: os.path.join('toppar/stream', stream, x) for x in streamfiles.keys()})
            toppar.update({os.path.basename(x): x for x in streamfiles.values()})

        self.filenamemap = {'par': par, 'top': top, 'toppar': toppar}
        self.all_topology_files = {x: v for x, v in self.filenamemap['top'].items()}
        self.all_topology_files.update({x: v for x, v in self.filenamemap['toppar'].items()})
        self.all_parameter_files = {x: v for x, v in self.filenamemap['par'].items()}
        self.all_parameter_files.update({x: v for x, v in self.filenamemap['toppar'].items()})
        self.custom_files = []
        self.custom_folder = self.charmmff_path / 'custom'
        if 'custom' in self.charmm_elements and self.custom_folder.exists():
            self._load_custom_files()
        self.patch_folder = self.charmmff_path / 'patches'
        if 'patches' in self.charmm_elements and self.patch_folder.exists():
            self._load_unified_patches(self.patch_folder)
        self.charmmstreamid = {f: CHARMMFFStreamID(f) for f in self.filenamemap['top'].keys()}
        self.charmmstreamid.update({f: CHARMMFFStreamID(f) for f in self.filenamemap['toppar'].keys()})

        logger.debug(f'Loaded {sum(len(v) for v in self.filenamemap.values())} files from CHARMM force field tarball {tarfilename}; subdir-streams: {self.streams}')

        for filetype in ['par', 'top', 'toppar']:
            logger.debug(f'  {filetype}: {len(self.filenamemap[filetype])} files:')
            for keyname, fullname in self.filenamemap[filetype].items():
                logger.debug(f'    {keyname} -> {fullname}')

        for shortname, fullname in self.all_topology_files.items():
            try:
                name_in_tarball = self.fs_resolver[shortname]
                with self.toppar_fs.open(name_in_tarball) as f:
                    contents = f.read().decode()
            except KeyError: # shortname is not in fs_resolver
                with open(fullname, 'r') as f:
                    contents = f.read()
            logger.debug(f'Extracting atom masses from topology file {shortname} ({len(contents)} bytes)')
            masses = CharmmMassList.from_cardlist(extract_mass_lines(contents)).to_dict()
            self.massdict.update(masses)

    def _load_custom_files(self):
        for f in self.custom_folder.iterdir():
            ext = CHARMMFFContent.charmmff_filetype(f.name)
            match ext:
                case 'par':
                    self.filenamemap['par'][f.name] = str(f)
                    self.custom_files.append(f.name)
                    self.all_parameter_files[f.name] = str(f)
                case 'top':
                    self.filenamemap['top'][f.name] = str(f)
                    self.all_topology_files[f.name] = str(f)
                    self.custom_files.append(f.name)
                case 'toppar':
                    self.filenamemap['toppar'][f.name] = str(f)
                    self.all_topology_files[f.name] = str(f)
                    self.all_parameter_files[f.name] = str(f)
                    self.custom_files.append(f.name)
                case _:
                    logger.debug(f'I do not recognize custom CHARMM file {f} in {self.custom_folder.name}')

    def _load_unified_patches(self, patch_folder: Path):
        for f in patch_folder.iterdir():
            if f.suffix in ['.patch', '.diff']:
                with open(f, 'r') as pf:
                    patchtext = pf.read()
                self.file_patches[f.stem] = patchtext

    def _initialize_resi_to_topfile_map(self):
        self.resi_to_topfile_map = {}
        for shortname, fullname in self.all_topology_files.items():
            try:
                name_in_tarball = self.fs_resolver[shortname]
                with self.toppar_fs.open(name_in_tarball) as f:
                    lines = f.read().decode().splitlines()
            except KeyError: # shortname is not in self.fs_resolver
                if not os.path.exists(fullname):
                    raise FileNotFoundError(f'File {fullname} not found in any CHARMM force field content')
                with open(fullname,'r') as f:
                    lines = f.read().splitlines()
            for line in lines:
                if line.startswith('RESI') or line.startswith('PRES'):
                    resi_name = line.split()[1]
                    self.resi_to_topfile_map[resi_name] = shortname

    def _report(self):
        logger.debug(f'Filename map:')
        my_logger(self.filenamemap, logger.debug)
        logger.debug(f'Unified patches:')
        my_logger(self.file_patches, logger.debug)

    @staticmethod
    def charmmff_filetype(filename: str) -> str | None:
        if filename.endswith('.prm'):
            return 'par'
        if filename.endswith('.rtf') or filename.endswith('.top'):
            return 'top'
        if filename.endswith('.str'):
            return 'toppar'
        # logger.debug(f'File {filename} is not a CHARMM force field file')
        return None

    def provision_pdbrepository(self, force_rebuild: bool = False, resnames: list[str] = []):
        self.pdbrepository = PDBRepository(os.path.join(self.charmmff_path, 'pdbrepository'), resnames=resnames, force_rebuild=force_rebuild)

    def provision_residueobjects(self, force_rebuild: bool = False, resnames: list[str] = []):
        is_custom = len(resnames) > 0
        if is_custom:
            logger.debug(f'Provisioning CHARMMFFContent with selected residues/patches:')
            my_logger(resnames, logger.debug)
            self.find_resis_and_patches(resnames=resnames)
        else:
            logger.debug(f'Provisioning CHARMMFFContent with all residues/patches')
            self.resitopcollection = CHARMMFFResiTopCollection(self.charmmff_path, resnames=resnames, force_rebuild=force_rebuild)
        # shortcuts
            self.residues.update(self.resitopcollection.residues)
            self.patches.update(self.resitopcollection.patches)

    def provision(self, force_rebuild: bool = False, resnames: list[str] = []):
        if self.provisioned:
            return
        self.provision_pdbrepository(force_rebuild=force_rebuild, resnames=resnames)
        self.provision_residueobjects(force_rebuild=force_rebuild, resnames=resnames)
        self.provisioned = True

    def deprovision(self):
        if not self.provisioned:
            return
        self.pdbrepository = None
        self.residues.clear()
        self.patches.clear()
        self.provisioned = False

    def get_filename(self, shortname):
        ext = CHARMMFFContent.charmmff_filetype(shortname)
        if ext is not None:
            return self.filenamemap[ext].get(shortname, None)
        return None

    def find_resis_and_patches(self, resnames: list[str] = []):
        """ 
        Find all residues in the CHARMM force field content and associate each with its topology file.
        This function scans all topology files for lines that start with ``RESI`` or ``PRES`` and extracts the residue names.
        It creates a dictionary mapping residue names to the topology files they are found in.
        The residues are stored in the :attr:`~CHARMMFFContent.residues` attribute and the patches in the
        :attr:`~CHARMMFFContent.patches` attribute.
        """
        logger.debug(f'Resnames {resnames}')
        for shortname, fullname in self.all_topology_files.items():
            try:
                name_in_tarball = self.fs_resolver[shortname]
                with self.toppar_fs.open(name_in_tarball) as f:
                    contents = f.read().decode()
            except KeyError: # shortname is not in fs_resolver
                with open(fullname, 'r') as f:
                    contents = f.read()
            # logger.debug(f'Processing topology file {shortname} ({len(contents)} bytes) for residue and patch objects...')
            charmmstreamid = CHARMMFFStreamID(shortname)
            blocks = extract_resi_pres_blocks(contents)
            resi, pres = CharmmResiDict.from_blockstring_list(blocks, metadata=dict(streamID=charmmstreamid.streamID, substreamID=charmmstreamid.substreamID, charmmfftopfile=shortname), resnames=resnames).to_resi_pres()
            self.residues.update(resi)
            self.patches.update(pres)
            # logger.debug(f' -> resis ({len(resi)}) {[r for r in resi.keys()]}')
            # logger.debug(f' -> pres ({len(pres)}) {[p for p in pres.keys()]}')
        self.residues.tally_masses(self.massdict)
        logger.debug(f'Processed {len(self.residues)} residues and {len(self.patches)} patches in CHARMM force field content')

    def get_topfile_of_resname(self, resname: str) -> str | None:
        """
        Given a residue name, return the name of the CHARMMFF topo or stream file that defines it
        """
        first_try = self.resi_to_topfile_map[resname] if resname in self.resi_to_topfile_map else None
        if first_try is None:
            resname = Labels.charmm_resname_of_pdb_resname.get(resname, resname)
            return self.resi_to_topfile_map[resname] if resname in self.resi_to_topfile_map else None
        return first_try

    def copy_charmmfile_local(self, basename: str) -> str:
        """
        Copy a NAMD-friendly version of a CHARMMFF file to the local directory.
        This function checks if the file already exists in the current working directory.
        If it does, it returns the basename. If the file is found in the custom files, it copies it from there.
        If the file is found in the tarball or any custom directory, it extracts it and writes it to the local directory, filtering out CHARMM commands that give NAMD trouble.
        If the file is not found in either location, it logs a warning.

        Parameters
        ----------
        basename : str
            The basename of the CHARMM file to copy. This should be a recognizable CHARMMFF rtf, prm, or str file name without any directory path.

        Returns
        -------
        str
            The basename of the copied file in the local directory.
        """
        # When copying a parameter file into a NAMD run directory, lines that begin with these keywords are removed
        comment_these_out = ['set', 'if', 'WRNLEV', 'BOMLEV', 'return', 'endif']

        if os.path.exists(basename):
            # logger.debug(f'{basename} already exists in {os.getcwd()}')
            return basename
        if os.sep in basename:
            # this is a path
            # logger.debug(f'expected a basename and got a path {basename}')
            raise ValueError(f'Expected a basename, but got a path: {basename}')
            # logger.debug(f'truncated to basename {basename}')
        ext = CHARMMFFContent.charmmff_filetype(basename)
        if basename in self.custom_files:
            with open(self.filenamemap[ext][basename]) as file:  # custom files are not part of the cache?
                lines = file.read().splitlines()
                # logger.debug(f'found {len(lines)} lines in {basename} in custom files')
                with open(basename, 'w') as f:
                    for l in lines:  # l will contain the newline character
                        is_comment = any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l + '\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}' + '\n')
        elif basename in self.filenamemap[ext]:
            logger.debug(f'found {basename} in at {self.filenamemap[ext][basename]} in tarball')
            stem, dum = os.path.splitext(basename)
            longname = self.filenamemap[ext][basename]
            with self.toppar_fs.open(self.fs_resolver[basename]) as f:
                # logger.debug(f'Opening {longname} in tarball')
                content = f.read().decode()
                if 'cholesterol' in longname:  # the cholesterol substream has two models, and it specifies the first one by default
                    # we will parse the conditional script to get the correct model
                    parsed_content_dict = parse_conditional_script(content)
                    parsed_content = parsed_content_dict['parsed']
                    content = parsed_content
                if stem in self.file_patches:
                    patchtext = self.file_patches[stem]
                    logger.debug(f'Applying unified patch {stem} to {basename}')
                    content = apply_unified_diff(content, patchtext)
            lines = content.splitlines()
            # logger.debug(f'type of lines is {type(lines)}')
            # logger.debug(f'found {len(lines)} lines in {basename} in tarfile')
            with open(basename, 'w') as f:
                for l in lines:  # l will NOT contain the newline character
                    is_comment = any([l.startswith(x) for x in comment_these_out])
                    if not is_comment:
                        f.write(l + '\n')
                    else:
                        f.write('! commented out by pestifer:\n')
                        f.write(f'! {l}\n')
        else:
            logger.warning(f'copy_charmmfile_local: {basename} not found in charmmff')
        return basename

    def add_custom_directory(self, user_custom_directory: str | Path):
        """ 
        Add a user custom directory to the :class:`CHARMMFFContent`.
        This directory should contain custom files that can be used in addition to the standard CHARMM force field files.

        Parameters
        ----------
        user_custom_directory : str
            The path to the user custom directory containing additional CHARMM files.

        Raises
        -------
        NotADirectoryError
            If the specified path is not a directory.
        """
        if not os.path.isdir(user_custom_directory):
            raise NotADirectoryError(f'Expected a directory at {user_custom_directory}, but it is not a directory')
        logger.debug(f'Adding user custom directory {user_custom_directory} to CHARMMFFContent')
        self.custom_files.extend(os.listdir(user_custom_directory))
        for f in self.custom_files:
            assert f not in self.filenamemap, f'user custom file {f} already exists in filenamemap'
            ext = CHARMMFFContent.charmmff_filetype(f)
            self.filenamemap[ext][f] = os.path.join(user_custom_directory, f)

    def clean_local_charmmff_files(self):
        """ 
        Remove all local CHARMM force field files.
        """
        for f in os.listdir('.'):
            logger.debug(f'Checking file {f} for removal')
            if CHARMMFFContent.charmmff_filetype(f):
                os.remove(f)
    
    def get_resi(self, resname: str) -> CharmmResi | None:
        if not self.provisioned:
            raise RuntimeError('CHARMMFFContent must be provisioned before accessing residues or patches')
        return self.residues.get_residue(resname)

    def get_pres(self, presname: str) -> CharmmResi | None:
        if not self.provisioned:
            raise RuntimeError('CHARMMFFContent must be provisioned before accessing residues or patches')
        return self.patches.get_residue(presname)

    def __contains__(self, resi_or_pres_name: str) -> bool:
        if not self.provisioned:
            raise RuntimeError('CHARMMFFContent must be provisioned before accessing residues or patches')
        return resi_or_pres_name in self.residues or resi_or_pres_name in self.patches
    
    def checkout_pdb(self, name: str) -> PDBInput | None:
        if not self.provisioned:
            raise RuntimeError('CHARMMFFContent must be provisioned before accessing PDB repository')
        return self.pdbrepository.checkout(name)

class CHARMMFFStreamID:
    """
    A class for handling the CHARMM force field stream ID and substream ID.
    
    This class parses the filename of a CHARMM force field file to extract the stream ID and substream ID.
    
    Parameters
    ----------
    charmmff_filename : str
        The name of the CHARMM force field file.

    Attributes
    ----------
    charmmff_filename : str
        The name of the CHARMM force field file.
    streamID : str
        The stream ID extracted from the filename.
    substreamID : str
        The substream ID extracted from the filename, if applicable.
    """
    def __init__(self, charmmff_filename: str):
        self.streamID = ''
        self.substreamID = ''
        self.charmmff_filename = os.path.basename(charmmff_filename)
        pref, ext = os.path.splitext(self.charmmff_filename)
        if ext == '.prm':
            tokens = pref.split('_')
            if tokens[0] == 'par' and (tokens[1] == 'all35' or tokens[1] == 'all36'):
                self.streamID = tokens[2]
                self.substreamID = ''
        elif ext == '.rtf':
            tokens = pref.split('_')
            if tokens[0] == 'top' and (tokens[1] == 'all35' or tokens[1] == 'all36'):
                self.streamID = tokens[2]
                self.substreamID = ''
        elif ext == '.str':
            tokens = pref.split('_')
            if tokens[0] == 'toppar':
                if len(tokens) == 2:
                    self.streamID = tokens[1]
                    self.substreamID = ''
                elif len(tokens) == 3:
                    self.streamID = '_'.join(tokens[1:3])
                    if self.streamID == 'all36_moreions':
                        # this is a special case for the all36_moreions stream
                        self.streamID = 'water_ions'
                        self.substreamID = ''
                    self.substreamID = ''
                elif len(tokens) >= 4:
                    self.streamID = tokens[2]
                    self.substreamID = '_'.join(tokens[3:])
        logger.debug(f'CHARMMFFStreamID: parsed {self.charmmff_filename} to streamID={self.streamID}, substreamID={self.substreamID}')


