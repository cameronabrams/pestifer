# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This module defines classes that facilitate the handling of CHARMM force field content.
"""
from dataclasses import dataclass
import logging
import os
import re

from pathlib import Path

from .charmmfftop import CharmmMassDict, CharmmMassList, CharmmResiDict, CharmmResi
from .pdbrepository import PDBRepository, PDBRepositoryData, PDBInput

from ..util.cacheable_object import CacheableObject, TarBytesFS

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

# @dataclass
class CHARMMFFContentData(CacheableObject):
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
    tarfile : TarBytesFS
        TarBytesFS object representing the CHARMM force field tarball.
    filenamemap : dict
        Maps file basenames to their full paths in the CHARMM force field content.
    custom_files : list
        List of custom files that can be added to the CHARMM force field content.
    all_topology_files : list
        List of all topology files in the CHARMM force field content.
    residues: CharmmResiDict
        Maps residue names to their corresponding CharmmResi objects
    patches: CharmmResiDict
        Maps patch names to their corresponding CharmmResi objects.
    """
    def _build_from_resources(self, charmmff_path: str = '.', tarfilename: str = 'toppar_c36_jul24.tgz'):
        """ Method to build the CHARMMFFContent object from resources, if the cache is stale. """
        self.tarfile = None
        self.tarfilename = tarfilename
        self.filenamemap = {}
        if not os.path.isdir(charmmff_path):
            raise NotADirectoryError(f'Expected a directory at {charmmff_path}, but it is not a directory')
        self.charmmff_path = os.path.abspath(charmmff_path)
        self.basename = os.path.basename(self.charmmff_path)
        self.parent_path = os.path.dirname(self.charmmff_path)
        cwd = os.getcwd()
        os.chdir(self.parent_path)
        self.dirtree = {a: [b, c] for a, b, c in os.walk(self.basename)}
        os.chdir(cwd)
        self.charmm_elements = self.dirtree[self.basename][0]
        logger.debug(f'Members of {self.charmmff_path}: {self.charmm_elements}')
        logger.debug(f'Initializing pdb repository for CHARMMFFContent at {self.charmmff_path}...')
        self.pdbrepository = PDBRepositoryData()
        if 'pdbrepository' in self.charmm_elements:
            os.chdir(os.path.join(self.charmmff_path, 'pdbrepository'))
            members = os.listdir('.')
            for m in members:
                self.pdbrepository.add_path(m)
            os.chdir(cwd)
        self.load_charmmff(tarfilename)
        logger.debug(f'Parsing all RESI and PRES objects...')
        self.all_topology_files = {x: v for x, v in self.filenamemap.items() if x.endswith('.str') or x.endswith('.rtf') or x.endswith('.top')}
        self.residues = CharmmResiDict({})
        self.patches = CharmmResiDict({})
        self.find_resis_and_patches()
        logger.debug(f'CHARMMFFContent initialized with {len(self.residues)} residues and {len(self.patches)} patches.')

    def load_charmmff(self, tarfilename='toppar_c36_jul24.tgz', skip_streams=['misc', 'cphmd']):
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
            return not 'history' in name and not 'all22' in name and not 'ljpme' in name and (name.endswith('.str') or name.endswith('.prm') or name.endswith('.rtf'))

        if not os.path.exists(os.path.join(self.charmmff_path, tarfilename)):
            raise FileNotFoundError(f'CHARMM force field tarball {tarfilename} not found in {self.charmmff_path}')
        logger.debug(f'Loading CHARMM force field tarball {tarfilename} from {self.charmmff_path}')
        self.tarmembers = []
        self.tarfile = TarBytesFS.from_file(os.path.join(self.charmmff_path, tarfilename), compression='gzip')
        self.toplevel_listing = [x['name'] for x in self.tarfile.ls('toppar')]
        self.contents = {}
        self.toplevel_par = {os.path.basename(x): x for x in self.toplevel_listing if x.startswith('toppar/par_') and okfilename(x)}
        self.toplevel_top = {os.path.basename(x): x for x in self.toplevel_listing if x.startswith('toppar/top_') and okfilename(x)}
        self.toplevel_toppar = {os.path.basename(x): x for x in self.toplevel_listing if x.startswith('toppar/toppar_') and okfilename(x)}
        self.filenamemap = {**self.toplevel_par, **self.toplevel_top, **self.toplevel_toppar}
        # logger.debug(f'filenamemap: {self.filenamemap}')
        self.stream_listing = [x['name'] for x in self.tarfile.ls('toppar/stream')]
        self.streams = [os.path.basename(x) for x in self.stream_listing if os.path.basename(x) not in skip_streams]
        self.streamfiles = {}
        for stream in self.streams:
            streamdir_listing = [x['name'] for x in self.tarfile.ls(f'toppar/stream/{stream}')]
            self.streamfiles[stream] = {os.path.basename(x): x for x in streamdir_listing if okfilename(x)}
            self.filenamemap.update(self.streamfiles[stream])
        for name, fullname in self.filenamemap.items():
            # logger.debug(f'Attempting open of {name}')
            with self.tarfile.open(fullname) as f:
                self.contents[name] = f.read().decode()
                if 'cholesterol' in fullname:  # the cholesterol substream has two models, and it specifies the first one by default
                # we will parse the conditional script to get the correct model
                    parsed_content_dict = parse_conditional_script(self.contents[name])
                    parsed_content = parsed_content_dict['parsed']
                    self.contents[name] = parsed_content
        check_basenames = list(self.filenamemap.keys())
        assert len(check_basenames) == len(set(check_basenames)), f'found duplicate basenames in charmmff tarball: {check_basenames}'
        logger.debug(f'Loaded {len(self.filenamemap)} files from CHARMM force field tarball {tarfilename}; subdir-streams: {self.streams}')
        self.custom_files = []
        if 'custom' in self.charmm_elements and os.path.exists(os.path.join(self.charmmff_path, 'custom')):
            self.custom_files = os.listdir(os.path.join(self.charmmff_path, 'custom'))
            # logger.debug(f'Found {len(self.custom_files)} custom files in {os.path.join(self.charmmff_path,"custom")}: {self.custom_files}')
        self.non_charmmff_custom_files = []
        for f in self.custom_files:
            logger.debug(f'Checking custom file {f}')
            if CHARMMFFContentData.is_charmmff_file(f):
                # logger.debug(f'Adding custom file {f} to CHARMMFFContentData filenamemap')
                self.filenamemap[f] = os.path.join(self.charmmff_path, 'custom', f)
                with open(self.filenamemap[f]) as custom_file:
                    self.contents[f] = custom_file.read()
            else:
                self.non_charmmff_custom_files.append(f)
        for nccf in self.non_charmmff_custom_files:
            self.custom_files.remove(nccf)

        total_bytes = 0
        for x, v in self.contents.items():
            total_bytes += len(v.encode())
        logger.debug(f'Total bytes of CHARMMFF content: {total_bytes}')

    @staticmethod
    def is_charmmff_file(filename: str) -> bool:
        if filename.startswith('par') or filename.startswith('top_') or \
           filename.startswith('toppar') or filename.startswith('charmm') or \
           filename.endswith('.str') or filename.endswith('.prm') or \
           filename.endswith('.rtf') or filename.endswith('.top'):
            return True
        # logger.debug(f'File {filename} is not a CHARMM force field file')
        return False

    def find_resis_and_patches(self):
        """ 
        Find all residues in the CHARMM force field content and associate each with its topology file.
        This function scans all topology files for lines that start with ``RESI`` or ``PRES`` and extracts the residue names.
        It creates a dictionary mapping residue names to the topology files they are found in.
        The residues are stored in the :attr:`~CHARMMFFContent.residues` attribute and the patches in the
        :attr:`~CHARMMFFContent.patches` attribute.
        """
        self.residues = CharmmResiDict({})
        self.patches = CharmmResiDict({})
        self.massdict = CharmmMassDict({})
        for topfile in self.all_topology_files.keys():
            logger.debug(f'Processing topology file {topfile} for residue and patch objects...')
            charmmstreamid = CHARMMFFStreamID(topfile)
            blocks = extract_resi_pres_blocks(self.contents[topfile])
            resi, pres = CharmmResiDict.from_blockstring_list(blocks, metadata=dict(streamID=charmmstreamid.streamID, substreamID=charmmstreamid.substreamID, charmmfftopfile=topfile)).to_resi_pres()
            masses = CharmmMassList.from_cardlist(extract_mass_lines(self.contents[topfile])).to_dict()
            self.residues.update(resi)
            self.patches.update(pres)
            self.massdict.update(masses)
            logger.debug(f' -> resis ({len(resi)}) {[r for r in resi.keys()]}')
            logger.debug(f' -> pres ({len(pres)}) {[p for p in pres.keys()]}')
        logger.debug(f'Found {len(self.residues)} residues, {len(self.patches)} patches, and {len(self.massdict)} masses in CHARMM force field content')

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

class CHARMMFFContent(CHARMMFFContentData):
    """
    A class for handling the content of CHARMM force field files.
    self, charmmff_path: str = '.', tarfilename: str = 'toppar_c36_jul24.tgz'
    """
    def __init__(self, *args, **kwargs):
        if len(args) == 1 and isinstance(args[0], CHARMMFFContentData):
            self.__dict__.update(args[0].__dict__)
        elif len(args) == 1 and isinstance(args[0], str):
            super().__init__(*args, **kwargs)

    def get_topfile_of_resname(self, resname: str) -> str | None:
        """
        Given a residue name, return the top file that contains it
        """
        resi = self.residues.get(resname, None)
        if resi is not None:
            return resi.metadata['charmmfftopfile']
        else:
            logger.warning(f'Residue {resname} not found in CHARMM force field content')
            return None

    def copy_charmmfile_local(self, basename):
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
            basename = os.path.split(basename)[1]
            # logger.debug(f'truncated to basename {basename}')
        if basename in self.custom_files:
            with open(self.filenamemap[basename]) as file:  # custom files are not part of the cache?
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
        elif basename in self.filenamemap:
            logger.debug(f'found {basename} in at {self.filenamemap[basename]} in tarball')
            lines = self.contents[basename].splitlines()
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
            self.filenamemap[f] = os.path.join(user_custom_directory, f)
            with open(self.filenamemap[f]) as custom_file:
                self.contents[f] = custom_file.read()

    def clean_local_charmmff_files(self):
        """ 
        Remove all local CHARMM force field files that start with ``par``, ``top_``, ``toppar``, ``charmm``, or end with ``.str``, ``.prm``, or ``.rtf``.
        This function is useful for cleaning up the local directory where CHARMM files are stored.
        It will remove files that match the specified patterns, ensuring that only relevant CHARMM files are kept.
        """
        for f in os.listdir('.'):
            # logger.debug(f'Checking file {f} for removal')
            if self.is_charmmff_file(f):
                os.remove(f)

    def get_topfile_of_resname(self, resname: str) -> str | None:
        """ 
        Given a residue name, return the top file that contains it 
        """
        resi = self.residues.get(resname, None)
        if resi is not None:
            return resi.metadata['charmmfftopfile']
        else:
            logger.warning(f'Residue {resname} not found in CHARMM force field content')
            return None

    def get_topfile_of_patchname(self, patchname: str) -> str | None:
        """ 
        Given a patch name, return the top file that contains it 
        """
        patch = self.patches.get(patchname, None)
        if patch is not None:
            return patch.metadata['charmmfftopfile']
        else:
            logger.warning(f'Patch {patchname} not found in CHARMM force field content')
            return None
    
    def get_resi(self, resname: str) -> CharmmResi | None:
        return self.residues.get_residue(resname)

    def get_pres(self, presname: str) -> CharmmResi | None:
        return self.patches.get_residue(presname)

    def __contains__(self, resi_or_pres_name: str) -> bool:
        return resi_or_pres_name in self.residues or resi_or_pres_name in self.patches
    
    def checkout_pdb(self, name: str) -> PDBInput | None:
        pdb_repo = PDBRepository(self.pdbrepository)
        return pdb_repo.checkout(name)