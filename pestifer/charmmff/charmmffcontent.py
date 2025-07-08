# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This module defines classes that facilitate the handling of CHARMM force field content.
It includes functionality to parse CHARMM force field files, extract relevant data such as masses and residues, and manage custom directories for additional CHARMM files.
"""
import logging
import os
import re
import tarfile
from .charmmtop import CharmmMassRecord, CharmmMasses, CharmmTopResi
from ..core.pdbrepository import PDBRepository
from ..core.stringthings import my_logger
from ..core.labels import Labels

logger=logging.getLogger(__name__)

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
        - 'parsed': a string containing the processed script text with comments and conditionals resolved.
        - 'vars': a dictionary of variables defined in the script.
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
    return dict(parsed="\n".join(output),vars=vars)

def extract_resi_pres_blocks(text, keywords=('RESI', 'PRES')):
    """ 
    Extract blocks of text starting with RESI or PRES and ending before the next RESI, PRES, ATOMS, or EOF.
    This function uses a regular expression to find blocks of text that start with the specified keywords
    and continue until the next occurrence of one of the keywords or the end of the file.
    
    Parameters
    ----------
    text : str
        The input text from which to extract the blocks.
    keywords : tuple of str, optional
        The keywords that indicate the start of a block. Default is ('RESI', 'PRES').   
    Returns
    -------
    list of str
        A list of strings, each containing a block of text that starts with one of the specified keywords.
    """
    
    keyword_pattern = '|'.join(re.escape(k) for k in keywords)
    # Include ATOMS and BONDS as additional block-end sentinels
    terminator_pattern = rf'^(?:{keyword_pattern}|ATOMS|BONDS|ANGLES|DIHEDRALS)\b'
    # Match starting with RESI or PRES, ending before the next RESI, PRES, ATOMS, or EOF
    pattern = rf'(?i)^(({keyword_pattern})\s.*?)(?={terminator_pattern}|\Z)'
    matches = re.finditer(pattern, text, flags=re.DOTALL | re.MULTILINE)
    return [m.group(1).strip() for m in matches]

def extract_mass_lines(file_contents):
    """
    Extract lines containing mass information from the CHARMM force field files.
    This function scans the contents of a CHARMM force field file and returns lines that start with "MASS".
    Parameters
    ----------
    file_contents : str
        The contents of the CHARMM force field file to scan for mass information.
    """
    return [line for line in file_contents.splitlines() if line.strip().upper().startswith("MASS")]

# When copying a parameter file into a NAMD run directory, lines that begin with these keywords are removed
comment_these_out=['set','if','WRNLEV','BOMLEV','return','endif']

class CHARMMFFContent:
    """ 
    A class for handling all CHARMM force field content.  
    
    The CHARMM force field is
    stored in a tarball downloaded directly from the MacKerell lab at the University of Michigan.
    https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
    
    Attributes
    ----------
    charmmff_path : str
        The path to the directory containing the CHARMM force field files.
    tarfilename : str
        The name of the tarball file containing the CHARMM force field files.
    tarfile : tarfile.TarFile
        The tarfile object representing the CHARMM force field tarball.
    filenamemap : dict
        A dictionary mapping file basenames to their full paths in the CHARMM force field content.
    custom_files : list
        A list of custom files that can be added to the CHARMM force field content.
    all_topology_files : list
        A list of all topology files in the CHARMM force field content.
    """
    def __init__(self,charmmff_path='',tarfilename='toppar_c36_jul24.tgz',user_custom_directory=None):
        self.tarfile=None
        self.tarfilename=tarfilename
        self.filenamemap={}
        if not charmmff_path:
            charmmff_path='.'
        if not os.path.isdir(charmmff_path):
            raise NotADirectoryError(f'Expected a directory at {charmmff_path}, but it is not a directory')
        self.charmmff_path=os.path.abspath(charmmff_path)
        self.basename=os.path.basename(self.charmmff_path)
        self.parent_path=os.path.dirname(self.charmmff_path)
        cwd=os.getcwd()
        os.chdir(self.parent_path)
        self.dirtree={a:[b,c] for a,b,c in os.walk(self.basename)}
        os.chdir(cwd)
        self.charmm_elements=self.dirtree[self.basename][0]
        logger.debug(f'Members of {self.charmmff_path}: {self.charmm_elements}')
        self.pdbrepository=PDBRepository()
        if 'pdbrepository' in self.charmm_elements:
            os.chdir(os.path.join(self.charmmff_path,'pdbrepository'))
            members=os.listdir('.')
            for m in members:
                self.pdbrepository.add_path(m)
            os.chdir(cwd)
        self.load_charmmff(tarfilename)
        self.custom_files=[]
        if 'custom' in self.charmm_elements and os.path.exists(os.path.join(self.charmmff_path,'custom')):
            self.custom_files=self.dirtree[f'{self.basename}/custom'][1]
        for f in self.custom_files:
            assert f not in self.filenamemap, f'custom file {f} already exists in filenamemap'
            logger.debug(f'Adding custom file {f} to CHARMMFFContent filenamemap ({len(self.filenamemap)} entries before adding)')
            self.filenamemap[f]=os.path.join(self.charmmff_path,'custom',f)
        if user_custom_directory is not None:
            if not os.path.isdir(user_custom_directory):
                raise NotADirectoryError(f'Expected a directory at {user_custom_directory}, but it is not a directory')
            logger.debug(f'Adding user custom directory {user_custom_directory} to CHARMMFFContent')
            self.custom_files.extend(os.listdir(user_custom_directory))
            for f in self.custom_files:
                assert f not in self.filenamemap, f'user custom file {f} already exists in filenamemap'
                self.filenamemap[f]=os.path.join(user_custom_directory,f)
        self.all_topology_files=[x for x in self.filenamemap.values() if x.endswith('.str') or x.endswith('.rtf') or x.endswith('.top')]
        self.find_patches()
        self.find_resis()
        logger.debug(f'filename map:')
        for k,v in self.filenamemap.items():
            logger.debug(f'  {k} -> {v}')

    def __del__(self):
        """ 
        Close the tarfile if it is open 
        """
        if self.tarfile is not None:
            self.tarfile.close()
            self.tarfile=None
            logger.debug('Closed CHARMM force field tarfile')
        else:
            logger.debug('No CHARMM force field tarfile to close')
        del self.pdbrepository

    def add_custom_directory(self,user_custom_directory):
        """ 
        Add a user custom directory to the CHARMMFFContent.
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
            self.filenamemap[f]=os.path.join(user_custom_directory,f)

    def load_charmmff(self,tarfilename='toppar_c36_jul24.tgz',skip_streams=['misc','cphmd']):
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

        if not os.path.exists(os.path.join(self.charmmff_path,tarfilename)):
            raise FileNotFoundError(f'CHARMM force field tarball {tarfilename} not found in {self.charmmff_path}')
        logger.debug(f'Loading CHARMM force field tarball {tarfilename} from {self.charmmff_path}')
        self.tarmembers=[]
        if self.tarfile is not None:
            self.tarfile.close()
            self.tarfile=None
        self.tarfile=tarfile.open(os.path.join(self.charmmff_path,tarfilename),'r:gz')
        self.tarmembers=self.tarfile.getmembers()
        self.toplevel_par={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/par_') and okfilename(x.name)}
        self.toplevel_top={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/top_') and okfilename(x.name)}
        self.toplevel_toppar={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/toppar_') and okfilename(x.name)}
        self.filenamemap={**self.toplevel_par,**self.toplevel_top,**self.toplevel_toppar}
        logger.debug(f'filenamemap: {self.filenamemap}')
        self.streams=[os.path.basename(x.name) for x in self.tarmembers 
                      if x.isdir() and x.name.startswith('toppar/stream/') and os.path.basename(x.name) not in skip_streams]
        self.streamfiles={}
        for stream in self.streams:
            self.streamfiles[stream]={os.path.basename(x.name):x.name 
                                      for x in self.tarmembers if x.isfile() 
                                      and x.name.startswith(f'toppar/stream/{stream}/') and okfilename(x.name)}
            self.filenamemap.update(self.streamfiles[stream])
        check_basenames=list(self.filenamemap.keys())
        assert len(check_basenames)==len(set(check_basenames)),f'found duplicate basenames in charmmff tarball: {check_basenames}'
        logger.debug(f'Loaded {len(self.filenamemap)} files from CHARMM force field tarball {tarfilename}; subdir-streams: {self.streams}')
        

    def copy_charmmfile_local(self,basename):
        """ 
        Copy a CHARMM file to the local directory.
        This function checks if the file already exists in the current working directory.
        If it does, it returns the basename. If the file is found in the custom files, it copies it from there.
        If the file is found in the tarball or any custom directory, it extracts it and writes it to the local directory, filtering out CHARMM commands that give NAMD trouble.
        If the file is not found in either location, it logs a warning.

        Parameters
        ----------
        basename : str
            The basename of the CHARMM file to copy. This should be a file name without any directory path.

        Returns
        -------
        str
            The basename of the copied file in the local directory.
        """
        if os.path.exists(basename):
            logger.debug(f'{basename} already exists in {os.getcwd()}')
            return basename
        if os.sep in basename:
            # this is a path
            logger.debug(f'expected a basename and got a path {basename}')
            basename=os.path.split(basename)[1]
            logger.debug(f'truncated to basename {basename}')
        if basename in self.custom_files:
            with open(self.filenamemap[basename]) as file:
                lines=file.read().splitlines()
                logger.debug(f'found {len(lines)} lines in {basename} in custom files')
                with open(basename,'w') as f:
                    for l in lines: # l will contain the newline character
                        is_comment=any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l+'\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}'+'\n')
        elif basename in self.filenamemap:
            logger.debug(f'found {basename} in at {self.filenamemap[basename]} in tarball')
            with self.tarfile.extractfile(self.filenamemap[basename]) as file:
                logger.debug(f' file has type {type(file)}')
                lines=file.read().decode().splitlines()
                logger.debug(f'type of lines is {type(lines)}')
                logger.debug(f'found {len(lines)} lines in {basename} in tarfile')
                with open(basename,'w') as f:
                    for l in lines: # l will NOT contain the newline character
                        is_comment=any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l+'\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}\n')
        else:
            logger.warning(f'copy_charmmfile_local: {basename} not found in charmmff')
        return basename
    
    def clean_local_charmmff_files(self):
        """ 
        Remove all local CHARMM force field files that start with 'par', 'top', 'toppar', 'charmm', or end with '.str', '.prm', or '.rtf'.
        This function is useful for cleaning up the local directory where CHARMM files are stored.
        It will remove files that match the specified patterns, ensuring that only relevant CHARMM files are kept.
        """
        for f in os.listdir('.'):
            if f.startswith('par') or f.startswith('top_') or f.startswith('toppar') or f.startswith('charmm') or f.endswith('.str') or f.endswith('.prm') or f.endswith('.rtf'):
                os.remove(f)
    
    def lines_from_topfile(self,topfile):
        """ 
        Extract the lines from a top file.
        This function reads the contents of a top file, either from the tarfile or from a mapped filename in the filenamemap.

        Parameters
        ----------
        topfile : str
            The name of the top file to extract lines from. This can be a full path or just the basename.

        Returns
        -------
        list of str
            A list of lines extracted from the top file. If the file is not found, an empty list is returned.
        """
        lines=self.contents_from_topfile(topfile).splitlines()
        return lines

    def contents_from_topfile(self,topfile):
        """ 
        Extract the contents of a top file.
        This function reads the contents of a top file, either from the tarfile or from a mapped filename in the filenamemap.  Applies the logic filter if this is a cholesterol substream.

        Parameters
        ----------
        topfile : str
            The name of the top file to extract contents from. This can be a full path or just the basename.

        Returns
        -------
        str
            The contents of the top file as a string. If the file is not found, an empty string is returned.
        """
        content=''
        for m in self.tarmembers:
            if topfile==m.name:
                logger.debug(f'Found {topfile} in tarfile member {m.name}')
                with self.tarfile.extractfile(m) as f:
                    content=f.read().decode()
                    break
        if not content:
            mapped_name=self.filenamemap.get(os.path.basename(topfile),None)
            if mapped_name is not None:
                logger.debug(f'Extracting lines from {topfile} using mapped name {mapped_name}')
                with open(mapped_name,'r') as f:
                    content=f.read()
            else:
                logger.warning(f'Could not find {topfile} in tarfile or filenamemap')
                return ''
        logger.debug(f'Extracted {len(content)} characters from {topfile}')

        parsed_content = content
        if 'cholesterol' in topfile:  # the cholesteral substream has two models, and it specifies the first one by default
            # we will parse the conditional script to get the correct model
            parsed_content_dict=parse_conditional_script(content)
            parsed_content=parsed_content_dict['parsed']
            logger.debug(f'Parsed {topfile} with conditional script based on {parsed_content_dict["vars"]}')

        return parsed_content
    
    def masses_from_topfile(self,topfile):
        """ 
        Extract the masses from a top file
        This function reads the contents of a top file and extracts the mass records.

        Parameters
        ----------
        topfile : str
            The name of the top file to extract masses from. This can be a full path or just the basename.

        Returns
        -------
        CharmmMasses
            An instance of CharmmMasses containing the mass records extracted from the top file.
        """
        masses=[]
        lines=extract_mass_lines(self.contents_from_topfile(topfile))
        for i in range(len(lines)):
            masses.append(CharmmMassRecord(lines[i]))
        # for m in masses:
        #     logger.debug(f'Found mass record: \'{str(m)}\'')
        return CharmmMasses(masses)

    def resis_from_topfile(self,topfile,metadata={}):
        """ 
        Extract the residues from a top file.
        This function reads the contents of a top file and extracts the residue blocks.
        It looks for blocks that start with 'RESI' or 'PRES' and creates instances of CharmmTopResi for each block.

        Parameters
        ----------
        topfile : str
            The name of the top file to extract residues from. This can be a full path or just the basename.
        metadata : dict, optional
            A dictionary containing metadata to be associated with each residue. Default is an empty dictionary.

        Returns
        -------
        list of CharmmTopResi
            A list of CharmmTopResi objects representing the residues extracted from the top file.
        """
        blocks=extract_resi_pres_blocks(self.contents_from_topfile(topfile))
        R=[]
        for block in blocks:
            # logger.debug(block)
            key=block.split()[0].upper()
            # if key!='RESI':
            #     logger.debug(f'Expected RESI block, but found {key} in {topfile}')
            resi=CharmmTopResi(block,key=key)
            resi.metadata=metadata
            # logger.debug(f'Found residue {resi.resname} in {topfile} with metadata {resi.metadata}')         
            R.append(resi)
        resname_list=[x.resname for x in R]
        if len(resname_list)!=len(set(resname_list)):
            logger.warning(f'Found duplicate residue names in {topfile}: {resname_list}')
        return R

    def find_resis(self):
        """ 
        Find all residues in the CHARMM force field content.
        This function scans all topology files for lines that start with 'RESI' or 'PRES' and extracts the residue names.
        It creates a dictionary mapping residue names to the topology file they are found in.
        The residues are stored in the `self.residues` attribute.
        """
        self.residues={}
        for topfile in self.all_topology_files:
            lines=self.lines_from_topfile(topfile)
            for line in lines:
                if line.upper().startswith('RESI'):
                    resname=line.split()[1]
                    if resname not in self.residues:
                        self.residues[resname]=os.path.basename(topfile)
                    alias=Labels.pdb_resname_of_charmm_resname.get(resname,None)
                    if alias is not None and alias not in self.residues:
                        self.residues[alias]=os.path.basename(topfile)
                        # logger.debug(f'Found residue {resname} in {topfile}')
        logger.debug(f'Found {len(self.residues)} residues in CHARMM force field content')

    def find_patches(self):
        """ 
        Find all patches in the CHARMM force field content.
        This function scans all topology files for lines that start with 'PRES' and extracts the patch names.
        It creates a dictionary mapping patch names to the topology file they are found in.
        The patches are stored in the `self.patches` attribute.
        """
        self.patches={}
        for topfile in self.all_topology_files:
            lines=self.lines_from_topfile(topfile)
            for line in lines:
                if line.upper().startswith('PRES'):
                    patchname=line.split()[1]
                    if patchname not in self.patches:
                        self.patches[patchname]=os.path.basename(topfile)
                        # logger.debug(f'Found patch {patchname} in {topfile}')
        logger.debug(f'Found {len(self.patches)} patches in CHARMM force field content')

    def get_topfile_of_patchname(self,patchname):
        """ 
        Given a patch name, return the top file that contains it """
        return self.patches.get(patchname, None)

    def get_topfile_of_resname(self,resname):
        """ 
        Given a residue name, return the top file that contains it.
        This function searches through all topology files and returns the first one that contains the specified residue name.

        Parameters
        ----------
        resname : str
            The name of the residue to search for.

        Returns
        -------
        str or None
            The name of the topology file containing the residue, or None if not found.
        """
        if resname in self.residues:
            return self.residues[resname]
        else:
            logger.warning(f'Residue {resname} not found in CHARMM force field content')
            return None

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
    charmff_filename : str
        The name of the CHARMM force field file.
    streamID : str
        The stream ID extracted from the filename.
    substreamID : str
        The substream ID extracted from the filename, if applicable.
    """
    def __init__(self,charmmff_filename):
        self.streamID=''
        self.substreamID=''
        self.charmmff_filename=os.path.basename(charmmff_filename)
        pref,ext=os.path.splitext(self.charmmff_filename)
        if ext=='.prm':
            tokens=pref.split('_')
            if tokens[0]=='par' and (tokens[1]=='all35' or tokens[1]=='all36'):
                self.streamID=tokens[2]
                self.substreamID=''
        elif ext=='.rtf':
            tokens=pref.split('_')
            if tokens[0]=='top' and (tokens[1]=='all35' or tokens[1]=='all36'):
                self.streamID=tokens[2]
                self.substreamID=''
        elif ext=='.str':
            tokens=pref.split('_')
            if tokens[0]=='toppar':
                if len(tokens)==2:
                    self.streamID=tokens[1]
                    self.substreamID=''
                elif len(tokens)==3:
                    self.streamID='_'.join(tokens[1:3])
                    if self.streamID=='all36_moreions':
                        # this is a special case for the all36_moreions stream
                        self.streamID='water_ions'
                        self.substreamID=''
                    self.substreamID=''
                elif len(tokens)>=4:
                    self.streamID=tokens[2]
                    self.substreamID='_'.join(tokens[3:])
        logger.debug(f'CHARMMFFStreamID: parsed {self.charmmff_filename} to streamID={self.streamID}, substreamID={self.substreamID}')

class CHARMMFFResiDatabase:
    """ 
    A class for handling a database of CHARMM force field residues.
    This class loads residues and patches from CHARMM force field topology/stream files, including toplevel files and specific streams.
    It maintains a collection of residues, patches, and masses, and provides methods to access and manipulate this data.
    
    Attributes
    ----------
    charmmff_content : CHARMMFFContent
        An instance of CHARMMFFContent containing the CHARMM force field data.
    streamIDs : list of str
        A list of stream IDs from which residues can be loaded.
    residues : dict
        A dictionary mapping residue names to their corresponding CharmmTopResi objects.
    patches : dict
        A dictionary mapping patch names to their corresponding CharmmTopResi objects.
    masses : CharmmMasses
        An instance of CharmmMasses containing the atom mass records extracted from the CHARMM force field files.
    overrides : dict
        A dictionary containing overrides for specific substreams, such as lipid classification.
    """
    def __init__(self,charmmff_content:CHARMMFFContent,streamIDs=[]):
        logger.debug(f'Initializing CHARMMFFResiDatabase with extra streamIDs {streamIDs}')
        self.charmmff_content=charmmff_content
        self.residues={}
        self.patches={}
        self.masses=CharmmMasses({})
        self.streamIDs=streamIDs
        self.load_from_toplevels()
        logger.debug(f'Loaded {len(self.residues)} residues from toplevels, streams: {self.streamIDs}')
        for streamID in streamIDs:
            logger.debug(f'Loading residues from stream {streamID}')
            self.load_from_stream(streamID)
            if not streamID in self.streamIDs:
                logger.debug(f'Adding stream {streamID} to streams')
                # if the stream is not already in the streams list, add it
                self.streamIDs.append(streamID)
        self.tally_masses()
        self.overrides={
            'substreams':{
                'C6DHPC':'lipid',
                'C7DHPC':'lipid',
                'C8DHPC':'lipid',
                'C9DHPC':'lipid',
                'C10DHPC':'lipid',
                'C11DHPC':'lipid',
                'C12DHPC':'lipid',
                'C13DHPC':'lipid',
                'C14DHPC':'lipid',
                'C15DHPC':'lipid',
                'C16DHPC':'lipid',
                'C17DHPC':'lipid',
                'C18DHPC':'lipid'
            }
        }
        logger.debug(f'Initialized CHARMMFFResiDatabase with extra streamIDs {streamIDs}')

    def load_from_toplevels(self):
        """ 
        Load residues and patches from the toplevel CHARMM force field files.
        This function iterates through the toplevel topology files and extracts residues and patches.
        It updates the `self.residues` and `self.patches` dictionaries with the extracted data.
        """
        for topfile in list(self.charmmff_content.toplevel_top.values())+list(self.charmmff_content.toplevel_toppar.values()):
            new_resis,new_patches=self.load_from_topfile(topfile)
            logger.debug(f'Loaded {len(new_resis)} residues from toplevel {topfile}')
            self.residues.update({x.resname:x for x in new_resis})
            logger.debug(f'Loaded {len(new_patches)} patches from toplevel {topfile}')
            self.patches.update({x.resname:x for x in new_patches})

    def load_from_stream(self,streamID):
        """ 
        Load residues and patches from a specific stream in the CHARMM force field content.
        This function checks if the streamID exists in the CHARMM force field content and loads the residues and patches from the corresponding topology files.

        Parameters
        ----------
        streamID : str
            The ID of the stream from which to load residues and patches.

        Raises
        -------
        Warning
            If the specified streamID is not found in the CHARMM force field content, a warning is logged.
        """
        if streamID not in self.charmmff_content.streams:
            logger.warning(f'load_from_stream: Stream {streamID} not found in CHARMM force field content')
            return
        logger.debug(f'Loading resis from stream {streamID}')
        for topfile in self.charmmff_content.streamfiles[streamID].values():
            this_residues,this_patches=self.load_from_topfile(topfile)
            logger.debug(f'Loaded {len(this_residues)} residues and {len(this_patches)} patches from {topfile} in stream {streamID}')
            self.residues.update({x.resname:x for x in this_residues})
            self.patches.update({x.resname:x for x in this_patches})

    def get_resnames_of_streamID(self,streamID,substreamID=None):
        """ 
        Get the residue names of a specific stream ID.
        This function retrieves the residue names associated with a given stream ID and optional substream ID.

        Parameters
        ----------
        streamID : str
            The ID of the stream from which to retrieve residue names.
        substreamID : str, optional
            The ID of the substream from which to retrieve residue names. If None, all residues in the stream are returned.

        Returns
        -------
        list of str
            A sorted list of residue names associated with the specified stream ID and optional substream ID.

        Raises
        -------
        Warning
            If the specified streamID is not found in the CHARMM force field residue database, a warning is logged and an empty list is returned.
        """
        if streamID not in self.streamIDs:
            logger.warning(f'get_resnames_of_streamID: Stream {streamID} not found in CHARMM force field residue database')
            return []
        resnames=[x.resname for x in self.residues.values() if (x.metadata['streamID']==streamID and ((substreamID is None) or (x.metadata.get('substreamID','')==substreamID)))]
        logger.debug(f'Found {len(resnames)} residues in stream {streamID}')
        resnames.sort()
        return resnames

    def load_from_topfile(self,topfile):
        """ 
        Load residues and patches from a specific topology file.
        This function reads the contents of a topology file, extracts residues and patches, and updates the `self.residues` and `self.patches` dictionaries.


        Parameters
        ----------
        topfile : str
            The name of the topology file to load residues and patches from. This can be a full path or just the basename.

        Returns
        -------
        tuple of list of CharmmTopResi
            A tuple containing two lists: the first list contains the residues extracted from the topology file, and the second list contains the patches.
        """
        logger.debug(f'CHARMMFFResiDatabase::load_from_topfile Loading resis from {topfile}')
        cstr=CHARMMFFStreamID(topfile)
        logger.debug(f'topfile {topfile} CHARMMFFStream: \'{cstr.streamID}\' \'{cstr.substreamID}\'')
        self.masses.update(self.charmmff_content.masses_from_topfile(topfile))
        all_resis=self.charmmff_content.resis_from_topfile(topfile,metadata=dict(
                streamID=cstr.streamID,
                substreamID=cstr.substreamID,
                charmmtopfile=topfile
            ))
        this_residues=[x for x in all_resis if x.key=='RESI']
        this_patches=[x for x in all_resis if x.key=='PRES']
        for resi in this_residues:
            if resi.resname in self.residues:
                logger.debug(f'Residue {resi.resname} found in {topfile} will overwrite already loaded {resi.resname} from {self.residues[resi.resname].metadata["charmmtopfile"]}')
        return this_residues,this_patches

    def tally_masses(self):
        """ 
        Tally the masses of all residues in the database.
        This function iterates through all residues in the database and sets their masses based on the atom mass records.
        It uses the `set_masses` method of each residue to assign the correct masses.
        This is necessary because the atom mass records are stored throughout the topology files, and we need to ensure that all residues have their masses set correctly after all topology files are read in.
        """
        for resi in self.residues.values():
            resi.set_masses(self.masses)

    def add_stream(self,streamID):
        """ 
        Add a stream to the CHARMM force field residue database.
        This function loads residues and patches from the specified stream ID and updates the `self.streamIDs` list.

        Parameters
        ----------
        streamID : str
            The ID of the stream to add to the CHARMM force field residue database.
        """
        self.load_from_stream(streamID)
        if not streamID in self.streamIDs:
            logger.debug(f'Adding stream {streamID} to streams')
            # if the stream is not already in the streams list, add it
            self.streamIDs.append(streamID)
        self.tally_masses()

    def add_topology(self,topfile,streamIDoverride=None):
        """ 
        Add a topology file to the CHARMM force field residue database.
        This function loads residues and patches from the specified topology file and updates the `self.residues` and `self.patches` dictionaries.
        If a `streamIDoverride` is provided, it will set the `streamID` and `substreamID` metadata for all residues and patches loaded from the topology file.

        Parameters
        ----------
        topfile : str
            The name of the topology file to load residues and patches from. This can be a full path or just the basename.
        streamIDoverride : str, optional
            An optional stream ID to override the default stream ID extracted from the topology file. If provided, it will set the `streamID` and `substreamID` metadata for all residues and patches loaded from the topology file.
        """
        new_resis,new_patches=self.load_from_topfile(topfile)
        if streamIDoverride is not None:
            for resi in new_resis+new_patches:
                resi.metadata['streamID']=streamIDoverride
                resi.metadata['substreamID']=''
        self.residues.update({x.resname:x for x in new_resis})
        self.patches.update({x.resname:x for x in new_patches})
        if streamIDoverride is not None and streamIDoverride not in self.streamIDs:
            self.streamIDs.append(streamIDoverride)
        self.tally_masses()

    def get_resi(self,resname):
        """ 
        Get a residue by its name from the CHARMM force field residue database.
        This function checks if the residue name exists in the database and returns the corresponding CharmmTopResi object.

        Parameters
        ----------
        resname : str
            The name of the residue to retrieve from the database.

        Returns
        -------
        CharmmTopResi or None
            The CharmmTopResi object corresponding to the residue name, or None if the residue is not found in the database.
        """
        if resname in self.residues:
            return self.residues[resname]
        else:
            logger.warning(f'Residue {resname} not found in CHARMM force field residue database')
            return None
    
    def __contains__(self,resname):
        """ 
        Check if a residue name exists in the CHARMM force field residue database.
        This function checks if the specified residue name is present in the ``self.residues`` dictionary.

        Parameters
        ----------
        resname : str
            The name of the residue to check for in the database.
            
        Returns
        -------
        bool
            True if the residue name exists in the database, False otherwise.
        """
        return resname in self.residues