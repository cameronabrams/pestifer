# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
import os
import re
import tarfile
from .charmmtop import CharmmMassRecord, CharmmMasses, CharmmTopResi
from .pdbrepository import PDBRepository
from .stringthings import my_logger

logger=logging.getLogger(__name__)

import re

def extract_resi_pres_blocks(text, keywords=('RESI', 'PRES')):
    keyword_pattern = '|'.join(re.escape(k) for k in keywords)
    # Include ATOMS and BONDS as additional block-end sentinels
    terminator_pattern = rf'^(?:{keyword_pattern}|ATOMS|BONDS|ANGLES|DIHEDRALS)\b'
    # Match starting with RESI or PRES, ending before the next RESI, PRES, ATOMS, or EOF
    pattern = rf'(?i)^(({keyword_pattern})\s.*?)(?={terminator_pattern}|\Z)'
    matches = re.finditer(pattern, text, flags=re.DOTALL | re.MULTILINE)
    return [m.group(1).strip() for m in matches]

def extract_mass_lines(file_contents):
    return [line for line in file_contents.splitlines() if line.strip().upper().startswith("MASS")]

comment_these_out=['set','if','WRNLEV','BOMLEV']

class CHARMMFFContent:
    """ A class for handling all CHARMM force field content.  The CHARMM force field is
    stored in a tarball downloaded directly from the MacKerell lab at the University of Michigan.
    https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
    """
    def __init__(self,charmmff_path='',tarfilename='toppar_c36_jul24.tgz'):
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
            self.filenamemap[f]=os.path.join(self.charmmff_path,'custom',f)

    def __del__(self):
        """ Close the tarfile if it is open """
        if self.tarfile is not None:
            self.tarfile.close()
            self.tarfile=None
            logger.debug('Closed CHARMM force field tarfile')
        else:
            logger.debug('No CHARMM force field tarfile to close')
        del self.pdbrepository

    def load_charmmff(self,tarfilename='toppar_c36_jul24.tgz',skip_streams=['cphmd','misc']):
        """ Load the entire CHARMM force field tarball """
        def okfilename(name):
            """ Check if a filename is ok to use in the tarfile """
            return not 'history' in name and not 'all22' in name and not 'ljpme' in name and (name.endswith('.str') or name.endswith('.prm') or name.endswith('.rtf'))

        """ Load the CHARMM force field tarball from the specified path """
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
        self.all_topology_files=[x for x in self.filenamemap.values() if x.endswith('.str') or x.endswith('.rtf')]

    def copy_charmmfile_local(self,basename):
        """ Given a basename for any charmmff file, extract from the existing unaltered charmmff installation
            and modify it for use directly in NAMD if necessary"""
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
                lines=file.readlines()
                logger.debug(f'found {len(lines)} lines in {basename} in custom files')
                with open(basename,'w') as f:
                    for l in lines:
                        is_comment=any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l+'\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}\n')
        elif basename in self.filenamemap:
            logger.debug(f'found {basename} in at {self.filenamemap[basename]} in tarball')
            with self.tarfile.extractfile(self.filenamemap[basename]) as file:
                logger.debug(f' file has type {type(file)}')
                lines=file.read().decode().splitlines()
                logger.debug(f'type of lines is {type(lines)}')
                logger.debug(f'found {len(lines)} lines in {basename} in tarfile')
                with open(basename,'w') as f:
                    for l in lines:
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
        """ remove all local charmmff files
        """
        for f in os.listdir('.'):
            if f.startswith('par') or f.startswith('top_') or f.startswith('toppar') or f.startswith('charmm') or f.endswith('.str') or f.endswith('.prm') or f.endswith('.rtf'):
                os.remove(f)
    
    def lines_from_topfile(self,topfile):
        """ Extract the lines from a top file """
        lines=self.contents_from_topfile(topfile).splitlines()
        return lines

    def contents_from_topfile(self,topfile):
        """ Extract the contents from a top file """
        content=''
        # logger.debug(f'Extracting content from {topfile}')
        for m in self.tarmembers:
            if topfile==m.name:
                # logger.debug(f'Found {topfile} in tarfile member {m.name}')
                with self.tarfile.extractfile(m) as f:
                    content=f.read().decode()
                    break
        if not content:
            mapped_name=self.filenamemap.get(topfile,None)
            if mapped_name is not None:
                # logger.debug(f'Extracting lines from {topfile} using mapped name {mapped_name}')
                with open(mapped_name,'r') as f:
                    content=f.read()
        return content
    
    def masses_from_topfile(self,topfile):
        masses=[]
        lines=extract_mass_lines(self.contents_from_topfile(topfile))
        for i in range(len(lines)):
            masses.append(CharmmMassRecord(lines[i]))
        # for m in masses:
        #     logger.debug(f'Found mass record: \'{str(m)}\'')
        return CharmmMasses(masses)

    def resis_from_topfile(self,topfile,metadata={}):
        blocks=extract_resi_pres_blocks(self.contents_from_topfile(topfile))
        R=[]
        for block in blocks:
            # logger.debug(block)
            key=block.split()[0].upper()
            # if key!='RESI':
            #     logger.debug(f'Expected RESI block, but found {key} in {topfile}')
            resi=CharmmTopResi(block,key=key)
            resi.metadata=metadata
            R.append(resi)
        return R

    def get_topfile_of_patchname(self,patchname):
        """ Given a patch name, return the top file that contains it """
        for topfile in self.all_topology_files:
            lines=self.lines_from_topfile(topfile)
            for line in lines:
                if line.upper().startswith('PRES') and patchname in line:
                    logger.debug(f'Found patch {patchname} in {topfile}')
                    return os.path.basename(topfile)
        return None

class CHARMMFFStreamID:
    # patterns: par_all3x_STREAM.prm (x=5,6), top_all3x_STREAM.rtf (x=5,6), toppar_STREAM.str, toppar_all36_STREAM_SUBSTREAM.str
    # SUBSTREAM can also be _-delimited
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
    """ A class for handling the CHARMM force field residue database.  This is a dictionary of residue names
    to their corresponding stream and substream.
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
        for topfile in list(self.charmmff_content.toplevel_top.values())+list(self.charmmff_content.toplevel_toppar.values()):
            new_resis,new_patches=self.load_from_topfile(topfile)
            logger.debug(f'Loaded {len(new_resis)} residues from toplevel {topfile}')
            self.residues.update({x.resname:x for x in new_resis})
            logger.debug(f'Loaded {len(new_patches)} patches from toplevel {topfile}')
            self.patches.update({x.resname:x for x in new_patches})

    def load_from_stream(self,streamID):
        """ Load residues from a specific stream """
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
        """ Get a list of residue names in a specific stream """
        if streamID not in self.streamIDs:
            logger.warning(f'get_resnames_of_streamID: Stream {streamID} not found in CHARMM force field residue database')
            return []
        resnames=[x.resname for x in self.residues.values() if (x.metadata['streamID']==streamID and ((substreamID is None) or (x.metadata.get('substreamID','')==substreamID)))]
        logger.debug(f'Found {len(resnames)} residues in stream {streamID}')
        resnames.sort()
        return resnames

    def load_from_topfile(self,topfile):
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
        # atom mass records are stored throughout the topfiles, so we will only set the masses once all tops are read in
        for resi in self.residues.values():
            resi.set_masses(self.masses)

    def add_stream(self,streamID):
        """ Add a stream to the database """
        self.load_from_stream(streamID)
        if not streamID in self.streamIDs:
            logger.debug(f'Adding stream {streamID} to streams')
            # if the stream is not already in the streams list, add it
            self.streamIDs.append(streamID)
        self.tally_masses()

    def add_topology(self,topfile,streamIDoverride=None):
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
        """ Get a residue by its name """
        if resname in self.residues:
            return self.residues[resname]
        else:
            logger.warning(f'Residue {resname} not found in CHARMM force field residue database')
            return None
    
    def __contains__(self,resname):
        """ Check if a residue is in the database """
        return resname in self.residues