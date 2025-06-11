# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PDBInput class for representing PDB files used as inputs for packmol
    Defines the PDBCollection class for managing the collection of said PDBs
"""
import os
import logging
import tarfile
import yaml

logger=logging.getLogger(__name__)

class PDBInput:
    """ A PDBInput object is what is returned by the PDBRepository.checkout() method. It contains """
    def __init__(self,name='',pdbcontents={},info={},opt_tags={}):
        self.name=name
        self.pdbcontents=pdbcontents # conformerID -> pdb contents
        self.info=info # info.yaml contents
        self.opt_tags=opt_tags

    @classmethod
    def from_local_filesystem(cls,name,streamname=''):
        """ Create a PDBInput object for RESI 'name' expected to be found in the local filesystem. """
        openname=os.path.join(streamname,name) if streamname else name
        logger.debug(f'Looking for PDBInput for {name} in {openname}')
        if os.path.exists(f'{openname}.pdb'):
            with open(f'{openname}.pdb','r') as f:
                pdbcontents={0: f.read()}
            return cls(name=name,pdbcontents=pdbcontents)
        if os.path.isdir(openname):
            opt_tags={}
            dircontents=[os.path.basename(x) for x in os.listdir(openname)]
            if not 'info.yaml' in dircontents:
                logger.warning(f'No info.yaml found in directory {openname}')
                return None
            os.chdir(openname)
            with open('info.yaml','r') as f:
                info=yaml.safe_load(f)
            conformerID=0
            while os.path.exists(f'{openname}-{conformerID:02d}.pdb'):
                with open(f'{openname}-{conformerID:02d}.pdb','r') as f:
                    pdbcontents[conformerID]=f.read()
                conformerID+=1
            conformerID=0
            expected_tags=['noh']
            for tag in expected_tags:
                while os.path.exists(f'{openname}-{conformerID:02d}-{tag}.pdb'):
                    with open(f'{openname}-{conformerID:02d}-{tag}.pdb','r') as f:
                        if not tag in opt_tags:
                            opt_tags[tag]={}
                        opt_tags[tag][conformerID]=f.read()
                    conformerID+=1
            os.chdir('..')
            return cls(name=name,pdbcontents=pdbcontents,info=info,opt_tags=opt_tags)

    @classmethod
    def from_tarfile(cls,name,tarfile):
        """ Create a PDBInput object for RESI 'name' expected to be found in the tarfile. """
        if not tarfile:
            return None
        # There is either a single PDB {name}.pdb or a directory {name}
        solo_member=None
        folder_member=None
        for member in tarfile.getmembers():
            if os.path.basename(member.name)==name+'.pdb':
                solo_member=member
                break
            elif member.isdir() and os.path.basename(member.name)==name:
                folder_member=member
                break
        if solo_member is None and folder_member is None:
            logger.warning(f'No members found in tarfile {tarfile.name} for {name}')
            return None
        pdbcontents={}
        info={}
        opt_tags={}
        conformerID=0

        if solo_member:
            # single PDB file case
            with tarfile.extractfile(solo_member) as f:
                pdbcontents[conformerID]=f.read().decode()
        elif folder_member:
            for member in tarfile.getmembers():
                if member.name.startswith(folder_member.name + '/'):
                    # member is inside the folder
                    if member.name.endswith('info.yaml'):
                        # extract info.yaml
                        with tarfile.extractfile(member) as f:
                            info=yaml.safe_load(f)
                    elif member.name.endswith('.pdb'):
                        # extract PDB file
                        tag=''
                        if '-' in member.name:
                            parts = member.name.split('-')
                            if len(parts)>1:
                                token=parts[-1].split('.')[0]
                                if len(token)>1 and token.isdigit():
                                    conformerID = int(token)
                                else:
                                    continue
                                if len(parts) > 2:
                                    tag=parts[-2]
                        with tarfile.extractfile(member) as f:
                            if tag=='':
                                pdbcontents[conformerID]=f.read().decode()
                            else:
                                if not tag in opt_tags:
                                    opt_tags[tag]={}
                                opt_tags[tag][conformerID]=f.read().decode()
        return cls(name=name,pdbcontents=pdbcontents,info=info,opt_tags=opt_tags)
    
    def get_pdb(self,conformerID=0,noh=False):
        if len(self.pdbcontents)==0:
            logger.warning('No PDB contents available to checkout.')
            return None
        if len(self.pdbcontents)==1:
            with open(f'{self.name}.pdb','w') as f:
                f.write(next(iter(self.pdbcontents.values())))
            modified_name=f'{self.name}.pdb'
        else:
            modified_name=f'{self.name}/{self.name}-{conformerID:02d}'
            if noh:
                modified_name += '-noh'
            modified_name += '.pdb'
            with open(os.path.basename(modified_name),'w') as f:
                if not noh:
                    if conformerID in self.pdbcontents:
                        f.write(self.pdbcontents[conformerID])
                    else:
                        logger.warning(f'Conformer ID {conformerID} ({modified_name}) not found in PDB contents.')
                else:
                    if 'noh' in self.opt_tags and conformerID in self.opt_tags['noh']:
                        f.write(self.opt_tags['noh'][conformerID])
                    else:
                        logger.warning(f'Conformer ID {conformerID} ({modified_name}) not found in PDB contents for noh.')
        return os.path.basename(modified_name)
    
    def get_charge(self):
        return self.info.get('charge',0.0)
    
    def get_conformer_data(self,conformerID=0):
        conformers=self.info.get('conformers',[])
        if len(conformers)==0:
            logger.warning('No conformers found in info.yaml.')
            return None
        if conformerID < 0 or conformerID >= len(conformers):
            logger.warning(f'Conformer ID {conformerID} is out of range; returning None')
            return None
        c=conformers[conformerID]
        return c
    
    def get_max_internal_length(self,conformerID=0):
        c=self.get_conformer_data(conformerID)
        if c is None:
            logger.warning('No conformer data available to get max internal length.')
            return 0.0
        return c.get('max-internal-length',0)

    def get_head_tail_length(self,conformerID=0):
        c=self.get_conformer_data(conformerID)
        if c is None:
            logger.warning('No conformer data available to get head-tail length.')
            return 0.0
        return c.get('head-tail-length',0.0)
    
    def get_ref_atoms(self):
        return self.info.get('reference-atoms',{})
    
    def get_parameters(self):
        return self.info.get('parameters',[])

class PDBRepository:
    """
    A PDBRepository is a set of _collections_, each of which respresents a CHARMMFF _stream_.  The base PDBRepository is the one that comes with pestifer, and it is located in PESTIFER/resources/charmmff/pdbrepository/. A user may register additional collections by specifying them in the yaml config file.

    A collection is a _single directory structure_ whose name is a STREAM and whose contents are either:
        - single PDB files with name format RESI.pdb, where RESI is the CHARMM resname
        - directories with the name format RESI, which contain PDB files with the name format RESI-00.pdb, RESI-01.pdb, etc., and a mandatory info.yaml file

    The base PDBRepository contains a 'lipid' collection and a 'water_ions' collection (as of v 1.13.1)

    A collection can be a directory resident in the filesystem or a tarball containing the _one_ such directory. If it is a tarball, the stream name is extracted from the name of the tarball in the format STREAM.tar.gz or STREAM.tgz, and it must be the case that the first directory in every member's name is the same as the stream name.

    A single conceptual stream (e.g., 'prot') can be distributed across multiple collections in a repository. Because all CHARMM resnames are unique, this means that resnames can also be found across multiple collections in a repository.  When a repository is queried for a resname, user-declared collections are searched first in the order they were registered, and then the base collection is searched last.
    """
    _user_count=0
    def __init__(self,basepath=''):
        self.collections={}
        # register the base collection that comes with pestifer
        self.registercollection(basepath,'base')

    def show(self,out_stream=print):
        out_stream('----------------------------------------------------------------------')
        out_stream('PDB Collections:')
        for cname,coll in self.collections.items():
            out_stream(f'\nCollection \'{cname}\' at {coll['path']}:')
            streams=coll.get('streams',[])
            for st,stream in streams.items():
                resis=stream.get('resnames',{})
                for r,info in resis.items():
                    syn=info.get('synonym','')  
                    out_stream(f'{st:>12s}: {r:>12s}{syn:<50s}')
        out_stream('----------------------------------------------------------------------')

    def registercollection(self,pathname,collection_key):
        """registercollection registers a collection into the PDBRepository.
           pathname: the path to the directory or tarball
           collection_key: a key to identify this collection
        """
        this_coll={}
        if not os.path.exists(pathname):
            logger.warning(f'{pathname} does not exist.')
            return
        if not os.path.isdir(pathname):
            logger.warning(f'{pathname} is not a directory.')
            return
        this_coll['path']=pathname
        this_coll['streams']={}
        dircontents=os.listdir(pathname)
        logger.debug(f'Registering collection {collection_key} at {pathname}. Contents: {dircontents}')
        streams_tarred=[x for x in dircontents if any([os.path.splitext(x)[1]==y for y in ['.tar.gz','.tgz']])]
        streams_in_filesystem=[x for x in dircontents if os.path.isdir(os.path.join(pathname,x)) and not x.startswith('.')]
        for streamball in streams_tarred:
            streamname=os.path.splitext(os.path.basename(streamball))[0]
            this_stream={}
            this_stream['resnames']={}
            this_stream['tarfile']=tarfile.open(os.path.join(pathname,streamball),'r:gz')
            members=this_stream['tarfile'].getmembers()
            leads=set([x.name.split('/')[0] for x in members])
            assert len(leads)==1, f'Multiple top-level directories found in tarball {streamball}: {leads}'
            assert streamname==leads.pop(), f'Tarball {streamball} does not match stream name {streamname}.'
            toplevels=[x for x in members if len(x.name.split('/'))==2]
            solos=[os.path.basename(x.name) for x in toplevels if x.name.endswith('.pdb')]
            subdirs=[os.path.basename(x.name) for x in toplevels if x.isdir()]
            for solo in solos:
                info={}
                resname=os.path.splitext(solo)[0]
                this_stream['resnames'][resname]=info
            for subdir in subdirs:
                info={}
                resname=subdir
                info_search=[x for x in members if x.name==f'{streamname}/{subdir}/info.yaml']
                if len(info_search)==0:
                    logger.warning(f'No info.yaml found for {resname} in tarball {streamball}.')
                    continue
                sd_member=info_search[0]
                with this_stream['tarfile'].extractfile(sd_member) as f:
                    info=yaml.safe_load(f)
                this_stream['resnames'][resname]=info
            this_coll['streams'][streamname]=this_stream
        for streamdir in streams_in_filesystem:
            this_stream={}
            streamname=streamdir
            if streamname in this_coll['streams']:
                logger.warning(f'Stream {streamname} already exists in collection {collection_key}. Skipping.')
                continue
            this_stream['directory']=streamdir
            this_stream['resnames']={}
            streampath=os.path.join(pathname,streamdir)
            cwd=os.getcwd()
            os.chdir(streampath)
            dircontents=os.listdir()
            solos=[x for x in dircontents if os.path.isfile(x) and x.endswith('.pdb')]
            subdirs=[x for x in dircontents if os.path.isdir(x) and not x.startswith('.')]
            for solo in solos:
                info={}
                resname=os.path.splitext(solo)[0]
                this_stream['resnames'][resname]=info
            for subdir in subdirs:
                if os.path.exists(os.path.join(subdir,'info.yaml')):
                    with open(os.path.join(subdir,'info.yaml'),'r') as f:
                        info=yaml.safe_load(f)
                    this_stream['resnames'][subdir]=info
            this_coll['streams'][streamname]=this_stream
            os.chdir(cwd)
        self.collections[collection_key]=this_coll

    def __del__(self):
        """ Close all tarballs in the collection """
        for colname in self.collections.keys():
            for stream in self.collections[colname]['streams'].values():
                if 'tarfile' in stream:
                    stream['tarfile'].close()

    def add_usercollection(self,userpath=''):
        self._user_count+=1
        self.registercollection(userpath,f'user{self._user_count}')

    def checkout(self,name):
        """ Given a name, return the PDBInput object for that name, or None if not found.
            If the name is a top-level PDB file, return the PDBInput object for that file.
            If the name is a subdirectory, return the PDBInput object for that subdirectory.
        """
        tar = None
        for c,coll in self.collections.items():
            logger.debug(f'looking for {name} in collection {c}')
            for st,stream in coll['streams'].items():
                logger.debug(f'    looking for {name} in stream {st}')
                if name in stream['resnames']:
                    logger.debug(f'    found {name} in stream {st}')
                    if 'tarfile' in stream:
                        tar=stream['tarfile']
                        return PDBInput.from_tarfile(name,tar)
                    else:
                        return PDBInput.from_local_filesystem(name,streamname=st)
        return None

class PDBRepositoryGenerator: 
    pass