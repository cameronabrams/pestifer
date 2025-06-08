# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PDBInput class for representing PDB files used as inputs for packmol
    Defines the PDBCollection class for managing the collection of said PDBs
"""
import os
import glob
import logging
import tarfile
import yaml

logger=logging.getLogger(__name__)

class PDBInput:
    def __init__(self,name='',pdbcontents={},info={},opt_tags={}):
        self.name=name
        self.pdbcontents=pdbcontents # conformerID -> pdb contents
        self.info=info # info.yaml contents
        self.opt_tags=opt_tags

    @classmethod
    def from_tarfile(cls,name,tarfile):
        """ Create a PDBInput object from a tarfile containing PDB files and an info.yaml file """
        if not tarfile:
            return None
        members=[]
        for member in tarfile.getmembers():
            if name in member.name:
                members.append(member)
        if len(members)==0:
            logger.warning(f'No members found in tarfile {tarfile.name} for {name}')
            return None
        pdbcontents={}
        info={}
        opt_tags={}
        conformerID=0
        for member in members:
            if member.name.endswith('.pdb'):
                # extract conformer ID from file name if possible
                tag=''
                if '-' in member.name:
                    parts = member.name.split('-')
                    if len(parts)>1: # e.g., 'name-optional_tag-01.pdb'
                        token=parts[-1].split('.')[0]
                        if len(token)>1 and token.isdigit():
                            conformerID = int(token)
                        else:
                            continue
                        if len(parts) > 2:  # e.g., 'name-tag-01.pdb'
                            tag=parts[-2]
                with tarfile.extractfile(member) as f:
                    if tag=='':
                        pdbcontents[conformerID]=f.read().decode()
                    else:
                        if not tag in opt_tags:
                            opt_tags[tag]={}
                        opt_tags[tag][conformerID]=f.read().decode()
            elif member.name.endswith('info.yaml'):
                with tarfile.extractfile(member) as f:
                    info=yaml.safe_load(f)
        return cls(name=name,pdbcontents=pdbcontents,info=info,opt_tags=opt_tags)
    
    def get_pdb(self,conformerID=0,noh=False):
        if len(self.pdbcontents)==0:
            logger.warning('No PDB contents available to checkout.')
            return None
        if len(self.pdbcontents)==1:
            with open(f'{self.name}.pdb','w') as f:
                f.write(next(iter(self.pdbcontents.values())))
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

    def get_charge(self):
        return self.info.get('charge',None)
    
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
    """ A PDB collection is defined as a directory containing one or more subdirectories, each of which
      corresponds to a 'stream' of PDB files; e.g., proteins, nucleic acids, lipids, etc., following
      the convention of the CHARMM force field.  The subdirectores can be either untarred directories or
      tarballs.  The subdirectories can contain PDB files or further subdirectories, which can also
      contain PDB files.  The PDB files should be named with the RESI name from the CHARMM force fiel."""
    _user_count=0
    def __init__(self,basepath=''):
        self.collections={}
        # register the base collection that comes with pestifer
        self.registercollection(basepath,'base')

    def show(self,out_stream=print):
        out_stream('----------------------------------------------------------------------')
        out_stream('PDB Collections:')
        infos={}
        for cname,coll in self.collections.items():
            out_stream(f'\n{cname} collection at {coll['path']}:')
            for st,L in coll['streams'].items():
                members=L.get('tarfile',None)
                toplevels=[]
                if tarfile:
                    members=tarfile.getmembers()
                    for m in members:
                        name=m.name.split('/')
                        if len(name)==1 and name[0].endswith('.pdb'):
                            toplevels.append(os.path.splitext(name[0])[0])
                        elif len(name)==2 and name[1]=='info.yaml':
                            toplevels.append(name[0])
                            with tarfile.extractfile(m) as f:
                                infos[name[0]]=yaml.safe_load(f)
                for name in toplevels:
                    syn=infos.get(name,{}).get('synonym','')
                    if syn:
                        syn=f' \'{syn}\''
                    out_stream(f'{st:>12s}: {name:>12s}{syn:<50s}')
        out_stream('----------------------------------------------------------------------')

    def registercollection(self,pathname,pathkey):
        this_coll={}
        if not os.path.exists(pathname):
            logger.warning(f'{pathname} does not exist.')
            return
        if not os.path.isdir(pathname):
            logger.warning(f'{pathname} is not a directory.')
            return
        this_coll['path']=pathname
        this_coll['streams']={}
        streams_tarred=[x for x in glob.glob(os.path.join(pathname,'*.tar.gz'))+glob.glob(os.path.join(pathname,'*.tgz')) if os.path.isfile(x) and (x.endswith('.tar.gz') or x.endswith('.tgz'))]
        if len(streams_tarred)==0:
            logger.warning(f'{pathname} is not a valid pdb collection; no archives found.')
            return
        for streamball in streams_tarred:
            this_stream={}
            this_stream['tarfile']=tarfile.open(streamball,'r:gz')
            this_stream['resnames']=[os.path.splitext(os.path.basename(x.name))[0] for x in this_stream['tarfile'].getmembers()]
            streamname=os.path.splitext(os.path.basename(streamball))[0]
            this_coll['streams'][streamname]=this_stream
        self.collections[pathkey]=this_coll

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
                    tar=stream['tarfile']
                    return PDBInput.from_tarfile(name,tar)
        return None

