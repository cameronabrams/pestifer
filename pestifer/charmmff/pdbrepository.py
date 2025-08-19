# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the PDBInput class for representing PDB files used as inputs for packmol.
Defines the PDBCollection class for managing the collection of said PDBs.
"""
import os
import logging
import yaml

from ..util.cacheable_object import TarBytesFS

logger=logging.getLogger(__name__)

class PDBInfo:
    """ 
    A ``PDBInfo`` object is a container for metadata about a PDB file. It is used to store information about the residue, such as its name, charge, conformers, and other metadata. It can be created from a YAML file or from a tarfile member that is expected to be an info.yaml file.   
    The metadata is stored in a dictionary, and the PDBInfo object provides methods to access the metadata, check equality with another ``PDBInfo`` object, and get items from the metadata dictionary.

    Parameters
    ----------
    metadata : dict
        A dictionary containing metadata about the PDB file.
    """
    def __init__(self,metadata={}):
        self.metadata=metadata

    @classmethod
    def from_yaml_file(cls,yamlfile='info.yaml'):
        metadata={}
        if os.path.exists(yamlfile):
            with open(yamlfile,'r') as f:
                metadata=yaml.safe_load(f)
        else:
            logger.warning(f'{yamlfile} not found; metadata will be empty.')
        return cls(metadata)

    def get(self,key,default=None):
        """ 
        Get a value from the metadata dictionary by key, returning a default value if the key is not found.

        Parameters
        ----------
        key : str
            The key to look up in the metadata dictionary.
            If the key is not found, the default value will be returned.
        default : any, optional
            The value to return if the key is not found.
        """
        return self.metadata.get(key,default)

    def __eq__(self,other):
        """ 
        Check if two PDBInfo objects are equal by comparing their metadata dictionaries.

        Parameters
        ----------
        other : PDBInfo
            The other PDBInfo object to compare against.

        Returns
        -------
        bool
            True if the metadata dictionaries are equal, False otherwise.
        """
        if not isinstance(other,PDBInfo):
            return False
        return self.metadata == other.metadata

    def __getitem__(self,key):
        """ Get an item from the metadata dictionary. """
        return self.metadata.get(key,None)
    def __setitem__(self,key,value):
        """ Set an item in the metadata dictionary. """
        self.metadata[key]=value
    def __contains__(self,key):
        """ Check if a key is in the metadata dictionary. """
        return key in self.metadata
    def __repr__(self):
        """ Return a string representation of the metadata dictionary. """
        return f'PDBInfo({self.metadata})'

class PDBInput:
    """ 
    A ``PDBInput`` object is what is returned by the ``PDBRepository.checkout()`` method. It contains at minimum the data needed to write a PDB file to the local filesystem for use by, e.g., packmol.

    Parameters
    ----------
    name : str
        The name of the residue, which is also the base name of the PDB file.
    pdbcontents : dict
        A dictionary containing the PDB contents for each conformer ID. The keys are conformer IDs and the values are the PDB contents as strings.
    info : PDBInfo
        The metadata for the residue.
    opt_tags : dict
        A dictionary containing optional tags for the residue.
    """
    def __init__(self,name='',pdbcontents={},info:PDBInfo=None,opt_tags={}):
        self.name=name
        self.pdbcontents = pdbcontents  # conformerID -> pdb contents
        self.info = info  # info.yaml contents (the metadata for the residue)
        self.opt_tags = opt_tags

    def get_pdb(self, conformerID=0, noh=False):
        """
        Get the PDB contents for a given conformer ID. If ``noh`` is True, it will return the PDB contents for the ``noh`` tag if available. 
        
        Parameters
        ----------
        conformerID : int
            The conformer ID to retrieve the PDB contents for. Defaults to 0.
        noh : bool
            If True, it will return the PDB contents for the ``noh`` tag if available. Defaults to False.

        Returns
        -------
        str | None
            The PDB contents for the specified conformer ID and tag, or None if not found.
        """
        if len(self.pdbcontents) == 0:
            logger.warning('No PDB contents available to checkout.')
            return None
        if len(self.pdbcontents) == 1:
            with open(f'{self.name}.pdb', 'w') as f:
                f.write(next(iter(self.pdbcontents.values())))
            modified_name = f'{self.name}.pdb'
        else:
            modified_name = f'{self.name}/{self.name}-{conformerID:02d}'
            if noh:
                modified_name += '-noh'
            modified_name += '.pdb'
            with open(os.path.basename(modified_name), 'w') as f:
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
        """
        Get the charge of the residue from the metadata. If not found, return 0.0.
        
        Returns
        -------
        float
            The charge of the residue, or 0.0 if not found."""
        return self.info.get('charge',0.0)
    
    def get_conformer_data(self,conformerID=0):
        """
        Get the conformer data for a given conformer ID. If no conformers are found, return None. If the conformer ID is out of range, return None.

        Parameters
        ----------
        conformerID : int
            The conformer ID to retrieve the data for. Defaults to 0.

        Returns
        -------
        dict | None
            The conformer data for the specified conformer ID, or None if not found.
        """
        conformers=self.info.get('conformers',[])
        if len(conformers)==0:
            logger.warning('No conformers found in metadata.')
            return None
        if conformerID < 0 or conformerID >= len(conformers):
            logger.warning(f'Conformer ID {conformerID} is out of range; returning None')
            return None
        c=conformers[conformerID]
        return c
    
    def get_max_internal_length(self,conformerID=0):
        """
        Get the maximum internal length for a given conformer ID. If no conformers are found, return 0.0. If the conformer ID is out of range, return 0.0.

        Parameters
        ----------
        conformerID : int
            The conformer ID to retrieve the maximum internal length for. Defaults to 0.
        
        Returns
        -------
        float
            The maximum internal length for the specified conformer ID, or 0.0 if not found.
        """
        c=self.get_conformer_data(conformerID)
        if c is None:
            logger.warning('No conformer data available to provide max internal length.')
            return 0.0
        return c.get('max-internal-length',0)

    def get_head_tail_length(self,conformerID=0):
        """
        Get the head-tail length for a given conformer ID. If no conformers are found, return 0.0. If the conformer ID is out of range, return 0.0.

        Parameters
        ----------
        conformerID : int
            The conformer ID to retrieve the head-tail length for. Defaults to 0.

        Returns
        -------
        float
            The head-tail length for the specified conformer ID, or 0.0 if not found.
        """
        c=self.get_conformer_data(conformerID)
        if c is None:
            logger.warning('No conformer data available to provide head-tail length.')
            return 0.0
        return c.get('head-tail-length',0.0)
    
    def get_ref_atoms(self):
        """
        Get the reference atoms for the residue from the metadata. If not found, return an empty dictionary.

        Returns
        -------
        dict
            The reference atoms for the residue, or an empty dictionary if not found.
        """
        return self.info.get('reference-atoms',{})
    
    def get_parameters(self):
        """
        Get the parameters for the residue from the metadata. If not found, return an empty list.
        
        Returns
        -------
        list
            The parameters for the residue, or an empty list if not found.
        """
        return self.info.get('parameters',[])

    def longname(self):
        """
        Get the long name of the residue from the metadata. If not found, return the resname.

        Returns
        -------
        str
            The long name of the residue, or the resname if not found.
        """
        return self.info.get('synonym','').strip() or self.name
        

class PDBCollection:
    """ 
    A ``PDBCollection`` is a collection of ``PDBInput`` objects, each of which represents a residue in the CHARMM force field.  It is used to manage the residues in a single stream, and it can be used to query for residues by name. 

    A collection can be initialized from a directory resident in the filesystem or a tarball containing _one_ such directory. If it is a tarball, the stream name is extracted from the name of the tarball in the format STREAM.tar.gz or STREAM.tgz, and it must be the case that the first directory in every member's name is the same as the stream name (i.e., the tarball must contain a single top-level directory that is the stream name).  The tarball may contain subdirectories with metadata and multiple conformers, but they must not be nested more than one level deep.  The tarball may also contain PDB files with the name format RESI.pdb, which will be treated as solo entries with no metadata.

    A single conceptual stream (e.g., ``prot``) can be distributed across multiple collections in a repository. Because all CHARMM resnames are unique, this means that resnames can also be found across multiple collections in a repository.  When a repository is queried for a resname, user-declared collections are searched first in the order they were registered, and then the base collection is searched last.
    
    Parameters
    ----------
    path : str
        The path to the collection directory or tarball relative to the current working directory. If it ends with `.tar.gz` or `.tgz`, it is treated as a tarball; otherwise, it is treated as a directory.
    streamID_override : str, optional
        If provided, this overrides the stream ID extracted from the path. This is useful if the path does not follow the expected naming convention or if you want to specify a different stream ID for the collection. 
        If not provided, the stream ID is derived from the path or the name of the tarball.
        The stream ID is used to identify the collection and is typically the name of the stream (e.g., ``prot``, ``lipid``, etc.) that the collection represents. 
        It is used to ensure that the collection is correctly associated with the residues it contains. 
    
    """

    def __init__(self, path: str = '', streamID_override: str = ''):
        self.info = {} # resname -> PDBInfo object
        self.contents = {} # resname -> PDBInput object
        self.tarfile = None
        self.path = path # path to the collection directory or tarball relative to CWD
        self.registration_place = 0
        if path.endswith('.tar.gz') or path.endswith('.tgz'):
            basename = os.path.basename(path)
            streamID = os.path.splitext(basename)[0]
            if streamID_override:
                streamID = streamID_override
            self.tarfile = TarBytesFS.from_file(path, compression='gzip')
            root_dir = self.tarfile.ls()[0].get('name', None)
            if not root_dir or root_dir != streamID:
                raise ValueError(f'Tarball {path}\'s top-level directory does not match stream name {streamID} (root_dir {root_dir}).')
            toplevel_solos = [x['name'] for x in self.tarfile.ls(root_dir) if x['type'] != 'directory']
            toplevel_subdirs = [x['name'] for x in self.tarfile.ls(root_dir) if x['type'] == 'directory']
            for solo in toplevel_solos:
                resname = os.path.splitext(os.path.basename(solo))[0]
                self.info[resname] = PDBInfo({})
                with self.tarfile.open(solo, 'r') as f:
                    self.contents[resname] = PDBInput(name=resname, pdbcontents={'0': f.read()}, info=self.info[resname])
            for subdir in toplevel_subdirs:
                resname = os.path.basename(subdir)
                info_search = [x for x in self.tarfile.ls(subdir) if 'info.yaml' in x['name']]
                if len(info_search) == 0:
                    logger.warning(f'No info.yaml found for {resname} in tarball {path}.')
                    continue
                sd_member = info_search[0]
                with self.tarfile.open(sd_member['name'], 'r') as f:
                    contents = yaml.safe_load(f)
                self.info[resname] = PDBInfo(contents)
                pdbcontents = {}
                for index, conformer in enumerate(self.info[resname].metadata['conformers']):
                    with self.tarfile.open(os.path.join(subdir, conformer['pdb']), 'r') as f:
                        pdbcontents[index] = f.read()
                self.contents[resname] = PDBInput(name=resname, pdbcontents=pdbcontents, info=self.info[resname])
        elif os.path.isdir(path):
            logger.debug(f'Scanning directory {path} as a PDBCollection')
            streamID = os.path.basename(path) if not streamID_override else streamID_override
            cwd = os.getcwd()
            os.chdir(path)
            dircontents = os.listdir()
            solos = [x for x in dircontents if os.path.isfile(x) and x.endswith('.pdb')]
            subdirs = [x for x in dircontents if os.path.isdir(x) and not x.startswith('.')]
            for solo in solos:
                resname = os.path.splitext(solo)[0]
                self.info[resname] = PDBInfo({})
                with open(solo, 'r') as f:
                    self.contents[resname] = PDBInput(name=resname, pdbcontents={0: f.read()}, info=self.info[resname])
            pdbcontents = {}
            for subdir in subdirs:
                resname = subdir
                if os.path.exists(os.path.join(subdir, 'info.yaml')):
                    self.info[resname] = PDBInfo.from_yaml_file(os.path.join(subdir, 'info.yaml'))
                for index, conformer in enumerate(self.info[resname].metadata['conformers']):
                    with open(os.path.join(subdir, conformer['pdb']), 'r') as f:
                        pdbcontents[index] = f.read()
            self.contents[resname] = PDBInput(name=resname, pdbcontents=pdbcontents, info=self.info[resname])
            os.chdir(cwd)
        self.streamID = streamID

    def show(self,fullnames=False,missing_fullnames={}):
        """
        Return a string representation of the PDBCollection.

        Parameters
        ----------
        fullnames : bool, optional
            If True, show the full names of the residues (synonyms) instead of just the resnames. Defaults to False.
        missing_fullnames : dict, optional
            A dictionary mapping resnames to their full names if they are not found in the collection.
            This is used to provide full names for residues that do not have a synonym in the metadata.
            Defaults to an empty dictionary.

        Returns
        -------
        str
            A string representation of the PDBCollection, including the stream ID, path, and a list of resnames.
            If `fullnames` is True, it will include the full names (synonyms) of the residues; otherwise, it will just list the resnames.

        """
        retstr=f'PDBCollection(registered_at={self.registration_place}, streamID={self.streamID}, path={self.path}, contains {len(self.info)} resnames)\n'
        if fullnames:
            for resname in sorted(self.info.keys()):
                fullname=self.info[resname].get('synonym','')
                if not fullname:
                    fullname=missing_fullnames.get(resname,'')
                retstr += f'  {fullname:>66s} ({resname})\n'
        else:
            for i,resname in enumerate(sorted(self.info.keys())):
                retstr += f'{resname:>6s}'
                if i<len(self.info)-1:
                    retstr+= ', '
                if (i+1) % 7==0:
                    retstr += '\n'
        return retstr

    def __contains__(self,resname):
        """ 
        Check if a resname is in the collection. 
        
        Parameters
        ----------
        resname : str
            The resname to check for in the collection.

        Returns
        -------
        bool
            True if the resname is found in the collection, False otherwise.
        """
        return resname in self.info

    def checkout(self,resname):
        """ 
        Given a resname, return the ``PDBInput`` object for that resname, or None if not found. 
        
        Parameters
        ----------
        resname : str
            The resname to check out from the collection.  

        Returns
        -------
        PDBInput | None
            The ``PDBInput`` object for the specified resname if found, or None if not found.   
        
        """
        if resname in self.contents:
            logger.debug(f'Found {resname} in collection {self.streamID}')
            return self.contents[resname]
        else:
            logger.debug(f'{resname} not found in collection {self.streamID}')
            return None

class PDBRepository:
    """
    A ``PDBRepository`` is a set of _collections_, each of which respresents a CHARMMFF _stream_.  The base ``PDBRepository`` is the one that comes with pestifer, and it is located in ``PESTIFER/resources/charmmff/pdbrepository/``. The base ``PDBRepository`` contains a ``lipid`` collection and a ``water_ions`` collection (as of v 1.13.1).  A user may register additional collections by specifying them in the yaml config file. 
    """
    def __init__(self):
        self.collections={} # streamID -> PDBCollection object
        self.registration_order=[]

    def add_path(self,path,streamID_override=''):
        """
        Add a new ``PDBCollection`` to the repository from a file path.

        Parameters
        ----------
        path : str
            The file path to the PDB collection.
        streamID_override : str
            An optional override for the stream ID.

        """
        c=PDBCollection(path,streamID_override=streamID_override)
        self.add_collection(c,collection_key=c.streamID)

    def add_collection(self,collection:PDBCollection,collection_key='generic'):
        """ 
        Add a ``PDBCollection`` to the repository. If the collection_key already exists, it will be overwritten. 
        
        Parameters
        ----------
        collection : PDBCollection
            The PDBCollection object to add to the repository.
        collection_key : str
            The key under which to register the collection in the repository. If it already exists, a warning will be logged and the collection will not be added again. If a collection with the same base name already exists, a numbered suffix will be added to the collection_key to avoid conflicts.
        
        """
        if not isinstance(collection,PDBCollection):
            raise TypeError('collection must be a PDBCollection object')
        if collection_key in self.registration_order:
            logger.warning(f'Collection {collection_key} already registered; will not add again.')
            tag=1
            while f'{collection_key}_{tag}' in self.registration_order:
                tag += 1
                if tag>10:
                    raise ValueError(f'Too many collections with the same base name {collection_key}; please choose a different name.')
            collection_key = f'{collection_key}_{tag}'
        self.collections[collection_key]=collection
        self.registration_order.append(collection_key)
        self.collections[collection_key].registration_place=len(self.registration_order)
        logger.debug(f'Added collection {collection_key} with {len(collection.info)} residues.')

    def show(self,out_stream=print,fullnames=False,missing_fullnames={}):
        """ 
        Show the contents of the PDBRepository, including all registered collections and their contents.

        Parameters
        ----------
        out_stream : callable, optional
            A callable that takes a string and outputs it. Defaults to print.
        fullnames : bool, optional
            If True, show the full names of the residues (synonyms) instead of just the resnames. Defaults to False.
        missing_fullnames : dict, optional
            A dictionary mapping resnames to their full names if they are not found in the collection.
            This is used to provide full names for residues that do not have a synonym in the metadata.
            Defaults to an empty dictionary."""
        out_stream('-'*75)
        out_stream('PDB Collections:')
        for cname in self.registration_order[::-1]:
            coll=self.collections[cname]
            out_stream(coll.show(fullnames=fullnames,missing_fullnames=missing_fullnames))
        out_stream('-'*75)

    def __contains__(self, resname: str):
        """
        Check if a resname is in any of the collections.

        Parameters
        ----------
        resname : str
            The resname to check for in the collections.

        """
        for c in self.registration_order[::-1]:
            coll: PDBCollection = self.collections[c]
            if resname in coll.info:
                return True
        return False

    def checkout(self, name: str) -> PDBInput | None:
        """
        Given a name, return the PDBInput object for that name, or None if not found.
        Search is conducted over collections in the order they were registered.

        Parameters
        ----------
        name : str
            The name of the residue to check out from the ``PDBRepository``.
        """
        for c in self.registration_order[::-1]:
            coll: PDBCollection = self.collections[c]
            result = coll.checkout(name)
            if result is not None:
                logger.debug(f'Found {name} in collection {c}')
                return result
        return None
