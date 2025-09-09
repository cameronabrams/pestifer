# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the PDBInput class for representing PDB files used as inputs for packmol.
Defines the PDBCollection class for managing the collection of said PDBs.
"""
from collections import UserDict
import os
import logging
from pathlib import Path
from typing import Callable
import yaml
from dataclasses import dataclass, field
from ..util.cacheable_object import TarBytesFS, CacheableObject
from ..util.util import countTime
from ..util.spinner_wrapper import with_spinner
from ..util.stringthings import my_logger, plu 

logger = logging.getLogger(__name__)

@dataclass
class PDBInput:
    """ 
    A ``PDBInput`` object represents the data needed to use a PDB file as input for packmol.
    """
    name: str = ''
    """ The name of the residue, which is also the base name of the PDB file. """
    pdbcontents: dict = field(default_factory=dict)
    """ A dictionary containing the PDB contents for each conformer ID. The keys are conformer IDs and the values are the PDB contents as strings. """
    info: dict = field(default_factory=dict)
    """ The metadata for the residue. """
    opt_tags: dict = field(default_factory=dict)
    """ A dictionary containing optional tags for the residue. """

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
        return self.info.get('charge', 0.0)

    def get_conformer_data(self, conformerID: int = 0):
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
        conformers = self.info.get('conformers', [])
        if len(conformers) == 0:
            logger.warning('No conformers found in metadata.')
            return None
        if conformerID < 0 or conformerID >= len(conformers):
            logger.warning(f'Conformer ID {conformerID} is out of range; returning None')
            return None
        c = conformers[conformerID]
        return c

    def get_max_internal_length(self, conformerID: int = 0):
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
        c = self.get_conformer_data(conformerID)
        if c is None:
            logger.warning('No conformer data available to provide max internal length.')
            return 0.0
        return c.get('max-internal-length', 0)

    def get_head_tail_length(self, conformerID: int = 0):
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
        c = self.get_conformer_data(conformerID)
        if c is None:
            logger.warning('No conformer data available to provide head-tail length.')
            return 0.0
        return c.get('head-tail-length', 0.0)

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

class PDBInputDict(UserDict[str, PDBInput]):
    """ A dictionary mapping residue names to their corresponding PDBInput objects. """
    pass

@dataclass
class PDBCollection:
    """ 
    A ``PDBCollection`` object is a collection of ``PDBInput`` objects.

    Collection data can be initialized from a directory resident in the filesystem or a tarball containing _one_ such directory. If it is a tarball, the stream name is extracted from the name of the tarball in the format STREAM.tar.gz or STREAM.tgz, and it must be the case that the first directory in every member's name is the same as the stream name (i.e., the tarball must contain a single top-level directory that is the stream name).  The tarball may contain subdirectories with metadata and multiple conformers, but they must not be nested more than one level deep.  The tarball may also contain PDB files with the name format RESI.pdb, which will be treated as solo entries with no metadata.

    A single conceptual stream (e.g., ``prot``) can be distributed across multiple collections in a repository. Because all CHARMM resnames are unique, this means that resnames can also be found across multiple collections in a repository.  When a repository is queried for a resname, user-declared collections are searched first in the order they were registered, and then the base collection is searched last.
    """
    path_or_tarball: str | Path = ''
    """ The path to the collection directory or tarball relative to the current working directory. If it ends with `.tar.gz` or `.tgz`, it is treated as a tarball; otherwise, it is treated as a directory. """
    streamID: str = ''
    """ Name of the CHARMMFF stream corresponding to this collection. """
    info: dict = field(default_factory=dict)
    """ Contents of each entry's info.yaml file (charge, conformerIDs, synonyms, etc.) """
    contents: PDBInputDict = field(default_factory=PDBInputDict)
    """ Mapping of residue names to their corresponding PDBInput objects. """
    registration_place: int = 0
    """ The registration place of the collection in the repository indicating when it registered. """

    @classmethod
    def build_from_resources(cls, path_or_tarball: str, resnames: list[str] = [], streamID_override: str = None):
        if not path_or_tarball:
            return cls()
        streamID = os.path.splitext(os.path.basename(path_or_tarball))[0] if not streamID_override else streamID_override
        # logger.debug(f'{path_or_tarball} streamID {streamID}')
        contents = PDBInputDict({})
        info = {}
        if path_or_tarball.endswith('.tar.gz') or path_or_tarball.endswith('.tgz'):
            logger.debug(f'Initializing PDBCollection from tarball {path_or_tarball}')
            pdbrepo_fs = TarBytesFS.from_file(path_or_tarball, compression='gzip')
            root_dir = pdbrepo_fs.ls()[0].get('name', None)
            logger.debug(f'Tarball {path_or_tarball} root directory: {root_dir}')
            if not root_dir or root_dir != streamID:
                raise ValueError(f'Tarball {path_or_tarball}\'s top-level directory does not match stream name {streamID} (root_dir {root_dir}).')
            toplevel_solos = [x['name'] for x in pdbrepo_fs.ls(root_dir) if x['type'] != 'directory']
            toplevel_subdirs = [x['name'] for x in pdbrepo_fs.ls(root_dir) if x['type'] == 'directory']
            for solo in toplevel_solos:
                # update info and contents for a solo PDB file (info is empty in this case)
                resname = os.path.splitext(os.path.basename(solo))[0]
                if not resnames or resname in resnames:
                    info[resname] = {}
                    with pdbrepo_fs.open(solo, 'r') as f:
                        contents[resname] = PDBInput(name=resname, pdbcontents={'0': f.read()}, info=info[resname])
            for subdir in toplevel_subdirs:
                resname = os.path.basename(subdir)
                if not resnames or resname in resnames:
                    info_search = [x for x in pdbrepo_fs.ls(subdir) if 'info.yaml' in x['name']]
                    if len(info_search) == 0:
                        logger.warning(f'No info.yaml found for {resname} in tarball {path_or_tarball}.')
                        continue
                    sd_member = info_search[0]
                    with pdbrepo_fs.open(sd_member['name'], 'r') as f:
                        info[resname] = yaml.safe_load(f)
                    pdbcontents = {}
                    for index, conformer in enumerate(info[resname]['conformers']):
                        with pdbrepo_fs.open(os.path.join(subdir, conformer['pdb']), 'r') as f:
                            pdbcontents[index] = f.read()
                    contents[resname] = PDBInput(name=resname, pdbcontents=pdbcontents, info=info[resname])
            if len(contents) == 0:
                logger.debug(f'No valid PDB contents for {"any resnames" if not resnames else resnames} found in tarball {path_or_tarball}.')
            del pdbrepo_fs
            return cls(path_or_tarball=path_or_tarball, streamID=streamID, info=info, contents=contents)
        elif os.path.isdir(path_or_tarball):
            logger.debug(f'Initializing PDBCollection from {path_or_tarball}')
            cwd = os.getcwd()
            os.chdir(path_or_tarball)
            dircontents = os.listdir()
            solos = [x for x in dircontents if os.path.isfile(x) and x.endswith('.pdb')]
            subdirs = [x for x in dircontents if os.path.isdir(x) and not x.startswith('.')]
            for solo in solos:
                resname = os.path.splitext(solo)[0]
                if not resnames or resname in resnames:
                    info[resname] = {}
                    with open(solo, 'r') as f:
                        contents[resname] = PDBInput(name=resname, pdbcontents={0: f.read()}, info=info[resname])
            pdbcontents = {}
            for subdir in subdirs:
                resname = subdir
                if not resnames or resname in resnames:
                    if os.path.exists(os.path.join(subdir, 'info.yaml')):
                        with open(os.path.join(subdir, 'info.yaml'), 'r') as f:
                            info[resname] = yaml.safe_load(f)
                    for index, conformer in enumerate(info[resname]['conformers']):
                        with open(os.path.join(subdir, conformer['pdb']), 'r') as f:
                            pdbcontents[index] = f.read()
                    contents[resname] = PDBInput(name=resname, pdbcontents=pdbcontents, info=info[resname])
            os.chdir(cwd)
            if len(contents) == 0:
                logger.debug(f'No valid PDB contents for {"any resnames" if not resnames else resnames} found in directory {path_or_tarball}.')
            return cls(path_or_tarball=path_or_tarball, streamID=streamID, info=info, contents=contents)

    def show(self, fullnames: bool = False, missing_fullnames: dict = None) -> str:
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
        retstr = f'PDBCollection(registered_at={self.registration_place}'
        retstr += f', streamID=\'{self.streamID}\''
        retstr += f', path=\'{self.path_or_tarball}\''
        retstr += f', contains {len(self.info)} resnames)\n'
        if fullnames:
            pass
            for resname in sorted(self.info.keys()):
                fullname = self.info[resname].get('synonym', '')
                if not fullname:
                    fullname = missing_fullnames.get(resname, '')
                retstr += f'  {fullname:>66s} ({resname})\n'
        else:
            pass
            for i, resname in enumerate(sorted(self.info.keys())):
                retstr += f'{resname:>6s}'
                if i < len(self.info) - 1:
                    retstr += ', '
                if (i + 1) % 7 == 0:
                    retstr += '\n'
        return retstr

    def __contains__(self, resname: str) -> bool:
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

    def checkout(self, resname: str) -> PDBInput | None:
        """ 
        Given a resname, return the ``PDBInput`` object for that resname, or None if not found. 
        
        Parameters
        ----------
        resname : str
            The resname to check out from the collection.  

        Returns
        -------
        PDBInputManager | None
            The ``PDBInputManager`` object for the specified resname if found, or None if not found.   

        """
        if resname in self.contents:
            # logger.debug(f'Found {resname} in collection {self.streamID}')
            return self.contents[resname]
        else:
            # logger.debug(f'{resname} not found in collection {self.streamID}')
            return None

class PDBCollectionDict(UserDict[str, PDBCollection]):
    pass

class PDBRepository(CacheableObject):
    """
    A ``PDBRepository`` is a set of _collections_, each of which respresents a CHARMMFF _stream_.  The base ``PDBRepository`` is the one that comes with pestifer, and it is located in ``PESTIFER/resources/charmmff/pdbrepository/``. The base ``PDBRepository`` contains a ``lipid`` collection and a ``water_ions`` collection (as of v 1.13.1).  A user may register additional collections by specifying them in the yaml config file. 
    """
    
    @countTime
    def __init__(self, *args, **kwargs):
        is_custom = 'resnames' in kwargs and len(kwargs['resnames']) > 0
        if is_custom:
            self.build_custom(*args, **kwargs)
        else:
            super().__init__(*args, **kwargs)

    @with_spinner('No cache yet -- building PDBRepository from package resources...')
    def _build_from_resources(self, charmmff_pdbrepository_path: str = '', **kwargs):
        """
        Build a full collection that represents the complete resource set.
        """
        self.collections: PDBCollectionDict = PDBCollectionDict({})
        self.registration_order: list[str] = []
        members = os.listdir(charmmff_pdbrepository_path)
        for m in members:
            logger.debug(f'Adding {m} to PDBRepository from {charmmff_pdbrepository_path}')
            datapath = os.path.join(charmmff_pdbrepository_path, m)
            self.add_resource(datapath)

    def build_custom(self, charmmff_pdbrepository_path: str = '', streamID_override: str = '', resnames: list[str] = [], **kwargs):
        """
        Build a custom collection that represents a user-defined set of residues.

        Parameters
        ----------
        charmmff_pdbrepository_path : str
            The path to the charmmff/pdbrepository directory.
        streamID_override : str
            An optional override for the stream ID.
        resnames: list[str]
            An optional list of resnames that the collection should minimally include.  This is for building custom, right-sized collections for any particular build.
        """
        self.collections: PDBCollectionDict = PDBCollectionDict({})
        self.registration_order: list[str] = []
        members = os.listdir(charmmff_pdbrepository_path)
        for m in members:
            logger.debug(f'Adding {m} to custom PDBRepository from {charmmff_pdbrepository_path}')
            datapath = os.path.join(charmmff_pdbrepository_path, m)
            self.add_resource(datapath, streamID_override=streamID_override, resnames=resnames)

    def add_resource(self, path_or_tarball: str = '', streamID_override: str = '', resnames: list[str] = []):
        """
        Add a new ``PDBCollection`` to the repository from a file path.

        Parameters
        ----------
        path : str
            The file path to the PDB collection.
        streamID_override : str
            An optional override for the stream ID.
        resnames: list[str]
            An optional list of resnames that the collection should minimally include.  This is for building custom, right-sized collections for any particular build.
        """
        c = PDBCollection.build_from_resources(path_or_tarball=path_or_tarball, streamID_override=streamID_override, resnames=resnames)
        # logger.debug(c.streamID)
        self.add_collection(c, collection_key=c.streamID)

    def add_collection(self, collection: PDBCollection, collection_key='generic'):
        """ 
        Add a ``PDBCollection`` to the repository. If the collection_key already exists, it will be overwritten. 
        
        Parameters
        ----------
        collection : PDBCollection
            The PDBCollection object to add to the repository.
        collection_key : str
            The key under which to register the collection in the repository. If it already exists, a warning will be logged and the collection will not be added again. If a collection with the same base name already exists, a numbered suffix will be added to the collection_key to avoid conflicts.
        
        """
        if not isinstance(collection, PDBCollection):
            raise TypeError('collection must be a PDBCollection object')
        # logger.debug(f'registration_order \'{self.registration_order}\'')
        if collection_key in self.registration_order:
            logger.warning(f'Collection {collection_key} already registered; will not add again.')
            tag=1
            while f'{collection_key}_{tag}' in self.registration_order:
                tag += 1
                if tag > 10:
                    raise ValueError(f'Too many collections with the same base name {collection_key}; please choose a different name.')
            collection_key = f'{collection_key}_{tag}'
        self.collections[collection_key] = collection
        self.registration_order.append(collection_key)
        # logger.debug(f' -> registration_order \'{self.registration_order}\' streamID {collection.streamID}')
        self.collections[collection_key].registration_place = len(self.registration_order)
        logger.debug(f'Added collection {collection_key} with {len(collection.info)} residue{plu(len(collection.info))}.')

    def show(self, out_stream: Callable = print, fullnames: bool = False, missing_fullnames: dict = {}):
        """ 
        Show the contents of the PDBRepository, including all registered collections and their contents.

        Parameters
        ----------
        out_stream : callable, optional
            A callable that takes a string and outputs it. Defaults to # logger.debug.
        fullnames : bool, optional
            If True, show the full names of the residues (synonyms) instead of just the resnames. Defaults to False.
        missing_fullnames : dict, optional
            A dictionary mapping resnames to their full names if they are not found in the collection.
            This is used to provide full names for residues that do not have a synonym in the metadata.
            Defaults to an empty dictionary.
        """
        out_stream('-'*75)
        out_stream('PDB Collections:')
        for cname in self.registration_order[::-1]:
            coll: PDBCollection = self.collections[cname]
            out_stream(coll.show(fullnames=fullnames, missing_fullnames=missing_fullnames))
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
        Given a name, return the PDBInputManager object for that name, or None if not found.
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
