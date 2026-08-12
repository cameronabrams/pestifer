# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the PDBInput class for representing PDB files used as inputs for grid-packed membrane building.
Defines the PDBCollection class for managing the collection of said PDBs.
"""
from collections import UserDict

# Suffix marking a residue's liquid-ordered conformer ensemble in a PDB collection, and the
# markers used to fold it onto the residue's own listing entry.
_LO_SUFFIX = '__Lo'
_LO_BOTH_MARK = '*'
_LO_ONLY_MARK = '^'
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
    A ``PDBInput`` object represents the data needed to use a PDB file as input for grid-packed membrane building.
    """
    name: str = ''
    """ The name of the residue, which is also the base name of the PDB file. """
    pdbcontents: dict = field(default_factory=dict)
    """ A dictionary containing the PDB contents for each conformer ID. The keys are conformer IDs and the values are the PDB contents as strings.  For a ``box`` entry there is a single entry holding the box PDB. """
    info: dict = field(default_factory=dict)
    """ The metadata for the residue. """
    opt_tags: dict = field(default_factory=dict)
    """ A dictionary containing optional tags for the residue. """
    psf_content: str = ''
    """ The PSF contents of a pre-equilibrated solvent box (``kind: box`` entries only; empty for ``molecule`` entries). """

    def is_box(self) -> bool:
        """True if this is a pre-equilibrated solvent *box* entry (``kind: box``) rather
        than a single-molecule conformer entry (``kind: molecule``, the default)."""
        return self.info.get('kind', 'molecule') == 'box'

    def get_box_edge(self):
        """The equilibrated cubic edge length (Å) of a box entry -- VMD ``solvate``'s
        ``-ws``. ``None`` for a molecule entry."""
        return self.info.get('box_edge')

    def get_key_atom(self):
        """The atom name occurring once per molecule of a box entry -- VMD ``solvate``'s
        ``-ks``. ``None`` for a molecule entry."""
        return self.info.get('key_atom')

    def get_box_pdb(self):
        """Write the box PDB to ``<name>-box.pdb`` in the cwd and return the filename
        (for VMD ``solvate -spdb``). ``None`` if this is not a box entry."""
        if not self.is_box() or not self.pdbcontents:
            return None
        fn = f'{self.name}-box.pdb'
        with open(fn, 'w') as f:
            f.write(next(iter(self.pdbcontents.values())))
        return fn

    def get_box_psf(self):
        """Write the box PSF to ``<name>-box.psf`` in the cwd and return the filename
        (for VMD ``solvate -spsf``). ``None`` if this is not a box entry."""
        if not self.is_box() or not self.psf_content:
            return None
        fn = f'{self.name}-box.psf'
        with open(fn, 'w') as f:
            f.write(self.psf_content)
        return fn

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

    def get_z_extent(self, conformerID: int = 0) -> float:
        """Coordinate z-span (``max_z - min_z``) of one conformer -- see :meth:`get_max_z_extent`."""
        text = self.pdbcontents.get(conformerID)
        if text is None and conformerID == 0 and self.pdbcontents:
            text = next(iter(self.pdbcontents.values()))   # solo entry keyed '0'/0
        if text is None:
            return 0.0
        zs = [float(ln[46:54]) for ln in text.splitlines() if ln.startswith(('ATOM', 'HETATM'))]
        return (max(zs) - min(zs)) if zs else 0.0

    def get_max_z_extent(self) -> float:
        """The largest coordinate z-span (``max_z - min_z``) over all conformers.

        Conformers are generated oriented head-up along z, so a conformer's z-span is the vertical
        thickness the lipid occupies in a leaflet -- the quantity the grid packer must reserve so
        lipids do not overflow into (and thin) the water chambers.  Unlike
        :meth:`get_head_tail_length` (head to tail *tip*), this does not collapse when a fluid
        conformer curls its tail tip back toward the head; for rod conformers the two agree.  The
        max over conformers covers the packer's per-lipid draw across the ensemble.  Returns 0.0 if
        there are no coordinates.
        """
        best = 0.0
        for text in self.pdbcontents.values():
            zs = [float(ln[46:54]) for ln in text.splitlines()
                  if ln.startswith(('ATOM', 'HETATM'))]
            if zs:
                best = max(best, max(zs) - min(zs))
        return best

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
            # Solo entries are single-PDB collection members; require the .pdb extension so
            # hidden bookkeeping files that ride along in the tarball -- notably the per-entry
            # `.<resname>.lock` advisory locks -- are not mistaken for entries (a `.PSM.lock`
            # would otherwise register a phantom resname `.PSM`, doubling the collection). This
            # mirrors the directory branch below, which already filters solos to `.pdb`.
            toplevel_solos = [x['name'] for x in pdbrepo_fs.ls(root_dir)
                              if x['type'] != 'directory' and x['name'].endswith('.pdb')]
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
                    if info[resname].get('kind', 'molecule') == 'box':
                        # a pre-equilibrated solvent box: one psf + one pdb, no conformers
                        with pdbrepo_fs.open(os.path.join(subdir, info[resname]['pdb']), 'r') as f:
                            box_pdb = f.read()
                        with pdbrepo_fs.open(os.path.join(subdir, info[resname]['psf']), 'r') as f:
                            box_psf = f.read()
                        contents[resname] = PDBInput(name=resname, pdbcontents={0: box_pdb},
                                                     info=info[resname], psf_content=box_psf)
                    else:
                        pdbcontents = {}
                        for index, conformer in enumerate(info[resname].get('conformers', [])):
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
            for subdir in subdirs:
                resname = subdir
                if not resnames or resname in resnames:
                    if os.path.exists(os.path.join(subdir, 'info.yaml')):
                        with open(os.path.join(subdir, 'info.yaml'), 'r') as f:
                            info[resname] = yaml.safe_load(f)
                    if info[resname].get('kind', 'molecule') == 'box':
                        with open(os.path.join(subdir, info[resname]['pdb']), 'r') as f:
                            box_pdb = f.read()
                        with open(os.path.join(subdir, info[resname]['psf']), 'r') as f:
                            box_psf = f.read()
                        contents[resname] = PDBInput(name=resname, pdbcontents={0: box_pdb},
                                                     info=info[resname], psf_content=box_psf)
                    else:
                        pdbcontents = {}
                        for index, conformer in enumerate(info[resname].get('conformers', [])):
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
        # A '<resname>__Lo' key is not a distinct residue: it is the same residue's
        # liquid-ordered conformer ensemble.  Listing both forms side by side nearly doubled the
        # lipid listing (266 entries for 219 residues) and interleaved every name with its own
        # near-duplicate, so the phase is folded into a marker on the residue instead.
        entries = self._display_entries()
        n_both = sum(1 for _, m in entries if m == _LO_BOTH_MARK)
        n_lo_only = sum(1 for _, m in entries if m == _LO_ONLY_MARK)
        counts = f'contains {len(entries)} residues'
        notes = []
        if n_both:
            notes.append(f'{n_both} with an Lo conformer ensemble')
        if n_lo_only:
            notes.append(f'{n_lo_only} Lo ensembles only')
        if notes:
            counts += ' (' + '; '.join(notes) + ')'

        retstr = f'PDBCollection(registered_at={self.registration_place}'
        retstr += f', streamID=\'{self.streamID}\''
        retstr += f', path=\'{self.path_or_tarball}\''
        retstr += f', {counts})\n'
        if n_both or n_lo_only:
            legend = []
            if n_both:
                legend.append(f'{_LO_BOTH_MARK} also has a liquid-ordered (Lo) conformer ensemble')
            if n_lo_only:
                legend.append(f'{_LO_ONLY_MARK} Lo conformer ensemble only (base conformers come from another collection)')
            retstr += '  ' + '; '.join(legend) + '\n'

        if fullnames:
            for resname, mark in entries:
                fullname = self.info.get(resname, {}).get('synonym', '')
                if not fullname:
                    fullname = self.info.get(resname + _LO_SUFFIX, {}).get('synonym', '')
                if not fullname:
                    fullname = missing_fullnames.get(resname, '')
                retstr += f'  {fullname:>66s} ({resname}{mark})\n'
        else:
            width = max((len(r) + len(m) for r, m in entries), default=6)
            for i, (resname, mark) in enumerate(entries):
                retstr += f'{resname + mark:>{width}s}'
                if i < len(entries) - 1:
                    retstr += ', '
                if (i + 1) % 7 == 0:
                    retstr += '\n'
        return retstr

    def _display_entries(self) -> list:
        """``(resname, marker)`` pairs for display, one per residue rather than one per key.

        Keys ending in ``__Lo`` carry a residue's liquid-ordered conformer ensemble, so they are
        folded onto the residue itself.  The two markers are kept distinct because a collection can
        hold an Lo ensemble *without* the base conformers -- the on-demand user cache is exactly
        that -- and marking those the same way would imply the collection provides a residue it
        does not.
        """
        keys = set(self.info.keys())
        lo_keys = {k for k in keys if k.endswith(_LO_SUFFIX)}
        base = keys - lo_keys
        lo_bases = {k[:-len(_LO_SUFFIX)] for k in lo_keys}
        out = []
        for resname in sorted(base | lo_bases):
            if resname in base and resname in lo_bases:
                out.append((resname, _LO_BOTH_MARK))
            elif resname in base:
                out.append((resname, ''))
            else:
                out.append((resname, _LO_ONLY_MARK))
        return out

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
    A ``PDBRepository`` is a set of _collections_, each of which respresents a CHARMMFF _stream_.  The base ``PDBRepository`` is the one that comes with pestifer, and it is located in ``PESTIFER/resources/charmmff/pdbrepository/``. The base ``PDBRepository`` contains a ``lipid`` collection and a ``solvent`` collection (water + ions, formerly ``water_ions``).  A user may register additional collections by specifying them in the yaml config file. 
    """
    
    @countTime
    def __init__(self, *args, **kwargs):
        is_custom = 'resnames' in kwargs and len(kwargs['resnames']) > 0
        if is_custom:
            self.build_custom(*args, **kwargs)
        else:
            if args and 'resource_label' not in kwargs:
                kwargs['resource_label'] = Path(args[0]).parent.name
            super().__init__(*args, **kwargs)

    @with_spinner('Building PDBRepository cache..')
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
        if c is None:
            logger.debug(f'Skipping {path_or_tarball}: not a recognized PDB collection format.')
            return
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
