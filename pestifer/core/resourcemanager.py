# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`ResourceManager` class for managing access to :mod:`pestifer.resources`.
"""

import logging
import os
import shutil

from pathlib import Path
from typing import Callable

from .errors import PestiferError
from .examplemanager import ExampleManager
from .labels import Labels

from .. import resources

from ..charmmff.charmmffcontent import CHARMMFFContent, charmmff_version_key
from ..util.gitutil import get_git_origin_url

logger = logging.getLogger(__name__)


USER_TOPPAR_DIR = Path("~/.pestifer/toppar").expanduser()
""" Well-known per-user cache directory for custom CHARMM stream files. When
this directory exists, it is implicitly prepended to ``charmmff.user_custom.searchpath``
so that a build picks up any ``.str``/``.rtf``/``.prm`` files the user has
dropped there without needing to declare it in YAML. Project-local searchpath
entries override anything in the cache (last-loaded wins). """


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
                raise PestiferError(f'Resource {r} not found at {self.resource_path[r]} -- your installation is likely incomplete')

        self._charmmff_config = charmmff_config
        self._charmmff_content = None
        self._example_manager = None
        self._resname_index_cache = None
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
                raise PestiferError(f'CHARMMFF version directory for {version_str!r} not found at {path}')
            return path
        version_dirs = self.charmmff_version_dirs()
        if not version_dirs:
            raise PestiferError(f'No CHARMMFF version directories found in {self.resource_path["charmmff"]}')
        return version_dirs[-1]  # newest by mtime

    @property
    def charmmff_content(self):
        if self._charmmff_content is None:
            charmmff_path = self.charmmff_version_path(self._charmmff_config.get('release', ''))
            logger.info(f'Using CHARMMFF version: {charmmff_path.name}')
            user_custom = self._charmmff_config.get('user_custom', {})
            user_custom_directories = list(user_custom.get('searchpath', []))
            if USER_TOPPAR_DIR.is_dir():
                # Load the cache first so that project-local searchpath entries
                # (listed in YAML) can override anything that lives there.
                user_custom_directories.insert(0, str(USER_TOPPAR_DIR))
                logger.debug(f'Including user toppar cache: {USER_TOPPAR_DIR}')
            self._charmmff_content = CHARMMFFContent(
                charmmff_path,
                user_custom_directories=user_custom_directories,
                user_custom_segtypes=user_custom.get('segtypes', {}),
                user_pdbrepository_paths=self._charmmff_config.get('pdbrepository', []),
                generate_missing_coordinates=self._charmmff_config.get('generate_missing_coordinates', True),
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
            elif c == 'tcl':
                path = self.get_tcldir()
                with open(os.path.join(path, '00PESTIFER-README.txt'), 'r') as f:
                    msg = f.read()
                out_stream(msg)

    def lookup_resname(self, resname: str) -> dict:
        """Look up a residue name in the CHARMM topologies/streams and in the built-in PDB
        repository.

        Returns a dict describing whether (and where) the residue is defined and whether
        coordinates for it exist:

        - ``resname``: the (upper-cased) name queried
        - ``in_topology``: whether it is defined as a CHARMM ``RESI`` or ``PRES``
        - ``kind``: ``'residue (RESI)'``, ``'patch (PRES)'``, or ``None``
        - ``is_patch``: True if it is a pure patch (the PDB repository is not searched)
        - ``topfile``: the topology/stream file that defines it (or ``None``)
        - ``source``: where the definition comes from -- ``'standard'`` (native to the CHARMM
          release), ``'custom'`` (a pestifer built-in custom file), ``'user'`` (a user-custom
          directory such as ``~/.pestifer/toppar``), or ``None``
        - ``segtype``: pestifer's segtype classification (or ``None``)
        - ``charmm_alias``: the CHARMM resname a PDB resname maps to, if different
        - ``in_pdbrepository``: whether the built-in PDB repository has coordinates
        - ``longname``: the descriptive name from the PDB repository (or ``None``)
        - ``nconformers``: number of stored conformers (0 if none)
        """
        cc = self.charmmff_content
        if cc.pdbrepository is None:
            cc.provision_pdbrepository()  # fast; also picks up any user PDB collections
        query = resname.upper()
        index = self._resname_index()
        entry = index.get(query)
        # user-custom topology residues are read cheaply when CHARMMFFContent is constructed
        # (they are not in the cached built-in index); treat them as residues
        user_custom = query in getattr(cc, 'user_custom_resnames', set())
        in_residues = (entry is not None and entry['kind'] == 'RESI') or (entry is None and user_custom)
        in_patches = entry is not None and entry['kind'] == 'PRES'
        alias = self.labels.charmm_resname_of_pdb_resname.get(query)
        # a pure patch (PRES) modifies a residue and has no coordinates of its own, so the
        # PDB repository is not applicable and is not searched
        is_patch = in_patches and not in_residues
        # provenance: a user-custom directory (~/.pestifer/toppar or a configured searchpath)
        # takes precedence; otherwise the cached index says 'standard' (native to the release)
        # or 'custom' (a pestifer built-in custom file)
        if user_custom:
            source = 'user'
        elif entry:
            source = entry.get('source')
        else:
            source = None
        pdbi = None
        if not is_patch and cc.pdbrepository and query in cc.pdbrepository:
            pdbi = cc.pdbrepository.checkout(query)
        head_tail = pdbi.get_head_tail_length(0) if pdbi else 0.0
        return {
            'resname': query,
            'in_topology': in_residues or in_patches,
            'kind': 'residue (RESI)' if in_residues else ('patch (PRES)' if in_patches else None),
            'is_patch': is_patch,
            'topfile': (entry['topfile'] if entry else None) or cc.get_topfile_of_resname(query),
            'source': source,
            'segtype': self.labels.segtype_of_resname.get(query),
            'charmm_alias': alias if (alias and alias != query) else None,
            'in_pdbrepository': pdbi is not None,
            'longname': pdbi.longname() if pdbi else None,
            'nconformers': len(pdbi.pdbcontents) if pdbi else 0,
            'charge': pdbi.get_charge() if pdbi else None,
            'head_tail_length': head_tail if head_tail else None,
        }

    def _resname_index(self) -> dict:
        """The cached name -> {kind, topfile} index of the built-in CHARMM residues/patches
        (see :class:`~pestifer.charmmff.charmmffcontent.ResnameIndex`).  Built once (loading
        the full residue collection), then cached, so residue lookups avoid that cost."""
        if self._resname_index_cache is None:
            from ..charmmff.charmmffcontent import ResnameIndex
            self._resname_index_cache = ResnameIndex(self.charmmff_content.charmmff_path).index
        return self._resname_index_cache

    def all_resnames(self) -> set:
        """Every residue/patch name pestifer knows about -- the union of the CHARMM
        ``RESI``/``PRES`` names (from the cached index, plus any user-custom residues) and the
        resnames in the PDB repository."""
        cc = self.charmmff_content
        if cc.pdbrepository is None:
            cc.provision_pdbrepository()
        names = set(self._resname_index().keys())
        names |= getattr(cc, 'user_custom_resnames', set())
        if cc.pdbrepository:
            for coll in cc.pdbrepository.collections.values():
                names.update(coll.info.keys())
        return names

    def search_resnames(self, substring: str) -> list:
        """Return the sorted residue names (see :meth:`all_resnames`) that contain
        ``substring`` (case-insensitive)."""
        needle = substring.upper()
        return sorted(name for name in self.all_resnames() if needle in name.upper())

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

    def add_custom_residue(self, source_file, segtype: str = 'ligand', force: bool = False) -> dict:
        """
        Install a CHARMM topology/stream file as a pestifer built-in *custom*
        residue definition.

        The file is copied into the active force field's ``custom/`` directory, the
        ``RESI`` names it defines are registered under ``segtype`` in
        :mod:`pestifer.core.labels` (so psfgen classifies them correctly), and the
        on-disk resource cache is cleared so the new residue is picked up on the
        next run.

        Parameters
        ----------
        source_file : str or Path
            Path to a ``.str``/``.rtf``/``.top``/``.prm`` file containing at least
            one ``RESI`` block.
        segtype : str, optional
            Segtype to classify the residue(s) under (default ``'ligand'``).
        force : bool, optional
            Overwrite an existing custom file and permit residue-name collisions
            with definitions already in the force field.

        Returns
        -------
        dict
            Summary with keys ``resnames``, ``patch_names``, ``segtype``,
            ``segtype_added``, ``dest``, ``labels_path``, and ``touched_paths``
            (the repository files changed, for staging into a commit).
        """
        from ..charmmff.charmmffcontent import extract_resi_pres_blocks, CHARMMFFContent
        from .labels import register_resnames_segtype, segtype_names

        src = Path(source_file).expanduser().resolve()
        if not src.is_file():
            raise PestiferError(f'{src}: no such file')
        if CHARMMFFContent.charmmff_filetype(src.name) is None:
            raise PestiferError(f'{src.name}: not a CHARMM force field file (expected .str, .rtf, .top, or .prm)')
        if segtype not in segtype_names():
            raise PestiferError(f'unknown segtype {segtype!r}; choose from {sorted(segtype_names())}')

        text = src.read_text()
        resnames = [b.split()[1].upper() for b in extract_resi_pres_blocks(text, keywords=('RESI',)) if len(b.split()) > 1]
        patch_names = [b.split()[1].upper() for b in extract_resi_pres_blocks(text, keywords=('PRES',)) if len(b.split()) > 1]
        if not resnames and not patch_names:
            raise PestiferError(f'{src.name}: no RESI or PRES blocks found; nothing to add')
        if not resnames:
            raise PestiferError(f'{src.name}: defines only PRES patch(es) {patch_names}; a custom RESI residue is required')

        if not force:
            clashes = sorted({n for n in resnames + patch_names if self.lookup_resname(n)['in_topology']})
            if clashes:
                raise PestiferError(f'residue name(s) already defined in the force field: {", ".join(clashes)}; use force=True to override')

        customdir = Path(self.get_charmmff_customdir())
        customdir.mkdir(parents=True, exist_ok=True)
        dest = customdir / src.name
        if dest.exists() and not force:
            raise PestiferError(f'{dest} already exists; use force=True to overwrite')
        shutil.copyfile(src, dest)

        labels_path, added = register_resnames_segtype(resnames, segtype)
        touched = [str(dest)] + ([str(labels_path)] if added else [])

        # the cached CHARMMFF content and resname index no longer reflect the new
        # custom file; clear them so the next run rebuilds
        from ..util.cacheable_object import CacheableObject
        CacheableObject.clear_cache()
        self._resname_index_cache = None

        return {
            'resnames': resnames,
            'patch_names': patch_names,
            'segtype': segtype,
            'segtype_added': added,
            'dest': str(dest),
            'labels_path': str(labels_path),
            'touched_paths': touched,
        }

    def _validate_pdb_entry(self, entry_dir):
        """
        Validate a single ``make-pdb-collection`` entry directory and return
        ``(entry_path, resname, info)``.  A ``kind: box`` entry must ship its psf/pdb; a
        molecule entry must list conformers whose PDBs all exist.
        """
        import yaml

        entry = Path(entry_dir).expanduser().resolve()
        if not entry.is_dir():
            raise PestiferError(f'{entry}: not a directory (expected a make-pdb-collection entry dir named after the residue)')
        resname = entry.name
        info_path = entry / 'info.yaml'
        if not info_path.is_file():
            raise PestiferError(f'{entry}: no info.yaml; is this a make-pdb-collection entry directory?')
        with open(info_path) as f:
            info = yaml.safe_load(f) or {}
        if info.get('kind') == 'box':
            # a pre-equilibrated solvent box (make-pdb-collection solvent): psf + pdb, no conformers
            for key in ('psf', 'pdb'):
                fname = info.get(key)
                if not fname:
                    raise PestiferError(f'{info_path}: kind:box entry is missing its {key!r} file reference')
                if not (entry / fname).is_file():
                    raise PestiferError(f'{info_path} references {key} {fname!r}, but {entry / fname} is missing')
        else:
            conformers = info.get('conformers', [])
            if not conformers:
                raise PestiferError(f'{info_path}: no conformers listed')
            for c in conformers:
                if not (entry / c['pdb']).is_file():
                    raise PestiferError(f'{info_path} references conformer {c["pdb"]!r}, but {entry / c["pdb"]} is missing')
        return entry, resname, info

    def _resolve_pdb_collection(self, resname, collection):
        """Return the target collection for ``resname`` (explicit ``collection`` wins)."""
        if collection:
            return collection
        seg = self.labels.segtype_of_resname.get(resname) or self.labels.segtype_of_resname.get(resname.upper())
        coll = {'ion': 'solvent', 'water': 'solvent'}.get(seg, seg)
        if not coll:
            raise PestiferError(f'cannot infer a collection for {resname} (its segtype is unknown); pass collection= explicitly')
        return coll

    def _install_pdb_entries(self, collection, entries, force):
        """
        Add ``entries`` (a list of ``(entry_path, resname, info)``) into the single
        ``<collection>.tgz`` tarball in one extract/repack pass, creating the tarball if
        absent.  Returns ``(tarball_path, created)``.
        """
        import tarfile
        import tempfile

        pdbrepo_dir = Path(self.charmmff_content.charmmff_path) / 'pdbrepository'
        pdbrepo_dir.mkdir(parents=True, exist_ok=True)
        tarball = pdbrepo_dir / f'{collection}.tgz'
        created = not tarball.exists()

        # extract the existing collection (if any) into a staging tree, add/replace each
        # <resname>/ entry, and re-tar with the single top-level directory <collection>
        with tempfile.TemporaryDirectory() as td:
            if tarball.exists():
                with tarfile.open(tarball, 'r:gz') as tf:
                    tf.extractall(td, filter='data')
            staging = Path(td) / collection
            staging.mkdir(parents=True, exist_ok=True)
            for entry_path, resname, _info in entries:
                dest = staging / resname
                if dest.exists():
                    if not force:
                        raise PestiferError(f'{resname} already present in collection {collection!r}; use force=True to overwrite')
                    shutil.rmtree(dest)
                shutil.copytree(entry_path, dest)
            with tarfile.open(tarball, 'w:gz') as tf:
                tf.add(staging, arcname=collection)
        return tarball, created

    def _clear_pdb_caches(self):
        """Invalidate the cached PDB repository / resname index after a repo mutation."""
        from ..util.cacheable_object import CacheableObject
        CacheableObject.clear_cache()
        self._resname_index_cache = None

    def add_pdb_entry(self, entry_dir, collection: str = None, force: bool = False) -> dict:
        """
        Install a single PDB-repository entry -- a residue's sampled coordinate set produced
        by ``make-pdb-collection`` -- into the built-in PDB repository, so that
        ``make_membrane_system`` can place the residue (or ``solvate`` can tile a box).

        A PDB repository is a set of *collections*, each a ``<stream>.tgz`` tarball whose
        single top-level directory is the stream name and which holds one ``<RESI>/``
        subdirectory per residue (its ``info.yaml`` plus conformer PDBs, or a ``kind: box``
        entry's psf/pdb).  This method adds (or, with ``force``, replaces) the ``<RESI>/``
        entry inside the target collection tarball, creating the tarball if it does not exist
        yet.  To install a whole directory of entries at once, use :meth:`add_pdb_collection`.

        Parameters
        ----------
        entry_dir : str or Path
            A directory named after the residue (its basename is the resname), containing
            ``info.yaml`` and the coordinate files it references.
        collection : str, optional
            Collection/stream to install into; the entry lands in
            ``pdbrepository/<collection>.tgz`` under ``<collection>/<RESI>/``.  Defaults to
            the residue's segtype (``ion``/``water`` map to ``solvent``).
        force : bool, optional
            Overwrite an entry already present for this resname in the collection.

        Returns
        -------
        dict
            Summary with keys ``resname``, ``collection``, ``tarball``, ``kind``,
            ``nconformers``, ``created_collection``, and ``touched_paths``.
        """
        entry, resname, info = self._validate_pdb_entry(entry_dir)
        collection = self._resolve_pdb_collection(resname, collection)
        tarball, created = self._install_pdb_entries(collection, [(entry, resname, info)], force)
        self._clear_pdb_caches()
        return {
            'resname': resname,
            'collection': collection,
            'tarball': str(tarball),
            'kind': info.get('kind', 'molecule'),
            'nconformers': len(info.get('conformers', [])),
            'created_collection': created,
            'touched_paths': [str(tarball)],
        }

    def add_pdb_collection(self, entries_dir, collection: str = None, force: bool = False) -> dict:
        """
        Install one or more PDB-repository entries from ``entries_dir`` into the built-in
        repository.  ``entries_dir`` may be either a single entry directory (one that
        contains ``info.yaml``) or a *container* directory holding one ``<RESI>/``
        subdirectory per entry -- e.g. the ``<streamID>/`` tree that
        ``make-pdb-collection --streamID <stream>`` writes.  This is the batch counterpart
        of :meth:`add_pdb_entry`: it validates every entry first, groups them by target
        collection, and repacks each affected tarball once.

        Parameters
        ----------
        entries_dir : str or Path
            A single entry directory, or a directory of entry subdirectories.
        collection : str, optional
            Force every entry into this collection; otherwise each entry goes to its
            residue's default collection (its segtype, with ``ion``/``water`` -> ``solvent``).
        force : bool, optional
            Overwrite entries already present for a resname in the target collection.

        Returns
        -------
        dict
            Summary with keys ``entries`` (list of ``(resname, collection, kind)``),
            ``collections`` (sorted list), ``created_collections``, ``tarballs``, and
            ``touched_paths``.
        """
        root = Path(entries_dir).expanduser().resolve()
        if not root.is_dir():
            raise PestiferError(f'{root}: not a directory')
        if (root / 'info.yaml').is_file():
            entry_dirs = [root]
        else:
            entry_dirs = [d for d in sorted(root.iterdir()) if d.is_dir() and (d / 'info.yaml').is_file()]
            if not entry_dirs:
                raise PestiferError(f'{root}: contains no make-pdb-collection entries (no <RESI>/info.yaml subdirectories) and is not itself an entry directory')

        # validate all entries up front, then group by resolved collection
        groups: dict = {}
        entries_meta = []
        for d in entry_dirs:
            entry, resname, info = self._validate_pdb_entry(d)
            coll = self._resolve_pdb_collection(resname, collection)
            groups.setdefault(coll, []).append((entry, resname, info))
            entries_meta.append((resname, coll, info.get('kind', 'molecule')))

        tarballs = []
        created_collections = []
        for coll in sorted(groups):
            tarball, created = self._install_pdb_entries(coll, groups[coll], force)
            tarballs.append(str(tarball))
            if created:
                created_collections.append(coll)
        self._clear_pdb_caches()
        return {
            'entries': entries_meta,
            'collections': sorted(groups),
            'created_collections': created_collections,
            'tarballs': tarballs,
            'touched_paths': tarballs,
        }

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
            f.write('## Instead, edit :mod:`pestifer.core.labels` and run ``pestifer modify-package charmmff update-atomselect-macros.``\n')
            f.write('## The Tcl proc ``update_atomselect_macro`` is defined in the PestiferUtil Tcl package.\n\n')
            self.labels.update_atomselect_macros(f)

    def regenerate_derived_segtypes(self, out_path=None) -> dict:
        """
        Regenerate the derived residue-name -> segtype classification from the installed
        CHARMM force field(s) and write it to the resource loaded by
        :mod:`pestifer.core.labels`.

        The classification is derived from the topology file each residue is defined in
        (see :mod:`pestifer.charmmff.segtype_classifier`), unioned across every installed
        CHARMMFF release for robustness, and excludes the curated names in
        :data:`pestifer.core.labels._segtypes` (so explicit classifications win).  Run this
        after updating the bundled force field or adding a residue whose classification
        should follow from its defining file.

        Returns
        -------
        dict
            The derived ``segtype -> sorted[resnames]`` mapping that was written.
        """
        import json
        from ..charmmff.segtype_classifier import derive_segtypes
        from .labels import curated_resname_set, water_resnames, _DERIVED_SEGTYPES_PATH

        merged: dict[str, str] = {}
        for vd in self.charmmff_version_dirs():
            for resname, topfile in CHARMMFFContent(vd).resi_to_topfile_map.items():
                merged.setdefault(resname, topfile)
        derived = derive_segtypes(merged, curated_names=curated_resname_set(),
                                  water_resnames=water_resnames())
        out_path = Path(out_path) if out_path else _DERIVED_SEGTYPES_PATH
        out_path.parent.mkdir(parents=True, exist_ok=True)
        doc = {'_comment': 'GENERATED by pestifer modify-package charmmff regenerate-segtypes; do not edit by hand',
               'segtypes': derived}
        with open(out_path, 'w') as f:
            json.dump(doc, f, indent=0, sort_keys=True)
        logger.info(f'Wrote {sum(len(v) for v in derived.values())} derived segtype classifications to {out_path}')
        return derived

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