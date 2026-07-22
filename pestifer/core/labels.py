# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Labels and segtypes for residues
This module defines the segment types and residue names used in
CHARMM and PDB files, along with their mappings. It also provides
a class for managing these labels and mappings.
"""
import ast
import json
import logging

from pathlib import Path

logger = logging.getLogger(__name__)

# Path to the *generated* residue-name -> segtype classification.  This file is derived
# from the CHARMM topology file each residue is defined in (see
# :mod:`pestifer.charmmff.segtype_classifier`) and regenerated with
# ``pestifer modify-package charmmff regenerate-segtypes``.  It is loaded and merged into
# :attr:`LabelMappers.segtype_of_resname` at import so that lipids, glycans, protein and
# nucleic-acid residues, CGenFF ligands, etc. are classified without being hand-listed
# below.  The ``_segtypes`` table here holds only the *curated* residue names: those not
# derivable from the force field (PDB-style sugar codes, cross-category exceptions,
# custom residues), the standard protein residues (used to expand psfgen pdbalias
# wildcards), plus per-segtype metadata (``macro`` flag, VMD version gate, rescodes).
_DERIVED_SEGTYPES_PATH = (
    Path(__file__).resolve().parent.parent / 'resources' / 'labels' / 'derived_segtypes.json'
)


def _load_derived_segtypes() -> dict:
    """Load the generated ``segtype -> [resnames]`` classification, or ``{}`` if absent."""
    try:
        with open(_DERIVED_SEGTYPES_PATH) as f:
            return json.load(f).get('segtypes', {})
    except (OSError, ValueError) as e:
        logger.warning(f'could not load derived segtypes from {_DERIVED_SEGTYPES_PATH}: {e}; '
                       'residue classification will be limited to the curated set')
        return {}


_segtypes = {
    'protein': {
        'macro': False,
        'resnames': [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'HSP', 'HSD', 'HSE', 'ILE', 'LEU', 'LYS', 'MET',
            'MSE', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
        'rescodes': {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
            'HIS': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F',
            'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
            'TYR': 'Y', 'VAL': 'V', 'HSP': 'H'},
        'invrescodes': {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HSD', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}},
    # --- curated resnames below: only names NOT derivable from the force field files
    # (PDB-style sugar codes, cross-category exceptions, custom residues).  CHARMM RESI
    # names for lipids/glycans/protein/nucleic/ligands come from the generated
    # derived_segtypes.json (merged at import); do not hand-add those here.
    'ion': {
        'macro': True,
        'resnames': ['CAD', 'CL', 'FE', 'H2PO', 'ZN']},
    'ligand': {
        'macro': True,
        # HEME is defined in a CHARMM prot_heme stream (would derive to 'protein'); pestifer
        # treats it (and the custom EIC/VCG) as a ligand.  ADP/ATP live in a nucleic-acid
        # stream but are used here as ligands.
        'resnames': ['ADP', 'ATP', 'EIC', 'HEM', 'HEME', 'VCG']},
    'cofactor': {
        'macro': True,
        # FMN/NAD/NADP are defined in prot/na streams; COA/PLP/TPP are aliases with no RESI
        'resnames': ['COA', 'FMN', 'NAD', 'NADP', 'PLP', 'TPP']},
    'nucleicacid': {
        'macro': True,
        # single-letter PDB deoxynucleotide codes + FAD (curated as nucleic to match legacy)
        'resnames': ['DA', 'DC', 'DG', 'DT', 'DU', 'FAD']},
    'glycan': {
        'macro': True,
        'vmd_version_lt': '2.0',
        # PDB/IUPAC sugar codes: CHARMM has no RESI of these names (its sugars are
        # alpha/beta-prefixed, e.g. AGLC/BMAN/BGLCNA), so they cannot be derived
        'resnames': [
            '4YS', 'ABE', 'ADA', 'AHR', 'ALL', 'ALT', 'API', 'BAC', 'BEM', 'BGC',
            'BMA', 'COL', 'DIG', 'FRU', 'FUC', 'FUL', 'GAL', 'GCS', 'GCU', 'GL0',
            'GLC', 'GMH', 'GUL', 'GUP', 'IDO', 'IDS', 'KDN', 'KDO', 'LGU', 'LXC',
            'LYX', 'MAL', 'MAN', 'MAV', 'MUR', 'NAG', 'NDG', 'NEU', 'NGA', 'OLI',
            'PAR', 'PSI', 'QUI', 'RAM', 'RIB', 'SGN', 'SIA', 'SOR', 'TAL', 'TYV',
            'WOO', 'XYL', 'XYP', 'XYS']},
    'lipid': {
        'macro': True,
        # lipids whose RESI is not captured by the topfile map (aliases like TOCL->TOCL1,
        # detergents, or names defined in cross-category streams, e.g. retinol BTE2)
        'resnames': [
            'BTE2', 'C6DH', 'C7DH', 'DIPG', 'DIPS', 'DLA', 'DLIPA', 'DLIPC', 'DLIPE',
            'DLIPI', 'DMP', 'MP_1', 'MP_2', 'PMSPE', 'PMSPG', 'SB3_10', 'SB3_12',
            'SB3_14', 'TOCL', 'YOPG']},
    'water': {
        'macro': False,
        # this list also distinguishes water from ions when both are defined in the CHARMM
        # ``water_ions`` stream, so the water model RESI names (TIP3) belong here, not just
        # the PDB aliases (HOH/WAT)
        'resnames': ['HOH', 'TIP3', 'WAT']},
    'other': {
        'macro': False,
        'resnames': ['CAF']},
    }

_atom_aliases = [
    "ILE CD1 CD",
    "MET SE SD",   # selenomethionine's selenium (SE) -> methionine sulfur (SD); pairs with the MSE->MET residue alias
    "BGLCNA C7 C",
    "BGLCNA O7 O",
    "BGLCNA C8 CT",
    "BGLCNA N2 N",
    "ANE5AC C10 C",
    "ANE5AC C11 CT",
    "ANE5AC N5 N",
    "ANE5AC O1A O11",
    "ANE5AC O1B O12",
    "ANE5AC O10 O",
    "VCG C01 C1",
    "VCG C01 C1",
    "VCG C02 C2",
    "VCG C03 C3",
    "VCG C04 C4",
    "VCG C05 C5",
    "VCG C06 C6",
    "VCG C07 C7",
    "VCG C08 C8",
    "VCG C09 C9",
    "TIP3 O OH2",
    # Chloride: the "CL CLA" residue alias below builds the CHARMM CLA residue, whose
    # single atom is named CLA, but RCSB PDBs name the atom CL.  Without this atom alias
    # psfgen cannot match the coordinate and leaves the ion at the origin.  (Zinc needs no
    # such alias: the ZN2 residue's atom is already named ZN, matching the PDB.)
    "CLA CL CLA"
]
_residue_aliases = [
    "HIS HSD",
    "MSE MET",   # selenomethionine -> methionine (SeMet is a crystallography phasing substitution)
    "PO4 H2PO4",
    "H2PO H2PO4",
    "MAN AMAN",
    "BMA BMAN",
    "BGLC BGLCNA",
    "NAG BGLCNA",
    "NDG AGLCNA",
    "FUC AFUC",
    "FUL BFUC",
    "GAL BGAL",
    "SIA ANE5AC",
    "ANE5 ANE5AC",
    "EIC LIN",
    "HOH TIP3",
    "ZN ZN2",
    "CL CLA",
    "C6DH C6DHPC",
    "C7DH C7DHPC",
    "DT THY",
    "DA ADE",
    "DC CYT",
    "DG GUA",
    "DU URA",
    "HEM HEME",
    "TOCL TOCL1",
]

_residue_fullnames = { # these are included because the CHARMMFF topology file does not include them
    "DIPA": "di-alpha-linoleoyl phosphatidic acid",
    "DTPA": "di-hexadecatrienoyl phosphatidic acid",
    "TIPA": "di-palmitoleoyl-2-linoleoyl phosphatidic acid"
}

class LabelMappers:
    """
    Class to hold label mappers for residue names and segment types.
    """
    def __init__(self):
        # process _data to initialize the label mappers
        self.aliases = {}
        self.aliases['atom'] = _atom_aliases
        self.aliases['residue'] = _residue_aliases
        self.residue_fullnames = _residue_fullnames
        self.segtypes = _segtypes
        self.segtype_of_resname = {}
        self.charmm_resname_of_pdb_resname = {}
        self.pdb_resname_of_charmm_resname = {}
        self.res_321 = self.segtypes['protein']['rescodes']
        self.res_123 = self.segtypes['protein']['invrescodes']
        # curated classifications are explicit and win over derived ones
        for segtype in self.segtypes:
            for resname in self.segtypes[segtype]['resnames']:
                self.segtype_of_resname[resname] = segtype
        # merge the force-field-derived classification (loaded from the generated
        # resource); ``setdefault`` keeps curated entries authoritative.  The derived
        # names populate ``segtype_of_resname`` only -- the per-segtype ``resnames``
        # lists in ``_segtypes`` stay curated (they drive psfgen pdbalias wildcards
        # and the atomselect macros).
        for segtype, resnames in _load_derived_segtypes().items():
            for resname in resnames:
                self.segtype_of_resname.setdefault(resname, segtype)
        self.update_alias_mappings()

    def update_alias_mappings(self):
        """
        Update the alias mappings for residues based on the current aliases.
        """
        for alias in self.aliases['residue']:
            parts = alias.split()
            resname, alias1 = parts
            self.charmm_resname_of_pdb_resname[resname] = alias1
            if len(resname) >= 3:
                self.pdb_resname_of_charmm_resname[alias1] = resname
            # PSF files use CHARMM names; propagate segtype so the alias is
            # also recognised by segtype_of_resname (e.g. HEM → HEME).
            if resname in self.segtype_of_resname and alias1 not in self.segtype_of_resname:
                self.segtype_of_resname[alias1] = self.segtype_of_resname[resname]
    
    def update_aliases(self, residue_aliases=[], atom_aliases=[]):
        """
        Update the aliases with new residue and atom aliases.

        Parameters
        ----------
        residue_aliases : list
            A list of strings containing new residue aliases to be added. Each string should be in the format "RESNAME ALIAS".
        atom_aliases : list
            A list of strings containing new atom aliases to be added. Each string should be in the format "RESNAME ATOMNAME ALIAS".
        """
        # logger.debug(f'Updating residue aliases with {residue_aliases}')
        # logger.debug(f'Updating atom aliases with {atom_aliases}')
        self.aliases['residue'].extend(residue_aliases)
        self.aliases['atom'].extend(atom_aliases)
        self.update_alias_mappings()

    def update_segtypes(self, new_segtypes):
        """
        Update the segment types with new resnames.

        Parameters
        ----------
        new_segtypes : dict
            A dictionary containing new segment types to be added or updated. Each key is a segment type, and the value is a list of residue names to *add* to that segment type.
        """
        logger.debug(f'Updating segtypes with {new_segtypes}')
        for segtype, data in new_segtypes.items():
            if segtype not in self.segtypes:
                self.segtypes[segtype] = {}
                self.segtypes[segtype]['resnames'] = data
            else:
                assert 'resnames' in self.segtypes[segtype], f'Segtype {segtype} does not have a "resnames" key.'
                self.segtypes[segtype]['resnames'].extend(data)
            # update segtype_of_resname mapping for this segtype
            for resname in self.segtypes[segtype]['resnames']:
                if resname not in self.segtype_of_resname:
                    self.segtype_of_resname[resname] = segtype

    def update_atomselect_macros(self, fp):
        """
        Update the atomselect macros in the file ``fp`` based on the ``segtypes`` dict.
        This is a developer-only feature.  Access to this method is provided by the ``pestifer modify-package`` command (see :func:`pestifer.cli.pestifer.modify_package`).
        """
        for segtype, data in self.segtypes.items():
            if data['macro']:
                resnames = data['resnames']
                macro_content = ' '.join(resnames)
                vmd_lt = data.get('vmd_version_lt', '')
                if vmd_lt:
                    vmd_major_lt = int(str(vmd_lt).split('.')[0])
                    fp.write(f"if {{[lindex [split [vmdinfo version] .] 0] < {vmd_major_lt}}} {{\n")
                    fp.write(f"    update_atomselect_macro {segtype} \"resname {macro_content}\" 0\n")
                    fp.write(f"}}\n")
                else:
                    fp.write(f"update_atomselect_macro {segtype} \"resname {macro_content}\" 0\n")

def segtype_names() -> list:
    """Return the list of segtype keys defined in :data:`_segtypes`."""
    return list(_segtypes.keys())


def curated_resname_set() -> set:
    """
    Return every residue name that is *curated* (explicitly classified in :data:`_segtypes`).
    The segtype generator excludes these when deriving from the force field, so an explicit
    classification is never overridden by a derived one.
    """
    return {name for seg in _segtypes.values() for name in seg['resnames']}


def water_resnames() -> set:
    """Residue names treated as water (used to split water from ions in ``water_ions``)."""
    return set(_segtypes['water']['resnames'])


def register_resnames_segtype(resnames, segtype: str, labels_path=None):
    """
    Insert one or more residue names into the ``_segtypes[<segtype>]['resnames']``
    list in this module's *source file*, so that a newly-contributed built-in
    custom residue is classified under ``segtype`` on the next pestifer run.

    The edit is made textually against the exact source span of the target list
    (located with :mod:`ast`), preserving the surrounding formatting; names
    already present are skipped.  The edited file is re-parsed to guarantee it is
    still valid Python before it is written back.

    Parameters
    ----------
    resnames : str or list of str
        Residue name(s) to register (upper-cased).
    segtype : str
        Target segtype; must be an existing key of :data:`_segtypes`.
    labels_path : str or Path, optional
        Path to the ``labels.py`` source to edit; defaults to this module's file.

    Returns
    -------
    tuple(Path, list)
        The path of the edited file and the list of names actually added.
    """
    if isinstance(resnames, str):
        resnames = [resnames]
    resnames = [r.strip().upper() for r in resnames if r.strip()]
    if segtype not in _segtypes:
        raise ValueError(f'unknown segtype {segtype!r}; choose from {sorted(_segtypes)}')
    already = set(_segtypes[segtype].get('resnames', []))
    to_add = [r for r in dict.fromkeys(resnames) if r not in already]
    labels_path = Path(labels_path) if labels_path else Path(__file__).resolve()
    if not to_add:
        return labels_path, []

    src = labels_path.read_text()
    tree = ast.parse(src)
    # locate the module-level `_segtypes = { ... }` assignment
    seg_dict = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == '_segtypes' for t in node.targets):
            seg_dict = node.value
            break
    if not isinstance(seg_dict, ast.Dict):
        raise RuntimeError(f'could not locate the _segtypes dict literal in {labels_path}')
    # find the value dict for our segtype, then its 'resnames' list node
    list_node = None
    for key, value in zip(seg_dict.keys, seg_dict.values):
        if isinstance(key, ast.Constant) and key.value == segtype and isinstance(value, ast.Dict):
            for k2, v2 in zip(value.keys, value.values):
                if isinstance(k2, ast.Constant) and k2.value == 'resnames' and isinstance(v2, ast.List):
                    list_node = v2
            break
    if list_node is None:
        raise RuntimeError(f'could not locate _segtypes[{segtype!r}]["resnames"] list in {labels_path}')

    lines = src.splitlines(keepends=True)
    # offset of the closing ']' : end_col_offset points just past it (0-based cols)
    end_line_idx = list_node.end_lineno - 1
    line = lines[end_line_idx]
    bracket_col = list_node.end_col_offset - 1  # index of ']'
    insert_text = ', '.join(repr(n) for n in to_add)
    injection = (', ' + insert_text) if list_node.elts else insert_text
    lines[end_line_idx] = line[:bracket_col] + injection + line[bracket_col:]
    new_src = ''.join(lines)
    ast.parse(new_src)  # guard: the result must still be valid Python
    labels_path.write_text(new_src)
    return labels_path, to_add


Labels = LabelMappers()
"""
Global instance of :class:`LabelMappers` class to access segment types and residue names.
This instance provides access to the segment types and residue names used in CHARMM and PDB files.
It allows for easy mapping between residue names and their corresponding segment types,
as well as providing access to the CHARMM residue names for PDB residue names.  Any module file
that imports ``Labels`` from ``pestifer.core.labels`` will have access to this instance.
This is useful for tasks that require residue name and segment type management, such as
preparing input files for molecular simulations or analyzing protein structures.
"""