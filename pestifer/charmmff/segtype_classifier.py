# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Derive a residue's *segtype* from the CHARMM topology/stream file it is defined in.

The CHARMM force field already organizes residues by category across its topology and
stream files (``*_lipid_*`` for lipids, ``*_carb_*`` for carbohydrates, ``*_prot_*`` for
protein, ``*_na_*`` for nucleic acids, ``*_cgenff*`` for CGenFF ligands, ``water_ions``
for water and ions, ...).  Classifying by the defining file therefore replaces a large
hand-maintained enumeration of residue names: a residue added to a CHARMM lipid stream
is a lipid automatically.

This module is deliberately dependency-free (it takes a plain ``resname -> topfile``
mapping) so it can be exercised in isolation and reused by the persist-time generator.
"""

# Ordered (substring, segtype) rules matched against the lowercased basename of the
# topology file a residue is defined in.  First match wins, so the list runs from the
# most specific/unambiguous file family to the most general.  ``water_ions`` is matched
# before the ion rules and split downstream (water resnames vs everything-else-is-an-ion).
_TOPFILE_RULES = [
    ('cgenff', 'ligand'),        # top_all36_cgenff*.rtf, <resi>-cgenff.str
    ('water_ions', 'water_ions'),  # toppar_water_ions.str -> split into water / ion
    ('moreions', 'ion'),         # toppar_all36_moreions.str
    ('lipid', 'lipid'),          # top_all36_lipid.rtf, toppar_all36_lipid_*.str
    ('carb', 'glycan'),          # top_all36_carb.rtf, toppar_all36_carb_*.str
    ('36_na', 'nucleicacid'),    # top_all36_na.rtf, toppar_all36_na_*.str
    ('prot', 'protein'),         # top_all36_prot.rtf, toppar_all36_prot_*.str
]


def segtype_of_topfile(topfile, water_resnames=frozenset(), resname=None):
    """
    Return the segtype implied by ``topfile`` (a topology/stream file basename), or
    ``None`` if no rule matches (e.g. ``topfile`` is ``None`` because the name is a
    PDB-style alias not defined directly in the force field).

    ``toppar_water_ions.str`` holds both water and ions; when a residue comes from it,
    ``resname`` is classified as ``water`` if it is in ``water_resnames`` and ``ion``
    otherwise.
    """
    if not topfile:
        return None
    t = topfile.lower()
    for substr, segtype in _TOPFILE_RULES:
        if substr in t:
            if segtype == 'water_ions':
                return 'water' if (resname and resname in water_resnames) else 'ion'
            return segtype
    return None


def derive_segtypes(resi_to_topfile_map, curated_names=frozenset(), water_resnames=frozenset()):
    """
    Classify every residue in ``resi_to_topfile_map`` by its defining topology file.

    Parameters
    ----------
    resi_to_topfile_map : dict
        ``resname -> topfile-basename`` (a CHARMMFFContent attribute).
    curated_names : set
        Residue names that are classified explicitly elsewhere (the curated core in
        :mod:`pestifer.core.labels`); these are skipped so an explicit classification is
        never overridden by a derived one.
    water_resnames : set
        Residue names treated as water when they come from ``water_ions``.

    Returns
    -------
    dict
        ``segtype -> sorted list of resnames`` for every residue that a rule classified
        and that is not already curated.
    """
    curated = set(curated_names)
    out: dict[str, list] = {}
    for resname, topfile in resi_to_topfile_map.items():
        if resname in curated:
            continue
        segtype = segtype_of_topfile(topfile, water_resnames=water_resnames, resname=resname)
        if segtype is None:
            continue
        out.setdefault(segtype, []).append(resname)
    return {segtype: sorted(names) for segtype, names in sorted(out.items())}
