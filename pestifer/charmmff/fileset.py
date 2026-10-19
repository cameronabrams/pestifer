# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""One ordered list of CHARMM force-field files, and the two views a build needs.

A build loads a *set* of CHARMM-format files, and the order it loads them in decides
two things that nothing downstream reports:

* **Topology load order fixes terminal patches.**  psfgen gives each residue the
  ``DEFA FIRS ... LAST ...`` in force when that residue's definition is read, not when
  the segment is built.  ``toppar_all36_lipid_prot.str`` has its ``DEFA`` line commented
  out, so ``GLYM`` inherits whatever the previously loaded file left behind.
* **Parameter load order decides which duplicate wins.**  :meth:`CharmmParamFile.merge`
  is last-wins, matching CHARMM's ``READ PARAM APPEND``.

The second is not hypothetical here.  Merging the shipped default set in its real order
versus "all .prm, then all .str" changes exactly one record out of 46,527 --

    ('NG2O1', 'CG2R61', 'CG2R61', 'NG2S3') n=2    Kchi 1.25  ->  3.1

-- which is `toppar_pestifer_conflict_resolutions.prm` deliberately overriding
`toppar_all36_carb_imlab.str`, annotated in the file as *"CGenFF v5.0 value, chosen over
imlab's analogy value"*.  Group the list by extension for tidiness and that arbitration
silently reverts.  So :data:`DEFAULT_FILES` is written in true load order, and both views
below derive from it without reordering.

Before 2026-10-19 this lived in three config groups named for a file's *provenance* --
``standard`` (in the CHARMM release), ``custom`` (shipped by pestifer, not in the
release), ``user_custom`` (the user's own) -- which named identity where the user needed
to express intent.  "Add one more stream from the release" had to go in ``standard``,
whose own docs said it was "not expected to be modified by the user"; "replace the
default set" was not expressible at all; and ``user_custom.{rtf,str,prm}`` were read by
no scripter whatsoever, so naming your own file there validated, built, and loaded
nothing.  See :func:`fold_legacy_groups`, which maps the old keys onto the new ones.
"""
import logging
import os

logger = logging.getLogger(__name__)

#: The shipped load set, **in load order**.  Reordering this is a force-field change:
#: see the module docstring.  Entries are basenames; a file from the CHARMM release is
#: found wherever it lives in the release tree, and a file in a ``charmmff.searchpath``
#: directory is found there.
DEFAULT_FILES = [
    # --- release topologies ---
    'top_all36_prot.rtf',
    'top_all35_ethers.rtf',
    'top_all36_cgenff.rtf',
    'top_all36_lipid.rtf',
    'top_all36_carb.rtf',
    'top_all36_na.rtf',
    # --- release parameters ---
    'par_all36m_prot.prm',
    'par_all36_lipid.prm',
    'par_all36_carb.prm',
    'par_all36_na.prm',
    'par_all36_cgenff.prm',
    # --- release streams (topology *and* parameters) ---
    'toppar_water_ions.str',
    'toppar_all36_carb_glycopeptide.str',
    'toppar_all36_carb_imlab.str',
    'toppar_all36_prot_modify_res.str',
    # --- pestifer's own corrections; last-wins, so they must follow the release ---
    'toppar_pestifer_dihedral_fills.prm',
    'toppar_pestifer_conflict_resolutions.prm',
    'toppar_all36_moreions.str',
]

#: psfgen wants ``top_all35_ethers.rtf`` read last among the topologies.
_TOPOLOGY_LAST = 'top_all35_ethers.rtf'


def _ext(name: str) -> str | None:
    """``'par'``, ``'top'``, ``'toppar'``, or None -- the same classification
    :meth:`CHARMMFFContent.charmmff_filetype` uses, duplicated here only so this module
    imports nothing from the content layer."""
    if name.endswith('.prm'):
        return 'par'
    if name.endswith('.rtf') or name.endswith('.top'):
        return 'top'
    if name.endswith('.str'):
        return 'toppar'
    return None


def resolve(charmmff_config: dict) -> list[str]:
    """Return the ordered list of files this build loads.

    ``files`` is appended to :data:`DEFAULT_FILES` unless ``replace_defaults`` is true,
    in which case ``files`` **is** the set.  A user addition therefore lands last and
    wins any duplicate -- including over pestifer's own corrections, which is the right
    reading of an explicit request.  Duplicates are dropped keeping the *first*
    occurrence, so naming a default file again does not move it.
    """
    cfg = charmmff_config or {}
    user = list(cfg.get('files') or [])
    base = [] if cfg.get('replace_defaults') else list(DEFAULT_FILES)
    out = []
    for name in base + user:
        if _ext(name) is None:
            raise ValueError(
                f'charmmff.files: {name!r} is not a CHARMM force-field file; expected '
                f'one of .rtf, .top, .str, .prm')
        if name not in out:
            out.append(name)
    return out


def topologies(files: list[str]) -> list[str]:
    """The ``topology`` lines psfgen emits: ``.rtf``/``.top`` in order, then ``.str`` in
    order, with :data:`_TOPOLOGY_LAST` moved to the end if present."""
    out = [f for f in files if _ext(f) == 'top'] + [f for f in files if _ext(f) == 'toppar']
    if _TOPOLOGY_LAST in out:
        out.remove(_TOPOLOGY_LAST)
        out.append(_TOPOLOGY_LAST)
    return out


def parameters(files: list[str]) -> list[str]:
    """The ``parameters`` lines NAMD gets: every ``.prm`` and ``.str``, **in list
    order**, because :meth:`CharmmParamFile.merge` is last-wins."""
    return [f for f in files if _ext(f) in ('par', 'toppar')]


#: Config groups replaced on 2026-10-19 by the flat ``files`` list.
_LEGACY_GROUPS = ('standard', 'custom', 'user_custom')


def fold_legacy_groups(charmmff_config: dict) -> dict:
    """Map a pre-2026-10-19 ``charmmff`` block onto the flat shape, in place.

    The old groups could only ever *grow*: ycleptic merges a user list into the schema
    default by appending, so an entry present in a group's default could not be removed
    from a config.  Folding is therefore exactly "append whatever is not already in
    :data:`DEFAULT_FILES`", preserving each group's own order.

    ``user_custom.{rtf,str,prm}`` are folded too, which is the first time naming a file
    there has ever loaded it.
    """
    cfg = charmmff_config
    if not isinstance(cfg, dict):
        return cfg
    extra, seen = [], set(DEFAULT_FILES) | set(cfg.get('files') or [])
    for group in _LEGACY_GROUPS:
        block = cfg.get(group)
        if not isinstance(block, dict):
            continue
        for key in ('rtf', 'str', 'prm'):
            for name in block.get(key) or []:
                if name not in seen:
                    seen.add(name)
                    extra.append((group, key, name))
    if extra:
        logger.warning(
            'charmmff.{%s} are deprecated in favor of the flat charmmff.files list '
            '(and charmmff.searchpath / charmmff.segtypes). Folding %d file(s) into '
            'charmmff.files: %s',
            ','.join(_LEGACY_GROUPS), len(extra),
            ', '.join(f'{n} (was {g}.{k})' for g, k, n in extra))
        cfg['files'] = list(cfg.get('files') or []) + [n for _, _, n in extra]
    legacy_uc = cfg.get('user_custom')
    if isinstance(legacy_uc, dict):
        if legacy_uc.get('searchpath') and not cfg.get('searchpath'):
            cfg['searchpath'] = list(legacy_uc['searchpath'])
        if legacy_uc.get('segtypes') and not cfg.get('segtypes'):
            cfg['segtypes'] = dict(legacy_uc['segtypes'])
    return cfg


def searchpath(charmmff_config: dict) -> list[str]:
    """Directories to index for the user's own CHARMM-format files, with ``~`` and
    ``$VAR`` expanded."""
    raw = (charmmff_config or {}).get('searchpath') or []
    return [os.path.expanduser(os.path.expandvars(str(d))) for d in raw]


def add(charmmff_config: dict, name: str, extension: str | None = None) -> None:
    """Append *name* to ``charmmff.files`` if it is not already in the resolved set.

    For pestifer's own generators (solvent boxes, PDB collections), which discover at
    run time that a residue's defining topology is not in the load set.  *extension* is
    accepted for call-site readability and only guards the file type; the real
    classification is by suffix.
    """
    if extension is not None and extension not in ('rtf', 'top', 'str', 'prm'):
        logger.debug(f'not adding {name}: extension {extension!r} is not a CHARMM file type')
        return
    if _ext(name) is None:
        logger.debug(f'not adding {name}: not a CHARMM force-field file')
        return
    if name in resolve(charmmff_config):
        return
    charmmff_config['files'] = list(charmmff_config.get('files') or []) + [name]
