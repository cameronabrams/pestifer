# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Citations owed by a build, determined from the build's own configuration.

Pestifer is an orchestrator: the science in a prepared system belongs to NAMD, VMD, psfgen and
the CHARMM force field, and -- depending on what the user asked for -- to CGenFF, to the
carbohydrate or lipid parameter sets, or to PDB2PQR and PROPKA.  A user cannot reasonably be
expected to know which of those a given YAML input pulled in, so pestifer works it out and says
so at the end of the run.

Every conditional citation carries the reason it applies, so the report is a derivation rather
than a list to be taken on faith.  Nothing here inspects the built system: the config is the
input the user wrote, and the citations follow from it.

Coordinates are cited too.  A build that starts from a PDB entry, grafts glycans from another,
or fuses a domain from a third owes the experimental papers behind all three, and only the config
knows which they were.  Those citations are resolved from the structure files the build already
downloaded -- no network call is made on their account -- and are reported as DOI and PMID rather
than as reconstructed reference strings, because the PDB format stores author names and titles in
upper case and the original capitalization is not recoverable from them (``A.B.MCDERMOTT`` cannot
be turned back into ``McDermott, A.B.``).  An identifier cannot be subtly wrong in that way.

Sources for the reference strings: the CHARMM/NAMD/VMD/psfgen entries match the bibliography of
the pestifer paper; the PDB2PQR and PROPKA entries are the citations those packages themselves
declare (``pdb2pqr.config.CITATIONS`` and ``propka.__doc__``).
"""

import glob
import logging
import os
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Citation:
    """One reference, what it is owed for, and (when conditional) why it applies here."""
    subject: str
    text: str
    doi: str = ''
    reason: str = ''

    def render(self):
        line = f'[{self.subject}] {self.text}'
        if self.doi:
            line += f' doi:{self.doi}'
        return line


# --- the catalog -----------------------------------------------------------------------------
# Unconditional: every build runs NAMD, drives VMD/psfgen, and uses the CHARMM36 protein set.

NAMD = Citation(
    'NAMD',
    'Phillips JC, et al. Scalable molecular dynamics with NAMD. '
    'J Comput Chem 26, 1781-1802 (2005).', '10.1002/jcc.20289')

NAMD3 = Citation(
    'NAMD 3',
    'Phillips JC, et al. Scalable molecular dynamics on CPU and GPU architectures with NAMD. '
    'J Chem Phys 153, 044130 (2020).', '10.1063/5.0014475',
    reason='this build ran NAMD 3')

VMD = Citation(
    'VMD',
    'Humphrey W, Dalke A, Schulten K. VMD: Visual molecular dynamics. '
    'J Mol Graph 14, 33-38 (1996).', '10.1016/0263-7855(96)00018-5')

PSFGEN = Citation(
    'psfgen',
    "Ribeiro J, Radak B, Stone J, Gullingsrud J, Saam J, Phillips J. "
    "psfgen User's Guide, v2.0 (2020).")

CHARMM36M = Citation(
    'CHARMM36m (protein)',
    'Huang J, et al. CHARMM36m: An improved force field for folded and intrinsically disordered '
    'proteins. Nat Methods 14, 71-73 (2016).', '10.1038/nmeth.4067')

CHARMM_OVERVIEW = Citation(
    'CHARMM force field',
    'Vanommeslaeghe K, MacKerell AD. CHARMM additive and polarizable force fields for biophysics '
    'and computer-aided drug design. Biochim Biophys Acta 1850, 861-871 (2015).',
    '10.1016/j.bbagen.2014.08.004')

CHARMM_LIPID = Citation(
    'CHARMM36 (lipids)',
    'Klauda JB, et al. Update of the CHARMM All-Atom Additive Force Field for Lipids: Validation '
    'on Six Lipid Types. J Phys Chem B 114, 7830-7843 (2010).', '10.1021/jp101759q')

CHARMM_CARB = Citation(
    'CHARMM36 (carbohydrates)',
    'Guvench O, et al. CHARMM Additive All-Atom Force Field for Glycosidic Linkages between '
    'Hexopyranoses. J Chem Theory Comput 5, 2353-2370 (2009).', '10.1021/ct900242e')

CGENFF = Citation(
    'CGenFF',
    'Vanommeslaeghe K, et al. CHARMM general force field: A force field for drug-like molecules '
    'compatible with the CHARMM all-atom additive biological force fields. '
    'J Comput Chem 31, 671-690 (2010).', '10.1002/jcc.21367')

CGENFF_PROGRAM_I = Citation(
    'CGenFF program',
    'Vanommeslaeghe K, MacKerell AD. Automation of the CHARMM General Force Field (CGenFF) I: '
    'Bond Perception and Atom Typing. J Chem Inf Model 52, 3144-3154 (2012).', '10.1021/ci300363c')

CGENFF_PROGRAM_II = Citation(
    'CGenFF program',
    'Vanommeslaeghe K, Raman EP, MacKerell AD. Automation of the CHARMM General Force Field '
    '(CGenFF) II: Assignment of Bonded Parameters and Partial Atomic Charges. '
    'J Chem Inf Model 52, 3155-3168 (2012).', '10.1021/ci3003649')

PDB2PQR = Citation(
    'PDB2PQR',
    'Dolinsky TJ, et al. PDB2PQR: expanding and upgrading automated preparation of biomolecular '
    'structures for molecular simulations. Nucleic Acids Res 35, W522-W525 (2007).')

APBS = Citation(
    'APBS/PDB2PQR',
    'Jurrus E, et al. Improvements to the APBS biomolecular solvation software suite. '
    'Protein Sci 27, 112-128 (2018).')

PROPKA = Citation(
    'PROPKA',
    'Olsson MHM, Sondergaard CR, Rostkowski M, Jensen JH. PROPKA3: Consistent Treatment of '
    'Internal and Surface Residues in Empirical pKa Predictions. '
    'J Chem Theory Comput 7, 525-537 (2011).')

ALWAYS = [NAMD, VMD, PSFGEN, CHARMM36M, CHARMM_OVERVIEW]

_WATER = {'TIP3', 'TIP3P', 'WATER', 'SPC', 'SPCE'}


@dataclass
class _Facts:
    """The handful of things about a config that decide which conditional citations apply."""
    tasks: list = field(default_factory=list)
    namd_version: int = 3
    nonwater_solvents: list = field(default_factory=list)
    user_custom: bool = False
    glycan_mods: bool = False


def _facts(config):
    """Reduce a config to the facts the predicates need, tolerating anything missing.

    A citation report must not be the thing that breaks a finished build, so every lookup here is
    defensive and an unreadable config simply yields no conditional citations.
    """
    f = _Facts()
    try:
        user = config['user']
    except Exception:
        return f
    tasks = user.get('tasks') or []
    for entry in tasks:
        if not isinstance(entry, dict):
            continue
        for name, spec in entry.items():
            f.tasks.append(name)
            spec = spec or {}
            if not isinstance(spec, dict):
                continue
            if name == 'solvate':
                solvent = str(spec.get('solvent', 'TIP3') or 'TIP3')
                if solvent.upper() not in _WATER:
                    f.nonwater_solvents.append(solvent)
            elif name == 'psfgen':
                # `mods` belongs to the psfgen *task*, not to the top-level `psfgen:` config
                # block (which holds only segtypes and aliases).  Test `grafts` for a non-empty
                # list rather than for truthiness of the enclosing blocks: the schema
                # materializes defaults for `source.sequence.glycans` on every psfgen task, so
                # that key is always present and always truthy and says nothing about this build.
                grafts = (spec.get('mods') or {}).get('grafts')
                if isinstance(grafts, (list, tuple)) and len(grafts):
                    f.glycan_mods = True
    try:
        f.namd_version = int(user.get('namd', {}).get('namd-version', 3))
    except Exception:
        pass
    custom = (user.get('charmmff') or {}).get('user_custom') or {}
    f.user_custom = any(custom.get(k) for k in ('searchpath', 'rtf', 'str', 'prm'))
    return f


def citations_for(config):
    """The citations this configuration owes, unconditional ones first.

    Returns a list of :class:`Citation`; conditional entries carry the ``reason`` they applied.
    """
    f = _facts(config)
    out = list(ALWAYS)

    if f.namd_version >= 3:
        out.append(NAMD3)

    membrane = [t for t in f.tasks if t in ('make_membrane_system', 'membrane_equilibrate')]
    if membrane:
        out.append(_because(CHARMM_LIPID, f'this build has a {membrane[0]} task'))

    if f.glycan_mods:
        out.append(_because(CHARMM_CARB, 'this build grafts residues (typically glycans) onto the structure'))

    cgenff_reason = ''
    if f.nonwater_solvents:
        cgenff_reason = f"solvated in {', '.join(sorted(set(f.nonwater_solvents)))}, " \
                        f"whose parameters come from CGenFF"
    elif f.user_custom:
        cgenff_reason = 'this build supplies its own CHARMM-format topology/parameter files'
    if cgenff_reason:
        out.append(_because(CGENFF, cgenff_reason))
        if f.user_custom:
            out.append(_because(CGENFF_PROGRAM_I, 'if those parameters came from the CGenFF program'))
            out.append(_because(CGENFF_PROGRAM_II, 'if those parameters came from the CGenFF program'))

    if 'pdb2pqr' in f.tasks:
        why = 'this build has a pdb2pqr task'
        out.append(_because(PDB2PQR, why))
        out.append(_because(APBS, why))
        out.append(_because(PROPKA, f'{why}, which assigns pKa with PROPKA'))

    return out


# --- coordinate sources ----------------------------------------------------------------------

_PIDIBBLE_CITATIONS_MIN = '1.10.0'
"""First pidibble release exposing ``PDBParser.citation_ids()``; below it, structures cannot be
cited and the report says so rather than quietly omitting them."""


@dataclass(frozen=True)
class StructureSource:
    """A set of coordinates a build took as input, and the role it played."""
    id: str
    origin: str
    role: str

    @property
    def citable(self):
        """Whether an experimental paper can be expected behind these coordinates.

        AlphaFold models and files the user supplied locally have no deposited citation of their
        own, so claiming one would be wrong.
        """
        return self.origin == 'rcsb'


def _looks_like_pdb_id(token):
    """PDB IDs are four characters beginning with a digit; anything else names a local file."""
    return len(token) == 4 and token[0].isdigit() and token.isalnum()


def structure_sources(config):
    """Every set of input coordinates this config names, deduplicated, in first-seen order.

    Covers the three places coordinates enter a build: the ``fetch`` task, graft donors, and
    C-terminal fusion donors.  Grafts and fusions are parsed with pestifer's own shortcode
    parsers rather than by splitting strings here, so this cannot drift from the syntax the
    build itself accepts.
    """
    from ..objs.graft import Graft
    from ..objs.cfusion import Cfusion

    found, seen = [], set()

    def add(token, origin, role):
        token = str(token or '').strip()
        if not token or token.lower() in seen:
            return
        seen.add(token.lower())
        found.append(StructureSource(token, origin, role))

    try:
        tasks = config['user'].get('tasks') or []
    except Exception:
        return found

    for entry in tasks:
        if not isinstance(entry, dict):
            continue
        for name, spec in entry.items():
            if not isinstance(spec, dict):
                continue
            if name == 'fetch':
                add(spec.get('sourceID'), str(spec.get('source', 'rcsb') or 'rcsb'),
                    'the structure this build starts from')
            elif name == 'psfgen':
                mods = spec.get('mods') or {}
                for shortcode in (mods.get('grafts') or []):
                    try:
                        pdbid = Graft._from_shortcode(str(shortcode))['source_pdbid']
                    except Exception as e:
                        logger.debug(f'could not read graft source from {shortcode!r}: {e}')
                        continue
                    add(pdbid, 'rcsb' if _looks_like_pdb_id(pdbid) else 'local',
                        'coordinates were grafted from it')
                for shortcode in (mods.get('Cfusions') or []):
                    try:
                        srcfile = Cfusion._from_shortcode(str(shortcode))['sourcefile']
                    except Exception as e:
                        logger.debug(f'could not read fusion source from {shortcode!r}: {e}')
                        continue
                    add(srcfile, 'rcsb' if _looks_like_pdb_id(srcfile) else 'local',
                        'a domain was fused from it')
    return found


def _pidibble_citation_ids(path):
    """DOI/PMID for every citation in the structure file at ``path``, or ``None`` if unavailable.

    Reads the file the build already downloaded rather than re-fetching it: the citation belongs
    to the coordinates that were actually used, and a citation report is no reason to put a
    network call in a build.  ``None`` distinguishes "pidibble cannot do this" from "this entry
    lists no citation".
    """
    try:
        from pidibble.pdbparse import PDBParser
    except Exception as e:
        logger.debug(f'pidibble unavailable: {e}')
        return None
    if not hasattr(PDBParser, 'citation_ids'):
        return None
    try:
        parser = PDBParser(filepath=path).parse()
    except Exception as e:
        logger.debug(f'could not parse {path} for citations: {e}')
        return None
    try:
        return list(parser.citation_ids(role='primary'))
    except Exception as e:
        logger.debug(f'could not read citations from {path}: {e}')
        return None


def _local_structure_file(source_id):
    """The downloaded file for ``source_id`` in the working directory, preferring mmCIF.

    mmCIF preserves the real capitalization of titles and author names; the PDB format does not.
    We report identifiers rather than reference strings either way, but preferring the richer
    format costs nothing and leaves the door open.
    """
    for ext in ('.cif', '.mmcif', '.pdb'):
        for candidate in glob.glob(f'{source_id}{ext}') + glob.glob(f'{source_id.upper()}{ext}') \
                + glob.glob(f'{source_id.lower()}{ext}'):
            if os.path.exists(candidate):
                return candidate
    return None


def structure_citations(config):
    """Citations owed to the input coordinates, plus a note when they cannot be produced.

    Returns ``(citations, notes)``.  A note rather than a silent omission is the point: a report
    that claims to be complete must say when it is not.
    """
    cites, notes = [], []
    sources = structure_sources(config)
    if not sources:
        return cites, notes

    citable = [s for s in sources if s.citable]
    for s in (s for s in sources if not s.citable):
        notes.append(f'{s.id} ({s.origin}) has no deposited citation of its own '
                     f'-- {s.role}')
    if not citable:
        return cites, notes

    if not hasattr(_pidibble_import(), 'citation_ids'):
        # Name the entries anyway.  Knowing *which* structures a build owes citations to is
        # most of the value, and silently dropping them would make the report claim a
        # completeness it does not have.
        notes.append(f'structure citations need pidibble >= {_PIDIBBLE_CITATIONS_MIN}; the '
                     f'installed version cannot resolve them, so the entries below are named '
                     f'but unresolved -- look them up at rcsb.org')
        for s in citable:
            cites.append(Citation(f'PDB {s.id.upper()}', '(citation unresolved)', reason=s.role))
        return cites, notes

    for s in citable:
        path = _local_structure_file(s.id)
        if path is None:
            notes.append(f'{s.id.upper()}: no downloaded structure file found in the working '
                         f'directory, so its citation could not be read')
            continue
        ids = _pidibble_citation_ids(path)
        if ids is None:
            notes.append(f'{s.id.upper()}: {os.path.basename(path)} could not be read for '
                         f'citations (see the debug log)')
            continue
        if not ids:
            notes.append(f'{s.id.upper()}: the structure file lists no citation identifier')
            continue
        for cid in ids:
            bits = []
            if getattr(cid, 'doi', None):
                bits.append(f'doi:{cid.doi}')
            if getattr(cid, 'pmid', None):
                bits.append(f'pmid:{cid.pmid}')
            cites.append(Citation(f'PDB {s.id.upper()}',
                                  ' '.join(bits) or '(no identifier)',
                                  reason=s.role))
    return cites, notes


def _pidibble_import():
    """The PDBParser class, or a stand-in with no citation API, for capability checks."""
    try:
        from pidibble.pdbparse import PDBParser
        return PDBParser
    except Exception:
        return object


def _because(citation, reason):
    """Copy of ``citation`` carrying the reason it applies to this build."""
    return Citation(citation.subject, citation.text, citation.doi, reason)


def log_citations(config, log=None):
    """Report the citations this build owes, and return them.

    Emitted at INFO, so it reaches the console and the diagnostics log alike, and every line
    carries the same prefix so the block can be recovered from a log with one grep.
    """
    log = log or logger.info
    try:
        cites = citations_for(config)
    except Exception as e:                      # never fail a finished build over a report
        logger.debug(f'citation report failed: {e}')
        return []
    log('citations: ---- please cite ----')
    log('citations: Pestifer orchestrates other people\'s software.  Work that uses a system')
    log('citations: prepared by this build should cite the following:')
    for c in cites:
        log(f'citations: {c.render()}')
        if c.reason:
            log(f'citations:     ^ included because {c.reason}')
    struct, notes = structure_citations(config)
    if struct or notes:
        log('citations: --- and the coordinates this build was given ---')
    for c in struct:
        log(f'citations: {c.render()}')
        if c.reason:
            log(f'citations:     ^ {c.reason}')
    for n in notes:
        log(f'citations: note: {n}')
    log('citations: Software citations are derived from your input config, so they cover what')
    log('citations: you asked for.  A structure that already carries glycans, nucleic acids or')
    log('citations: ligands may owe further CHARMM parameter-set citations.')
    log('citations: ------------------------')
    return cites + struct
