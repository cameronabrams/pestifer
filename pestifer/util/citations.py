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

Sources for the reference strings: the CHARMM/NAMD/VMD/psfgen entries match the bibliography of
the pestifer paper; the PDB2PQR and PROPKA entries are the citations those packages themselves
declare (``pdb2pqr.config.CITATIONS`` and ``propka.__doc__``).
"""

import logging
from dataclasses import dataclass, field
from typing import Callable, Optional

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
    log('citations: This list is derived from your input config, so it covers what you asked')
    log('citations: for.  A structure that already carries glycans, nucleic acids or ligands')
    log('citations: may owe further CHARMM parameter-set citations; see the force-field files.')
    log('citations: ------------------------')
    return cites
