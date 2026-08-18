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

import contextlib
import glob
import logging
import os
import tarfile
import tempfile
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Citation:
    """One reference, what it is owed for, and (when conditional) why it applies here.

    ``bibtex`` is the entry a Methods draft emits.  It is written out by hand for the software
    citations, whose bibliographic details are fixed and known; coordinate citations get a
    minimal ``@misc`` carrying the DOI instead, because a PDB-format file cannot tell us the real
    capitalization of a title or an author name and a guess would be wrong in a bibliography.
    """
    subject: str
    text: str
    doi: str = ''
    reason: str = ''
    key: str = ''
    bibtex: str = ''

    def render(self):
        line = f'[{self.subject}] {self.text}'
        if self.doi:
            line += f' doi:{self.doi}'
        return line

    def bib_entry(self):
        """This citation as a BibTeX entry, falling back to a DOI-only ``@misc``.

        The fallback is deliberately thin: an entry that carries only an identifier can be
        resolved correctly by the author's reference manager, whereas one padded with guessed
        metadata looks complete and is not.
        """
        if self.bibtex:
            return self.bibtex.strip()
        key = self.key or _slug(self.subject)
        fields = [f'  note = {{{self.subject}}}']
        if self.doi:
            fields.insert(0, f'  doi = {{{self.doi}}}')
            fields.insert(0, f'  howpublished = {{\\url{{https://doi.org/{self.doi}}}}}')
        return '@misc{' + key + ',\n' + ',\n'.join(fields) + '\n}'


def _slug(text):
    """A BibTeX-safe key fragment."""
    return ''.join(c if c.isalnum() else '-' for c in text.lower()).strip('-')


# --- the catalog -----------------------------------------------------------------------------
# Unconditional: every build runs NAMD, drives VMD/psfgen, and uses the CHARMM36 protein set.

NAMD = Citation(
    'NAMD',
    'Phillips JC, et al. Scalable molecular dynamics with NAMD. '
    'J Comput Chem 26, 1781-1802 (2005).', '10.1002/jcc.20289',
    key='namd2005',
    bibtex=r'''@article{namd2005,
  author  = {Phillips, James C. and Braun, Rosemary and Wang, Wei and Gumbart, James and
             Tajkhorshid, Emad and Villa, Elizabeth and Chipot, Christophe and Skeel, Robert D.
             and Kal{\'e}, Laxmikant and Schulten, Klaus},
  title   = {Scalable molecular dynamics with {NAMD}},
  journal = {Journal of Computational Chemistry},
  volume  = {26}, pages = {1781--1802}, year = {2005},
  doi     = {10.1002/jcc.20289}
}''')

NAMD3 = Citation(
    'NAMD 3',
    'Phillips JC, et al. Scalable molecular dynamics on CPU and GPU architectures with NAMD. '
    'J Chem Phys 153, 044130 (2020).', '10.1063/5.0014475',
    reason='this build ran NAMD 3',
    key='namd2020',
    bibtex=r'''@article{namd2020,
  author  = {Phillips, James C. and Hardy, David J. and Maia, Julio D. C. and Stone, John E. and
             Ribeiro, Jo{\~a}o V. and Bernardi, Rafael C. and Buch, Ronak and Fiorin, Giacomo and
             H{\'e}nin, J{\'e}r{\^o}me and Jiang, Wei and McGreevy, Ryan and Melo, Marcelo C. R.
             and Radak, Brian K. and Skeel, Robert D. and Singharoy, Abhishek and Wang, Yi and
             Roux, Beno{\^i}t and Aksimentiev, Aleksei and Luthey-Schulten, Zaida and
             Kal{\'e}, Laxmikant V. and Schulten, Klaus and Chipot, Christophe and
             Tajkhorshid, Emad},
  title   = {Scalable molecular dynamics on {CPU} and {GPU} architectures with {NAMD}},
  journal = {The Journal of Chemical Physics},
  volume  = {153}, number = {4}, pages = {044130}, year = {2020},
  doi     = {10.1063/5.0014475}
}''')

VMD = Citation(
    'VMD',
    'Humphrey W, Dalke A, Schulten K. VMD: Visual molecular dynamics. '
    'J Mol Graph 14, 33-38 (1996).', '10.1016/0263-7855(96)00018-5',
    key='vmd1996',
    bibtex=r'''@article{vmd1996,
  author  = {Humphrey, William and Dalke, Andrew and Schulten, Klaus},
  title   = {{VMD}: Visual molecular dynamics},
  journal = {Journal of Molecular Graphics},
  volume  = {14}, pages = {33--38}, year = {1996},
  doi     = {10.1016/0263-7855(96)00018-5}
}''')

PSFGEN = Citation(
    'psfgen',
    "Ribeiro J, Radak B, Stone J, Gullingsrud J, Saam J, Phillips J. "
    "psfgen User's Guide, v2.0 (2020).",
    key='psfgen2020',
    bibtex=r'''@misc{psfgen2020,
  author       = {Ribeiro, Jo{\~a}o and Radak, Brian and Stone, John and Gullingsrud, Justin and
                  Saam, Jan and Phillips, Jim},
  title        = {psfgen User's Guide, v2.0},
  year         = {2020},
  howpublished = {Theoretical and Computational Biophysics Group, University of Illinois}
}''')

CHARMM36M = Citation(
    'CHARMM36m (protein)',
    'Huang J, et al. CHARMM36m: An improved force field for folded and intrinsically disordered '
    'proteins. Nat Methods 14, 71-73 (2016).', '10.1038/nmeth.4067',
    key='charmm36m',
    bibtex=r'''@article{charmm36m,
  author  = {Huang, Jing and Rauscher, Sarah and Nawrocki, Grzegorz and Ran, Ting and Feig, Michael
             and de Groot, Bert L. and Grubm{\"u}ller, Helmut and MacKerell, Alexander D.},
  title   = {{CHARMM36m}: An improved force field for folded and intrinsically disordered proteins},
  journal = {Nature Methods},
  volume  = {14}, pages = {71--73}, year = {2017},
  doi     = {10.1038/nmeth.4067}
}''')

CHARMM_OVERVIEW = Citation(
    'CHARMM force field',
    'Vanommeslaeghe K, MacKerell AD. CHARMM additive and polarizable force fields for biophysics '
    'and computer-aided drug design. Biochim Biophys Acta 1850, 861-871 (2015).',
    '10.1016/j.bbagen.2014.08.004',
    key='charmm-ff-overview',
    bibtex=r'''@article{charmm-ff-overview,
  author  = {Vanommeslaeghe, K. and MacKerell, A. D.},
  title   = {{CHARMM} additive and polarizable force fields for biophysics and computer-aided drug design},
  journal = {Biochimica et Biophysica Acta - General Subjects},
  volume  = {1850}, pages = {861--871}, year = {2015},
  doi     = {10.1016/j.bbagen.2014.08.004}
}''')

CHARMM_LIPID = Citation(
    'CHARMM36 (lipids)',
    'Klauda JB, et al. Update of the CHARMM All-Atom Additive Force Field for Lipids: Validation '
    'on Six Lipid Types. J Phys Chem B 114, 7830-7843 (2010).', '10.1021/jp101759q',
    key='charmm36-lipid',
    bibtex=r'''@article{charmm36-lipid,
  author  = {Klauda, Jeffery B. and Venable, Richard M. and Freites, J. Alfredo and
             O'Connor, Joseph W. and Tobias, Douglas J. and Mondragon-Ramirez, Carlos and
             Vorobyov, Igor and MacKerell, Alexander D. and Pastor, Richard W.},
  title   = {Update of the {CHARMM} All-Atom Additive Force Field for Lipids: Validation on Six Lipid Types},
  journal = {The Journal of Physical Chemistry B},
  volume  = {114}, pages = {7830--7843}, year = {2010},
  doi     = {10.1021/jp101759q}
}''')

CHARMM_CARB = Citation(
    'CHARMM36 (carbohydrates)',
    'Guvench O, et al. CHARMM Additive All-Atom Force Field for Glycosidic Linkages between '
    'Hexopyranoses. J Chem Theory Comput 5, 2353-2370 (2009).', '10.1021/ct900242e',
    key='charmm36-carb',
    bibtex=r'''@article{charmm36-carb,
  author  = {Guvench, Olgun and Hatcher, Elizabeth and Venable, Richard M. and Pastor, Richard W.
             and MacKerell, Alexander D.},
  title   = {{CHARMM} Additive All-Atom Force Field for Glycosidic Linkages between Hexopyranoses},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {5}, pages = {2353--2370}, year = {2009},
  doi     = {10.1021/ct900242e}
}''')

CGENFF = Citation(
    'CGenFF',
    'Vanommeslaeghe K, et al. CHARMM general force field: A force field for drug-like molecules '
    'compatible with the CHARMM all-atom additive biological force fields. '
    'J Comput Chem 31, 671-690 (2010).', '10.1002/jcc.21367',
    key='cgenff2010',
    bibtex=r'''@article{cgenff2010,
  author  = {Vanommeslaeghe, K. and Hatcher, E. and Acharya, C. and Kundu, S. and Zhong, S. and
             Shim, J. and Darian, E. and Guvench, O. and Lopes, P. and Vorobyov, I. and
             MacKerell, A. D.},
  title   = {{CHARMM} general force field: A force field for drug-like molecules compatible with
             the {CHARMM} all-atom additive biological force fields},
  journal = {Journal of Computational Chemistry},
  volume  = {31}, pages = {671--690}, year = {2010},
  doi     = {10.1002/jcc.21367}
}''')

CGENFF_PROGRAM_I = Citation(
    'CGenFF program',
    'Vanommeslaeghe K, MacKerell AD. Automation of the CHARMM General Force Field (CGenFF) I: '
    'Bond Perception and Atom Typing. J Chem Inf Model 52, 3144-3154 (2012).', '10.1021/ci300363c',
    key='cgenff-program-1',
    bibtex=r'''@article{cgenff-program-1,
  author  = {Vanommeslaeghe, K. and MacKerell, A. D.},
  title   = {Automation of the {CHARMM} General Force Field ({CGenFF}) {I}: Bond Perception and Atom Typing},
  journal = {Journal of Chemical Information and Modeling},
  volume  = {52}, pages = {3144--3154}, year = {2012},
  doi     = {10.1021/ci300363c}
}''')

CGENFF_PROGRAM_II = Citation(
    'CGenFF program',
    'Vanommeslaeghe K, Raman EP, MacKerell AD. Automation of the CHARMM General Force Field '
    '(CGenFF) II: Assignment of Bonded Parameters and Partial Atomic Charges. '
    'J Chem Inf Model 52, 3155-3168 (2012).', '10.1021/ci3003649',
    key='cgenff-program-2',
    bibtex=r'''@article{cgenff-program-2,
  author  = {Vanommeslaeghe, K. and Raman, E. Prabhu and MacKerell, A. D.},
  title   = {Automation of the {CHARMM} General Force Field ({CGenFF}) {II}: Assignment of Bonded
             Parameters and Partial Atomic Charges},
  journal = {Journal of Chemical Information and Modeling},
  volume  = {52}, pages = {3155--3168}, year = {2012},
  doi     = {10.1021/ci3003649}
}''')

PDB2PQR = Citation(
    'PDB2PQR',
    'Dolinsky TJ, et al. PDB2PQR: expanding and upgrading automated preparation of biomolecular '
    'structures for molecular simulations. Nucleic Acids Res 35, W522-W525 (2007).',
    key='pdb2pqr2007',
    bibtex=r'''@article{pdb2pqr2007,
  author  = {Dolinsky, Todd J. and Czodrowski, Paul and Li, Hui and Nielsen, Jens E. and
             McCammon, J. Andrew and Klebe, Gerhard and Baker, Nathan A.},
  title   = {{PDB2PQR}: expanding and upgrading automated preparation of biomolecular structures
             for molecular simulations},
  journal = {Nucleic Acids Research},
  volume  = {35}, pages = {W522--W525}, year = {2007},
  doi     = {10.1093/nar/gkm276}
}''')

APBS = Citation(
    'APBS/PDB2PQR',
    'Jurrus E, et al. Improvements to the APBS biomolecular solvation software suite. '
    'Protein Sci 27, 112-128 (2018).',
    key='apbs2018',
    bibtex=r'''@article{apbs2018,
  author  = {Jurrus, Elizabeth and Engel, Dave and Star, Keith and Monson, Kyle and Brandi, Juan
             and Felberg, Lisa E. and Brookes, David H. and Wilson, Leighton and Chen, Jiahui and
             Liles, Karina and Chun, Minju and Li, Peter and Gohara, David W. and Dolinsky, Todd
             and Konecny, Robert and Koes, David R. and Nielsen, Jens Erik and Head-Gordon, Teresa
             and Geng, Weihua and Krasny, Robert and Wei, Guo-Wei and Holst, Michael J. and
             McCammon, J. Andrew and Baker, Nathan A.},
  title   = {Improvements to the {APBS} biomolecular solvation software suite},
  journal = {Protein Science},
  volume  = {27}, pages = {112--128}, year = {2018},
  doi     = {10.1002/pro.3280}
}''')

PROPKA = Citation(
    'PROPKA',
    'Olsson MHM, Sondergaard CR, Rostkowski M, Jensen JH. PROPKA3: Consistent Treatment of '
    'Internal and Surface Residues in Empirical pKa Predictions. '
    'J Chem Theory Comput 7, 525-537 (2011).',
    key='propka2011',
    bibtex=r'''@article{propka2011,
  author  = {Olsson, Mats H. M. and S{\o}ndergaard, Chresten R. and Rostkowski, Michal and
             Jensen, Jan H.},
  title   = {{PROPKA3}: Consistent Treatment of Internal and Surface Residues in Empirical
             p{K}a Predictions},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {7}, pages = {525--537}, year = {2011},
  doi     = {10.1021/ct100578z}
}''')

ALWAYS = [NAMD, VMD, PSFGEN, CHARMM36M, CHARMM_OVERVIEW]

#: Every catalog entry that carries a stable BibTeX key, indexed by that key.  A run record
#: stores only the key, not the entry -- so a Methods draft rendered by a later pestifer picks up
#: any corrections made to the bibliography since the build ran, and records stay small.
CATALOG = {c.key: c for c in (
    NAMD, NAMD3, VMD, PSFGEN, CHARMM36M, CHARMM_OVERVIEW, CHARMM_LIPID, CHARMM_CARB,
    CGENFF, CGENFF_PROGRAM_I, CGENFF_PROGRAM_II, PDB2PQR, APBS, PROPKA) if c.key}


def by_key(key):
    """The catalog entry for ``key``, or ``None``."""
    return CATALOG.get(key)

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


def _candidate_names(source_id):
    """Filenames the downloaded structure for ``source_id`` could have, mmCIF first.

    mmCIF preserves the real capitalization of titles and author names; the PDB format does not.
    We report identifiers rather than reference strings either way, but preferring the richer
    format costs nothing and leaves the door open.
    """
    names = []
    for ext in ('.cif', '.mmcif', '.pdb'):
        for stem in (source_id, source_id.lower(), source_id.upper()):
            if f'{stem}{ext}' not in names:
                names.append(f'{stem}{ext}')
    return names


def _local_structure_file(source_id):
    """The downloaded file for ``source_id`` in the working directory, or ``None``."""
    for name in _candidate_names(source_id):
        if os.path.exists(name):
            return name
    return None


def _archived_structure_member(source_id):
    """``(tarball, member)`` for ``source_id`` inside a run's artifacts archive, or ``None``.

    The ``terminate`` task tidies the run directory by sweeping every intermediate -- the fetched
    structures included -- into ``<basename>-artifacts.tar.gz``.  Since the citation report is
    emitted after the run, the originals are typically gone by the time it looks, so it reads the
    archive instead.  That is still the file the build used, and it needs no network.
    """
    wanted = {n.lower() for n in _candidate_names(source_id)}
    for tarball in sorted(glob.glob('*-artifacts.tar.gz')):
        try:
            with tarfile.open(tarball) as tf:
                for member in tf.getnames():
                    if os.path.basename(member).lower() in wanted:
                        return tarball, member
        except Exception as e:
            logger.debug(f'could not inspect {tarball}: {e}')
    return None


@contextlib.contextmanager
def _structure_file(source_id):
    """Yield a readable path to ``source_id``'s structure file, or ``None``.

    Prefers a file still in the working directory; falls back to extracting it from the run's
    artifacts archive into a temporary directory that is cleaned up on exit.
    """
    path = _local_structure_file(source_id)
    if path is not None:
        yield path
        return
    found = _archived_structure_member(source_id)
    if found is None:
        yield None
        return
    tarball, member = found
    with tempfile.TemporaryDirectory() as td:
        try:
            with tarfile.open(tarball) as tf:
                tf.extract(member, path=td, filter='data')
            yield os.path.join(td, member)
        except Exception as e:
            logger.debug(f'could not extract {member} from {tarball}: {e}')
            yield None


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
        with _structure_file(s.id) as path:
            if path is None:
                notes.append(f'{s.id.upper()}: its structure file is neither in the working '
                             f'directory nor in a *-artifacts.tar.gz here, so its citation '
                             f'could not be read')
                continue
            ids = _pidibble_citation_ids(path)
        if ids is None:
            notes.append(f'{s.id.upper()}: the structure file could not be read for citations '
                         f'(see the debug log)')
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
