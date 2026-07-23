# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Structure inspection for ``pestifer new-system --inspect``.

Given a PDB/mmCIF/UniProt id, fetch and parse the structure (via pidibble, no ``Molecule`` build
needed) and discover the sequence features a user typically has to read out of the header by hand:

- **Missing residues** (``REMARK 465`` / mmCIF ``pdbx_unobs_or_zero_occ_residues``), grouped into
  contiguous runs and classified -- by comparison with each chain's resolved residue span -- as an
  interior gap (both anchors resolved), or an N-/C-terminal tail (a free end).
- **Engineered mutations / conflicts** (``SEQADV`` / mmCIF ``struct_ref_seq_dif``), where the
  modeled residue differs from the sequence-database reference.
- **Cloning artifacts / expression tags** (``SEQADV`` typed ``cloning`` / ``expression tag``) --
  extra residues not part of the native protein.

The findings are rendered as commented mod stubs the user can uncomment in the generated config
(see :meth:`Findings.annotation_lines`); interior gaps additionally justify an active ``ligate``
task, since interior missing loops are built and closed automatically.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class MissingRun:
    """A contiguous run of unresolved residues in one chain."""
    chain: str
    start: int
    end: int
    kind: str          # 'N', 'C', or 'interior'

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        span = f'{self.start}' if self.start == self.end else f'{self.start}-{self.end}'
        return f'chain {self.chain}: {span} ({self.length} residue{"s" if self.length != 1 else ""})'


@dataclass
class MutationFinding:
    """A modeled residue that differs from the sequence-database reference (SEQADV)."""
    chain: str
    resid: int
    resname: str       # residue present in the structure
    dbres: str         # database reference residue
    typekey: str       # 'engineered mutation' | 'conflict'

    @property
    def revert_shortcode(self) -> str:
        """The ``mutations`` shortcode that reverts this residue to the database reference."""
        return f'{self.chain}:{self.resname},{self.resid},{self.dbres}'


@dataclass
class ExcisionRun:
    """A contiguous run of cloning-artifact / expression-tag residues (SEQADV)."""
    chain: str
    start: int
    end: int
    typekey: str       # 'cloning' | 'expression tag'

    @property
    def delete_shortcode(self) -> str:
        span = f'{self.start}' if self.start == self.end else f'{self.start}-{self.end}'
        return f'{self.chain}:{span}'


@dataclass
class ChainIdentity:
    """What a single chain is: its dominant segment type, size, and a resname sample."""
    chain: str
    segtype: str
    n_residues: int
    resnames: list = field(default_factory=list)
    auth: str = None          # mmCIF author chain id (label chain is `chain`); None for PDB

    def describe(self) -> str:
        sample = ', '.join(self.resnames[:4]) + ('...' if len(self.resnames) > 4 else '')
        if self.segtype == 'water':
            return 'water'
        if self.segtype == 'ion':
            return f'ion ({sample})'
        if self.segtype == 'glycan':
            return f'glycan ({sample})'
        if self.segtype in ('protein', 'nucleicacid'):
            name = 'protein' if self.segtype == 'protein' else 'nucleic acid'
            return f'{name} ({self.n_residues} residues)'
        return f'{self.segtype} ({sample})'


@dataclass
class AssemblyInfo:
    """A biological assembly: its 1-based index, how many copies it generates, and its chains."""
    index: int
    n_copies: int
    chains: list = field(default_factory=list)


@dataclass
class Findings:
    db_id: str
    source_format: str = 'pdb'
    missing_runs: list = field(default_factory=list)
    mutations: list = field(default_factory=list)
    excisions: list = field(default_factory=list)
    chains: list = field(default_factory=list)
    assemblies: list = field(default_factory=list)

    def is_empty(self) -> bool:
        return not (self.missing_runs or self.mutations or self.excisions
                    or len(self.chains) > 1 or len(self.assemblies) > 1)

    def interior_gaps(self) -> list:
        return [r for r in self.missing_runs if r.kind == 'interior']

    def terminal_tail_chains(self) -> dict:
        """``{'n': [chains], 'c': [chains]}`` for chains with an unresolved N-/C-terminal run."""
        out = {'n': [], 'c': []}
        for r in self.missing_runs:
            if r.kind == 'N' and r.chain not in out['n']:
                out['n'].append(r.chain)
            elif r.kind == 'C' and r.chain not in out['c']:
                out['c'].append(r.chain)
        return out

    def annotation_lines(self, indent: str = '    ') -> list:
        """
        Render the findings as YAML comment lines (each prefixed with ``indent`` + ``# ``), meant to
        be injected at the **psfgen-task level** (a sibling of ``source:``). Each block lists the
        detected items and a ready-to-uncomment stub, with the correct nesting shown: ``sequence:``
        blocks go *under* ``source:``; ``mods:`` is a *sibling* of ``source:`` (both under psfgen).
        """
        L = []

        def c(text=''):
            L.append(f'{indent}#{" " + text if text else ""}')

        c('=== pestifer new-system --inspect: detected structure features ===')
        c('Uncomment/paste the block(s) you want; nesting is shown relative to the psfgen task.')

        if self.assemblies:
            c()
            c('Biological assemblies (set `biological_assembly:` under `source:`; 0 = asymmetric unit):')
            for a in self.assemblies:
                c(f'  {a.index}: {a.n_copies} cop{"y" if a.n_copies == 1 else "ies"} of '
                  f'chain(s) [{", ".join(a.chains)}]')
            c(f'The template uses `biological_assembly: 1`; change it to another index or to 0.')

        if self.chains:
            c()
            c('Chains present (reference these ids in `exclude:`; label ids for mmCIF):')
            for ch in self.chains:
                auth = f' (author {ch.auth})' if ch.auth and ch.auth != ch.chain else ''
                c(f'  {ch.chain}{auth}: {ch.describe()}')
            c('To OMIT chains, add under `source:` -')
            c('  exclude:')
            c(f"    - chainID in ['X', 'Y']")

        tails = self.terminal_tail_chains()
        term_runs = [r for r in self.missing_runs if r.kind in ('N', 'C')]
        if term_runs:
            c()
            c('Missing TERMINAL residues (unresolved; often disordered ends or tags):')
            for r in term_runs:
                c(f'  {r} [{"N" if r.kind == "N" else "C"}-terminus]')
            c('To BUILD any as a modeled tail (dropped by default), add under `source:` -')
            c('  sequence:')
            c('    terminal_tails:')
            if tails['n']:
                c(f'      n: [{", ".join(tails["n"])}]')
            if tails['c']:
                c(f'      c: [{", ".join(tails["c"])}]')

        gaps = self.interior_gaps()
        if gaps:
            c()
            c('Interior missing LOOPS (built and closed automatically; a `ligate` task was added):')
            for r in gaps:
                c(f'  {r}')
            c('Short gaps are skipped unless >= `source.sequence.loops.min_loop_length` (default 4).')

        revertible = [m for m in self.mutations if m.dbres]
        if self.mutations:
            c()
            c('Engineered mutations / conflicts (structure differs from the sequence database):')
            for m in self.mutations:
                c(f'  {m.chain}:{m.resname}{m.resid} (db {m.dbres or "?"}) [{m.typekey}]')
            if revertible:
                c('To REVERT specific residues, add a `mods:` block (a SIBLING of `source:`) -')
                c('  mods:')
                c(f'    mutations: [{", ".join(m.revert_shortcode for m in revertible)}]')
                c('Or revert ALL engineered mutations at once, add under `source:` -')
                c('  sequence:')
                c('    fix_engineered_mutations: true')

        if self.excisions:
            c()
            c('Cloning artifacts / expression tags (extra residues not in the native protein):')
            for e in self.excisions:
                span = f'{e.start}' if e.start == e.end else f'{e.start}-{e.end}'
                c(f'  chain {e.chain}: {span} [{e.typekey}]')
            c('To EXCISE them, add a `mods:` block (a SIBLING of `source:`) -')
            c('  mods:')
            c(f'    deletions: [{", ".join(e.delete_shortcode for e in self.excisions)}]')

        c('=== end detected features ===')
        return L


def make_prompter():
    """A default yes/no prompter reading from stdin; EOF/blank returns the question's default."""
    def ask(question, default=False):
        suffix = '[Y/n]' if default else '[y/N]'
        try:
            resp = input(f'{question} {suffix} ').strip().lower()
        except EOFError:
            return default
        if not resp:
            return default
        return resp in ('y', 'yes')
    return ask


def interactive_select(findings: 'Findings', ask=None, say=None) -> dict:
    """
    Walk the user through the findings, prompting per item, and return the chosen mods as
    ``{'sequence': {...}, 'mods': {...}, 'add_ligate': bool}`` -- the *active* (uncommented) blocks
    to inject into the generated psfgen task.

    ``ask(question, default) -> bool`` and ``say(message)`` are injectable (default: stdin prompter
    and ``print``) so the walkthrough is testable without a TTY. Declining everything yields empty
    ``sequence``/``mods`` and ``add_ligate`` False.
    """
    ask = ask or make_prompter()
    say = say or (lambda m: print(m))
    sequence, mods, source = {}, {}, {}

    # biological assembly
    if findings.assemblies:
        say('\nBiological assemblies:')
        for a in findings.assemblies:
            say(f'  {a.index}: {a.n_copies} cop{"y" if a.n_copies == 1 else "ies"} of '
                f'chain(s) [{", ".join(a.chains)}]')
        if len(findings.assemblies) == 1:
            a = findings.assemblies[0]
            if a.n_copies > 1 and not ask(f'  Build biological assembly {a.index} '
                                          f'({a.n_copies} copies)? (No = asymmetric unit only)', True):
                source['biological_assembly'] = 0
            else:
                source['biological_assembly'] = a.index
        else:
            chosen = 0
            for a in findings.assemblies:
                if ask(f'  Build biological assembly {a.index} '
                       f'({a.n_copies} copies of [{", ".join(a.chains)}])?', False):
                    chosen = a.index
                    break
            source['biological_assembly'] = chosen
    else:
        source['biological_assembly'] = 0        # no assemblies defined -> asymmetric unit

    # chain omission
    if len(findings.chains) > 1:
        say('\nChains present:')
        omit = []
        for ch in findings.chains:
            auth = f' (author {ch.auth})' if ch.auth and ch.auth != ch.chain else ''
            if ask(f'  Omit chain {ch.chain}{auth} [{ch.describe()}]?', False):
                omit.append(ch.chain)
        if omit:
            source['exclude'] = [f"chainID in [{', '.join(repr(c) for c in omit)}]"]

    tails = {'n': [], 'c': []}
    term_runs = [r for r in findings.missing_runs if r.kind in ('N', 'C')]
    if term_runs:
        say('\nMissing terminal residues (unresolved; dropped unless you build them):')
        for r in term_runs:
            end = 'N' if r.kind == 'N' else 'C'
            if ask(f'  Build a modeled tail for chain {r.chain} {end}-terminus '
                   f'({r.start}-{r.end}, {r.length} res)?', False):
                key = 'n' if r.kind == 'N' else 'c'
                if r.chain not in tails[key]:
                    tails[key].append(r.chain)
    tt = {k: v for k, v in tails.items() if v}
    if tt:
        sequence['terminal_tails'] = tt

    gaps = findings.interior_gaps()
    add_ligate = False
    if gaps:
        say(f'\n{len(gaps)} interior missing loop(s) will be built and closed:')
        for r in gaps:
            say(f'  {r}')
        add_ligate = ask('  Add a `ligate` task to close them?', True)

    revertible = [m for m in findings.mutations if m.dbres]
    if revertible:
        say('\nEngineered mutations / conflicts (structure differs from the sequence database):')
        muts = []
        for m in revertible:
            if ask(f'  Revert {m.chain}:{m.resname}{m.resid} -> db {m.dbres} [{m.typekey}]?', False):
                muts.append(m.revert_shortcode)
        if muts:
            mods['mutations'] = muts

    if findings.excisions:
        say('\nCloning artifacts / expression tags (extra residues not in the native protein):')
        dels = []
        for e in findings.excisions:
            span = f'{e.start}' if e.start == e.end else f'{e.start}-{e.end}'
            if ask(f'  Excise {e.typekey} in chain {e.chain} ({span})?', False):
                dels.append(e.delete_shortcode)
        if dels:
            mods['deletions'] = dels

    return {'source': source, 'sequence': sequence, 'mods': mods, 'add_ligate': add_ligate}


def _group_runs(resseqnums):
    """Yield ``(start, end)`` for each maximal run of consecutive integers in ``resseqnums``."""
    s = sorted(set(resseqnums))
    if not s:
        return
    start = prev = s[0]
    for n in s[1:]:
        if n - prev <= 1:
            prev = n
        else:
            yield (start, prev)
            start = prev = n
    yield (start, prev)


def inspect_structure(db_id: str, source_format: str = 'pdb', source_db: str = 'rcsb') -> Findings:
    """
    Fetch and inspect a structure, returning the discovered :class:`Findings`.

    ``source_format`` is ``'pdb'`` or ``'cif'``; ``source_db`` is ``'rcsb'`` (default) or
    ``'alphafold'`` (UniProt id). Uses only the lightweight pidibble parse -- no ``Molecule``.
    """
    from pidibble.pdbparse import PDBParser
    from ..molecule.residue import ResiduePlaceholderList
    from ..molecule.atom import AtomList
    from ..objs.seqadv import SeqadvList

    kw = {'source_id': db_id, 'source_db': source_db}
    if source_format == 'cif':
        kw['input_format'] = 'mmCIF'
    parsed = PDBParser(**kw).parse().parsed

    is_cif = source_format == 'cif'
    missings = ResiduePlaceholderList.from_cif(parsed) if is_cif else ResiduePlaceholderList.from_pdb(parsed)
    seqadvs = SeqadvList.from_cif(parsed) if is_cif else SeqadvList.from_pdb(parsed)
    atoms = AtomList.from_cif(parsed) if is_cif else AtomList.from_pdb(parsed)

    from ..core.labels import Labels

    resolved = defaultdict(set)
    for a in atoms.data:
        resolved[a.chainID].add(a.resid.resseqnum)
    missing_pairs = [(m.chainID, m.resid.resseqnum) for m in missings.data]
    seqadv_tuples = [(s.chainID, getattr(s.resid, 'resseqnum', None), s.resname, s.dbRes or '',
                      (s.typekey or '').strip()) for s in seqadvs.data]

    # chain identities (segtype by resname; chainID is the label id for mmCIF, author for PDB)
    from collections import Counter
    chain_data = {}
    for a in atoms.data:
        d = chain_data.setdefault(a.chainID, {'residues': set(), 'resnames': [],
                                              'auth': getattr(a, 'auth_asym_id', None)})
        d['residues'].add((a.resid.resseqnum, getattr(a.resid, 'insertion', '')))
        if a.resname not in d['resnames']:
            d['resnames'].append(a.resname)
    chains = []
    for ch, d in chain_data.items():
        segs = Counter(Labels.segtype_of_resname.get(rn, 'other') for rn in d['resnames'])
        chains.append(ChainIdentity(ch, segs.most_common(1)[0][0], len(d['residues']),
                                    d['resnames'], d['auth']))
    chains.sort(key=lambda c: c.chain)

    # biological assemblies
    assemblies = []
    try:
        from ..molecule.bioassemb import BioAssembList
        for ba in BioAssembList(parsed):
            tr = getattr(ba, 'transforms', None)
            chs = list(tr.data[0].applies_chainIDs) if tr and len(tr.data) else []
            assemblies.append(AssemblyInfo(ba.index, len(tr.data) if tr else 0, chs))
    except Exception as e:                        # never let assembly parsing break inspection
        logger.debug(f'--inspect: biological-assembly enumeration failed: {e}')
    assemblies.sort(key=lambda a: a.index)

    return build_findings(db_id, source_format, missing_pairs, seqadv_tuples,
                          {c: (min(v), max(v)) for c, v in resolved.items()},
                          chains=chains, assemblies=assemblies)


def build_findings(db_id, source_format, missing_pairs, seqadv_tuples, resolved_ranges,
                   chains=(), assemblies=()) -> Findings:
    """
    Assemble :class:`Findings` from already-extracted primitives -- the pure, network-free core of
    :func:`inspect_structure`.

    Parameters
    ----------
    missing_pairs : list[(chain, resseqnum)]
        Every unresolved residue.
    seqadv_tuples : list[(chain, resseqnum, resname, dbres, typekey)]
        Every SEQADV entry (``resseqnum`` may be ``None`` -> skipped).
    resolved_ranges : dict[chain, (lo, hi)]
        The min/max resolved sequence number per chain (empty/absent -> a run is treated interior).
    """
    findings = Findings(db_id=db_id, source_format=source_format,
                        chains=list(chains), assemblies=list(assemblies))

    by_chain = defaultdict(list)
    for chain, num in missing_pairs:
        by_chain[chain].append(num)
    for chain, nums in by_chain.items():
        rng = resolved_ranges.get(chain)
        for start, end in _group_runs(nums):
            if not rng:
                kind = 'interior'            # no resolved residues to anchor against; treat as loop
            elif end < rng[0]:
                kind = 'N'
            elif start > rng[1]:
                kind = 'C'
            else:
                kind = 'interior'
            findings.missing_runs.append(MissingRun(chain, start, end, kind))

    tag_nums = defaultdict(list)
    for chain, resid, resname, dbres, tk in seqadv_tuples:
        if resid is None:
            continue
        if tk in ('engineered mutation', 'conflict'):
            findings.mutations.append(MutationFinding(chain, resid, resname, dbres, tk))
        elif tk in ('cloning', 'expression tag'):
            tag_nums[(chain, tk)].append(resid)
    for (chain, tk), nums in tag_nums.items():
        for start, end in _group_runs(nums):
            findings.excisions.append(ExcisionRun(chain, start, end, tk))

    findings.missing_runs.sort(key=lambda r: (r.chain, r.start))
    findings.mutations.sort(key=lambda m: (m.chain, m.resid))
    findings.excisions.sort(key=lambda e: (e.chain, e.start))
    return findings
