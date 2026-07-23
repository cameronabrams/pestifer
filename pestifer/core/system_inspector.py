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
class Findings:
    db_id: str
    source_format: str = 'pdb'
    missing_runs: list = field(default_factory=list)
    mutations: list = field(default_factory=list)
    excisions: list = field(default_factory=list)

    def is_empty(self) -> bool:
        return not (self.missing_runs or self.mutations or self.excisions)

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

        c('=== pestifer new-system --inspect: detected sequence features ===')
        c('Uncomment/paste the block(s) you want; nesting is shown relative to the psfgen task.')

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

    resolved = defaultdict(set)
    for a in atoms.data:
        resolved[a.chainID].add(a.resid.resseqnum)
    missing_pairs = [(m.chainID, m.resid.resseqnum) for m in missings.data]
    seqadv_tuples = [(s.chainID, getattr(s.resid, 'resseqnum', None), s.resname, s.dbRes or '',
                      (s.typekey or '').strip()) for s in seqadvs.data]
    return build_findings(db_id, source_format, missing_pairs, seqadv_tuples,
                          {c: (min(v), max(v)) for c, v in resolved.items()})


def build_findings(db_id, source_format, missing_pairs, seqadv_tuples, resolved_ranges) -> Findings:
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
    findings = Findings(db_id=db_id, source_format=source_format)

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
