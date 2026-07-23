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
    """What a single chain is: its dominant segment type, size, molecule name, and resname sample."""
    chain: str
    segtype: str
    n_residues: int
    resnames: list = field(default_factory=list)
    auth: str = None          # mmCIF author chain id (label chain is `chain`); None for PDB
    molecule: str = None      # molecule name from COMPND / mmCIF entity description, if available
    attached_to: str = None   # for a glycan chain: the protein chain it is glycosidically linked to

    def describe(self) -> str:
        sample = ', '.join(self.resnames[:4]) + ('...' if len(self.resnames) > 4 else '')
        if self.segtype == 'water':
            base = 'water'
        elif self.segtype == 'ion':
            base = f'ion ({sample})'
        elif self.segtype == 'glycan':
            base = f'glycan ({sample})'
        elif self.segtype in ('protein', 'nucleicacid'):
            name = 'protein' if self.segtype == 'protein' else 'nucleic acid'
            base = f'{name} ({self.n_residues} residues)'
        else:
            base = f'{self.segtype} ({sample})'
        if self.molecule and self.segtype in ('protein', 'nucleicacid', 'other'):
            return f'{base} — {self.molecule}'
        return base


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
                att = f' [glycan on chain {ch.attached_to}]' if ch.attached_to else ''
                c(f'  {ch.chain}{auth}: {ch.describe()}{att}')
            c('To OMIT chains, add under `source:` (also exclude any glycan chain attached to an')
            c('omitted protein chain) -')
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
            c('To replace a long/disordered loop with a short built stub (not the full sequence),')
            c('add a `mods:` block (a SIBLING of `source:`) -')
            c('  mods:')
            c(f'    substitutions: [{gaps[0].chain}:{gaps[0].start}-{gaps[0].end},GGG]')

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


def make_text_prompter():
    """A default free-text prompter reading from stdin; EOF/blank returns the shown default."""
    def prompt_text(question, default=''):
        try:
            resp = input(f'{question} [{default}] ').strip()
        except EOFError:
            return default
        return resp or default
    return prompt_text


def interactive_select(findings: 'Findings', ask=None, say=None, prompt_text=None) -> dict:
    """
    Walk the user through the findings, prompting per item, and return the chosen mods as
    ``{'source', 'sequence', 'mods', 'add_ligate'}`` -- the *active* (uncommented) blocks to inject
    into the generated psfgen task.

    ``ask(question, default) -> bool``, ``prompt_text(question, default) -> str``, and
    ``say(message)`` are injectable (defaults: stdin prompters and ``print``) so the walkthrough is
    testable without a TTY.
    """
    ask = ask or make_prompter()
    prompt_text = prompt_text or make_text_prompter()
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

    # chain inclusion (default INCLUDE). A glycan chain glycosidically attached to a protein chain
    # follows it: it is not queried separately, and is omitted iff its protein chain is omitted --
    # and when a chain is omitted, its missing loops/tails and mutations are not queried either.
    included = {c.chain for c in findings.chains} if findings.chains else None   # None = all
    if included is not None and len(findings.chains) > 1:
        say('\nChains present:')
        attached = {c.chain: c.attached_to for c in findings.chains if c.attached_to}
        for ch in findings.chains:
            if ch.attached_to:                       # glycan chain -> follows its protein chain
                continue
            auth = f' (author {ch.auth})' if ch.auth and ch.auth != ch.chain else ''
            glys = [g for g, prot in attached.items() if prot == ch.chain]
            extra = f' (+ attached glycan chain(s) {", ".join(glys)})' if glys else ''
            if not ask(f'  Include chain {ch.chain}{auth} [{ch.describe()}]{extra}?', True):
                included.discard(ch.chain)
                for g in glys:
                    included.discard(g)
        omitted = sorted({c.chain for c in findings.chains} - included)
        if omitted:
            source['exclude'] = [f"chainID in [{', '.join(repr(c) for c in omitted)}]"]

    def _keep(ch):
        return included is None or ch in included

    tails = {'n': [], 'c': []}
    term_runs = [r for r in findings.missing_runs if r.kind in ('N', 'C') and _keep(r.chain)]
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

    gaps = [r for r in findings.interior_gaps() if _keep(r.chain)]
    add_ligate = False
    if gaps:
        say('\nInterior missing loops (built and closed automatically; a `ligate` task is added):')
        stubs = []
        for r in gaps:
            if ask(f'  Build interior loop {r.chain} {r.start}-{r.end} ({r.length} res) in full?', True):
                continue
            # An interior loop cannot currently be left genuinely unbuilt (no capped-break path);
            # the supported alternative is a short built stub in its place.
            if ask(f'    Replace it with a short built stub sequence instead? (No = keep full)', True):
                seq = ''.join(prompt_text(f'      Stub sequence (one-letter codes)?', 'GGG').split()).upper()
                stubs.append(f'{r.chain}:{r.start}-{r.end},{seq or "GGG"}')
        if stubs:
            mods.setdefault('substitutions', []).extend(stubs)
        add_ligate = True     # interior loops (full or stubbed) are built and need closing

    def _built(chain, resid):
        # a mutation must target a residue that will exist in the built structure: a resolved
        # residue, an interior loop (always built), or a terminal tail the user chose to build.
        for r in findings.missing_runs:
            if r.chain == chain and r.start <= resid <= r.end:
                if r.kind == 'interior':
                    return True
                return chain in (tails['n'] if r.kind == 'N' else tails['c'])
        return True   # not in any missing run -> resolved

    revertible = [m for m in findings.mutations if m.dbres and _keep(m.chain) and _built(m.chain, m.resid)]
    if revertible:
        say('\nEngineered mutations / conflicts (structure differs from the sequence database):')
        muts = []
        for m in revertible:
            if ask(f'  Revert {m.chain}:{m.resname}{m.resid} -> db {m.dbres} [{m.typekey}]?', False):
                muts.append(m.revert_shortcode)
        if muts:
            mods['mutations'] = muts

    excisions = [e for e in findings.excisions if _keep(e.chain)]
    if excisions:
        say('\nCloning artifacts / expression tags (extra residues not in the native protein):')
        dels = []
        for e in excisions:
            span = f'{e.start}' if e.start == e.end else f'{e.start}-{e.end}'
            if ask(f'  Excise {e.typekey} in chain {e.chain} ({span})?', False):
                dels.append(e.delete_shortcode)
        if dels:
            mods['deletions'] = dels

    return {'source': source, 'sequence': sequence, 'mods': mods, 'add_ligate': add_ligate}


def interactive_pipeline(db_id: str = 'system', ask=None, say=None) -> list:
    """
    Walk the user through the post-psfgen build pipeline and return the chosen task dicts (to append
    after ``fetch``/``psfgen``/``ligate``): vacuum minimization, vacuum MD, solvation and its
    solvated minimization / NVT warm-up / self-terminating NPT density equilibration / optional
    production, and a standard ``terminate`` (package) task. Each stage is a yes/no prompt with a
    sensible default; ``ask`` is injectable for testing. nsteps are conservative defaults the user
    can edit in the YAML.
    """
    ask = ask or make_prompter()
    say = say or (lambda m: print(m))
    tasks = []
    say('\nPost-psfgen build pipeline (each stage optional; edit nsteps in the YAML as needed):')
    if ask('  Vacuum minimization?', True):
        tasks.append({'md': {'ensemble': 'minimize'}})
    if ask('  Vacuum MD (NVT) equilibration?', False):
        tasks.append({'md': {'ensemble': 'NVT', 'nsteps': 1000}})
    if ask('  Solvate (water box + neutralizing ions)?', True):
        tasks.append({'solvate': None})
        if ask('    Solvated minimization?', True):
            tasks.append({'md': {'ensemble': 'minimize'}})
        if ask('    NVT warm-up before density equilibration?', True):
            tasks.append({'md': {'ensemble': 'NVT', 'nsteps': 1000}})
        # Self-terminating NPT density equilibration (replaces the fixed NPT ladder): runs
        # stability-bounded short restarts and stops when the box density has converged. Writes
        # its own convergence report + density-vs-time plot, so no separate mdplot stage is needed.
        if ask('    NPT density equilibration (self-terminating)?', True):
            tasks.append({'density_equilibrate': None})
        if ask('    Longer NPT production run?', False):
            tasks.append({'md': {'ensemble': 'NPT', 'nsteps': 100000}})
    if ask('  Terminate task (package the built system)?', True):
        tasks.append({'terminate': {'basename': f'my_{db_id.lower()}',
                                    'package': {'basename': f'prod_{db_id.lower()}',
                                                'namd': {'ensemble': 'NPT'}}}})
    return tasks


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


def _chain_molecule_names(parsed) -> dict:
    """
    ``{author-chain-id: molecule name}`` from the ``COMPND`` record.

    pidibble normalizes both PDB ``COMPND`` (``MOLECULE``/``CHAIN`` tokens) and mmCIF entity data
    (``_entity.pdbx_description`` joined to ``_entity_poly.pdbx_strand_id``) into this record, so one
    reader handles both. The chain ids here are *author* ids; the caller maps them onto its label
    ids for mmCIF via each chain's ``auth``.
    """
    out = {}
    comp = parsed.get('COMPND')
    if comp is None:
        return out
    compound = getattr(comp, 'compound', None)
    if isinstance(compound, list):                      # PDB: 'KEY: VALUE' tokens grouped by MOL_ID
        cur = None
        for tok in compound:
            if ':' not in tok:
                continue
            key, val = tok.split(':', 1)
            key = key.strip().upper()
            val = val.strip().rstrip(';')
            if key == 'MOLECULE':
                cur = val
            elif key == 'CHAIN' and cur:
                for ch in val.split(','):
                    out.setdefault(ch.strip(), cur)
        return out
    for r in (comp if isinstance(comp, list) else [comp]):   # mmCIF: per-entity records
        mol = getattr(r, 'molName', None)
        chs = getattr(r, 'chains', None)
        if mol and chs:
            for ch in (chs if isinstance(chs, list) else [chs]):
                out.setdefault(str(ch), mol)
    return out


def _glycan_attachment(parsed) -> dict:
    """``{glycan-chain: protein-chain}`` from cross-chain protein<->glycan ``LINK`` records -- so a
    separate glycan chain can follow the protein chain it is glycosidically attached to."""
    from .labels import Labels
    out = {}
    link = parsed.get('LINK')
    if link is None:
        return out
    recs = getattr(link, 'data', link)
    recs = recs if isinstance(recs, list) else [recs]

    def seg(rn):
        return Labels.segtype_of_resname.get(rn, 'other')

    for r in recs:
        try:
            c1, rn1 = r.residue1.chainID, r.residue1.resName
            c2, rn2 = r.residue2.chainID, r.residue2.resName
        except AttributeError:
            continue
        if c1 == c2:
            continue
        s1, s2 = seg(rn1), seg(rn2)
        if s1 == 'protein' and s2 == 'glycan':
            out.setdefault(c2, c1)
        elif s2 == 'protein' and s1 == 'glycan':
            out.setdefault(c1, c2)
    return out


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

    # Resolved *polymer* residues per chain (protein/nucleic only). A missing residue is a polymer
    # residue, so its interior/terminal classification must be against the polymer span -- NOT the
    # whole chain, whose range can be inflated by in-chain glycans/ions at high resids (e.g. a
    # C-terminal protein gap would otherwise look interior because NAGs sit past it in the same chain).
    resolved = defaultdict(set)
    for a in atoms.data:
        if Labels.segtype_of_resname.get(a.resname, 'other') in ('protein', 'nucleicacid'):
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
    molnames = _chain_molecule_names(parsed)
    attach = _glycan_attachment(parsed)
    chains = []
    for ch, d in chain_data.items():
        segs = Counter(Labels.segtype_of_resname.get(rn, 'other') for rn in d['resnames'])
        mol = molnames.get(ch) or (molnames.get(d['auth']) if d['auth'] else None)
        chains.append(ChainIdentity(ch, segs.most_common(1)[0][0], len(d['residues']),
                                    d['resnames'], d['auth'], mol, attach.get(ch)))
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
