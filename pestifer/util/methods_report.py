# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A drafted Methods section, generated from what a build actually did.

Analogue of ``gmx report-methods``, with more to say: a pestifer build knows not only the final
system but how it was produced -- which entry it started from, what was modelled in, which force
field, what equilibration actually ran, under which software versions, from which seed.

**The register here is deliberately machine-generated.** Short declarative sentences, no
connective flourish, no rationale. That is a design decision, not a limitation: fluent prose reads
as authoritative and invites being pasted unread, while text that visibly came from a tool invites
the editing it needs. The author supplies the scientific voice; this supplies the facts they would
otherwise transcribe by hand and occasionally get wrong.

Three rules hold throughout:

1. Every sentence traces to a fact in the run record. Nothing is inferred from the config where
   the run could have diverged from it -- and the adaptive equilibrations routinely do.
2. Pestifer never describes what it did not do. Production simulation, ensemble choice for the
   science, and analysis are emitted as visible ``\\TODO`` placeholders, never as prose.
3. The output says it is a draft, in a form that has to be deleted deliberately.

See ``docs/design/methods-report.md``.
"""

import logging
import os

from .citations import Citation, by_key, _slug

logger = logging.getLogger(__name__)

DRAFT_BANNER = 'PESTIFER-GENERATED DRAFT -- READ AND EDIT BEFORE USE'


def _fmt_int(n):
    return f'{int(n):,}' if isinstance(n, (int, float)) else str(n)


def _tex_escape(text):
    """Escape the LaTeX specials that appear in pestifer's own vocabulary.

    Task names carry underscores (``density_equilibrate``, ``make_membrane_system``); left raw
    they are a math-mode subscript and the fragment does not compile.
    """
    out = str(text)
    for ch in ('\\', '&', '%', '$', '#', '_', '{', '}'):
        out = out.replace(ch, '\\' + ch)
    return out


def _first(records, *path, default=None):
    """First non-empty value at ``path`` across ``records``."""
    for r in records:
        node = r
        for k in path:
            node = (node or {}).get(k) if isinstance(node, dict) else None
        if node not in (None, '', {}, []):
            return node
    return default


def _agree(records, *path):
    """The shared value at ``path`` if every record agrees, else ``None``.

    Replicas of one system should agree on composition and software; a disagreement is a fact the
    reader needs, not something to average away.
    """
    values = []
    for r in records:
        node = r
        for k in path:
            node = (node or {}).get(k) if isinstance(node, dict) else None
        values.append(node)
    return values[0] if values and all(v == values[0] for v in values) else None



def _box_spread(records):
    """Min and max of each box edge across records, or ``None`` if any lacks a box."""
    edges = []
    for r in records:
        box = ((r.get('system') or {}).get('box') or {})
        if not all(k in box for k in 'abc'):
            return None
        edges.append((box['a'], box['b'], box['c']))
    if not edges:
        return None
    return {'min': tuple(min(e[i] for e in edges) for i in range(3)),
            'max': tuple(max(e[i] for e in edges) for i in range(3))}


def _seeds(records):
    return [r.get('seed') for r in records if r.get('seed') is not None]


def _protocol_summary(record):
    """The equilibration ladder a single run actually performed."""
    lines = []
    for step in record.get('protocol', []):
        ens = str(step.get('ensemble', '')).upper()
        steps = step.get('steps')
        bits = [f"{step.get('task')}"]
        if ens:
            bits.append(f'({ens})')
        if steps:
            bits.append(f'{_fmt_int(steps)} steps')
        if step.get('adaptive'):
            bits.append('run to convergence' if step.get('converged') else 'ran to its step ceiling')
        lines.append(' '.join(bits))
    return lines


def collect_citations(records):
    """Citations owed across every record, deduplicated, in first-seen order.

    Records store citations as plain dicts (that is what JSON allows), so they are rehydrated
    here into :class:`~pestifer.util.citations.Citation` objects to reuse ``bib_entry``.
    """
    seen, out = set(), []
    for r in records:
        for entry in (r.get('citations') or {}).get('entries', []):
            text = entry.get('text', '')
            subject = entry.get('subject', '')
            ident = (subject, text)
            if ident in seen:
                continue
            seen.add(ident)
            key = entry.get('key', '')
            # Prefer the catalog's entry over anything frozen into the record: a bibliography
            # corrected since the build ran should reach the draft.
            catalog = by_key(key) if key else None
            out.append(Citation(subject=subject, text=text, doi=entry.get('doi', ''),
                                reason=entry.get('reason', ''), key=key,
                                bibtex=(catalog.bibtex if catalog else entry.get('bibtex', ''))))
    return out


def render_bib(records):
    """``methods.bib`` covering every citation the runs owe."""
    lines = [f'% {DRAFT_BANNER}',
             '% Generated by pestifer from run-record.json.',
             '% Entries without full bibliographic detail carry a DOI: a reference manager can',
             '% resolve those correctly, whereas guessed metadata would look complete and be wrong.',
             '']
    for c in collect_citations(records):
        lines.append(c.bib_entry())
        lines.append('')
    return '\n'.join(lines)


def _citation_keys(records):
    return [c.key or _slug(c.subject) for c in collect_citations(records)]


def render_tex(records):
    """``methods.tex`` -- a fragment, no preamble, ``\\input``-able into a manuscript."""
    n = len(records)
    L = []
    add = L.append

    add(f'% {DRAFT_BANNER}')
    add('% Generated by pestifer from run-record.json. Every number below is what the build')
    add('% actually did, not what its input requested. Sentences are deliberately terse and')
    add('% machine-shaped; rewrite them in your own voice.')
    add('%')
    add('% Requires \\newcommand{\\TODO}[1]{\\textbf{[TODO: #1]}} or similar in your preamble.')
    add('')
    add('\\subsection*{System preparation \\TODO{retitle or merge into your Methods}}')
    add('')

    version = _agree(records, 'pestifer_version') or _first(records, 'pestifer_version')
    charmmff = _agree(records, 'environment', 'charmmff')

    # --- provenance of the coordinates
    coord = [c for c in collect_citations(records) if c.subject.startswith('PDB ')]
    if coord:
        ids = ', '.join(sorted({c.subject.replace('PDB ', '') for c in coord}))
        cites = ','.join(sorted({c.key or _slug(c.subject) for c in coord}))
        add(f'Starting coordinates were taken from Protein Data Bank entries {ids}~\\cite{{{cites}}}.')
    else:
        add('Starting coordinates \\TODO{state the source; the run record names none}.')

    # --- system composition
    atoms = _agree(records, 'system', 'atoms')
    box = _agree(records, 'system', 'box')
    charge = _agree(records, 'system', 'total_charge')
    if atoms:
        sentence = f'The prepared system contains {_fmt_int(atoms)} atoms'
        if box:
            sentence += (f' in a periodic box of '
                         f'{box.get("a")} $\\times$ {box.get("b")} $\\times$ {box.get("c")} \\AA')
        add(sentence + '.')
        if box is None and n > 1:
            # Replicas equilibrate to slightly different cells.  Report the spread rather than
            # dropping the fact: a reader needs the box, and a range is the honest form of it.
            spread = _box_spread(records)
            if spread:
                add(f'Final box dimensions across replicas ranged from '
                    f'{spread["min"][0]} $\\times$ {spread["min"][1]} $\\times$ {spread["min"][2]} to '
                    f'{spread["max"][0]} $\\times$ {spread["max"][1]} $\\times$ {spread["max"][2]} \\AA.')
        if charge is not None:
            add(f'The net charge of the built system is {charge:g}~e.')
    elif n > 1:
        add('\\TODO{Replica system sizes differ; state them individually or explain why.}')

    # --- force field
    ff_keys = [c.key or _slug(c.subject) for c in collect_citations(records)
               if 'CHARMM' in c.subject or c.subject == 'CGenFF']
    if ff_keys:
        ff = f'~\\cite{{{",".join(ff_keys)}}}'
        release = f' (release {charmmff})' if charmmff else ''
        add(f'The CHARMM36 force field{release} was used{ff}.')

    # --- what actually ran
    add('')
    add('\\subsection*{Equilibration \\TODO{retitle}}')
    add('')
    ladder = _protocol_summary(records[0])
    if ladder:
        add('The following was performed, in order:')
        add('\\begin{enumerate}')
        for item in ladder:
            add(f'  \\item \\texttt{{{_tex_escape(item)}}}')
        add('\\end{enumerate}')
        adaptive = [s for s in records[0].get('protocol', []) if s.get('adaptive')]
        for step in adaptive:
            why = step.get('stopped_because')
            if why:
                add(f'% {step.get("task")}: {why}')
        if adaptive:
            converged = [a for a in adaptive if a.get('converged')]
            ceilinged = [a for a in adaptive if not a.get('converged')]
            if converged:
                add('Adaptive stages ran until their convergence criterion was met rather than '
                    'for a fixed number of steps; the step counts above are what this build '
                    'performed.')
            if ceilinged:
                names = ', '.join(f'\\texttt{{{_tex_escape(a.get("task"))}}}' for a in ceilinged)
                add(f'{names} stopped at its step ceiling without meeting its convergence '
                    f'criterion.')
                add('\\TODO{Report this honestly or re-run to convergence: the system may not '
                    'have settled.}')
        if n > 1:
            add('\\TODO{Step counts are from replica 1; adaptive stages may differ between '
                'replicas. State the range if it matters.}')
    else:
        add('\\TODO{The run record lists no simulation stages.}')

    # --- replicas and reproducibility
    add('')
    add('\\subsection*{Software and reproducibility \\TODO{retitle or move}}')
    add('')
    seeds = _seeds(records)
    if n > 1:
        add(f'{n} independent replicas were prepared from the same input, differing only in '
            f'their random number seeds ({", ".join(str(s) for s in seeds)}).')
    elif seeds:
        add(f'The random number seed was {seeds[0]}.')

    env = _first(records, 'environment') or {}
    exes = env.get('executables') or {}
    named = []
    for name in ('namd3', 'namd3gpu', 'vmd'):
        info = exes.get(name) or {}
        exe_version = info.get('version')
        # A probe that could not answer records 'unknown'; that is meaningful in a provenance
        # log and meaningless in a sentence, so it is dropped rather than published.
        if exe_version and exe_version != 'unknown':
            named.append(exe_version)
    if named:
        add('Simulations and system construction used ' + '; '.join(named) + '.')
    py = env.get('python') or {}
    if version:
        add(f'Systems were prepared with pestifer {version}'
            + (f' under Python {py.get("version")}' if py.get('version') else '') + '.')
    cfg = _agree(records, 'config_file')
    if cfg:
        add(f'The complete specification of the build is the input file \\texttt{{{cfg}}}, '
            f'distributed with this work.')
    add('\\TODO{Cite the archived pestifer release (Zenodo DOI) here.}')

    add('')
    add('\\subsection*{Production simulation}')
    add('')
    add('\\TODO{Pestifer prepared the system; it did not run the production simulation. '
        'Describe the ensemble, length, timestep, thermostat/barostat settings and analysis here.}')
    add('')
    return '\n'.join(L) + '\n'


def render_markdown(records):
    """The same facts as Markdown, for authors not writing LaTeX."""
    n = len(records)
    L = [f'<!-- {DRAFT_BANNER} -->', '', '## System preparation', '']
    coord = [c for c in collect_citations(records) if c.subject.startswith('PDB ')]
    if coord:
        L.append('Starting coordinates: ' +
                 ', '.join(sorted({c.subject.replace('PDB ', '') for c in coord})) + '.')
    atoms = _agree(records, 'system', 'atoms')
    box = _agree(records, 'system', 'box')
    if atoms:
        s = f'System size: {_fmt_int(atoms)} atoms'
        if box:
            s += f'; box {box.get("a")} x {box.get("b")} x {box.get("c")} A'
        L.append(s + '.')
    charmmff = _agree(records, 'environment', 'charmmff')
    if charmmff:
        L.append(f'Force field: CHARMM36, release {charmmff}.')
    L += ['', '## Equilibration', '']
    for item in _protocol_summary(records[0]):
        L.append(f'- {item}')
    seeds = _seeds(records)
    L += ['', '## Reproducibility', '']
    if n > 1:
        L.append(f'{n} replicas, seeds {", ".join(str(s) for s in seeds)}.')
    elif seeds:
        L.append(f'Seed: {seeds[0]}.')
    L.append('TODO: describe the production simulation; pestifer did not run it.')
    L.append('')
    return '\n'.join(L)


def write_report(records, outdir='.', formats=('tex', 'bib')):
    """Write the requested renderings; returns the paths written."""
    renderers = {'tex': ('methods.tex', render_tex),
                 'bib': ('methods.bib', render_bib),
                 'md': ('methods.md', render_markdown)}
    written = []
    for fmt in formats:
        if fmt not in renderers:
            logger.warning(f'unknown methods-report format {fmt!r}; skipping')
            continue
        name, render = renderers[fmt]
        path = os.path.join(outdir, name)
        with open(path, 'w') as f:
            f.write(render(records))
        logger.info(f'wrote {path}')
        written.append(path)
    return written
