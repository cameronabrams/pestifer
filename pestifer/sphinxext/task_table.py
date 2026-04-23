# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Sphinx directive ``task-table`` that reads a pestifer YAML config file and
renders a compact pipeline-summary table in the documentation.

Consecutive ``md`` tasks are collapsed into a single row showing the ensemble
sequence and total step counts, keeping the table concise even for workflows
with many equilibration phases.

Usage in rst::

    .. task-table:: ../../../../pestifer/resources/examples/ex01/inputs/bpti1.yaml

The path is resolved relative to the rst source file (same convention as
``.. literalinclude::``).
"""
import os
import yaml
from docutils import nodes
from docutils.statemachine import ViewList
from sphinx.util.docutils import SphinxDirective


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def _fetch_details(specs: dict) -> str:
    source = (specs or {}).get('source', 'rcsb')
    sid = (specs or {}).get('sourceID', '?')
    if source == 'alphafold':
        return f'AlphaFold model, UniProt {sid}'
    url = f'https://www.rcsb.org/structure/{sid}'
    return f'`PDB {sid.upper()} <{url}>`_'


def _psfgen_details(specs: dict) -> str:
    specs = specs or {}
    src = specs.get('source', {}) or {}
    if isinstance(src, list):
        src = src[0] if src else {}
    mods = specs.get('mods', {}) or {}
    tags = []
    # exclusions live under source
    if src.get('exclude'):
        tags.append('exclusions')
    # modifications live under mods (mutations, bonds, grafts, etc.)
    if mods.get('mutations'):
        tags.append('mutations')
    if mods.get('ssbonds'):
        tags.append('new disulfides')
    if mods.get('ssbondsdelete'):
        tags.append('delete disulfides')
    if mods.get('substitutions'):
        tags.append('loop substitutions')
    grafts = mods.get('grafts') or []
    if grafts:
        tags.append(f'{len(grafts)} glycan graft(s)')
    if mods.get('loops') or src.get('loops'):
        tags.append('loop modeling')
    if mods.get('crotations'):
        tags.append('crotations')
    return ', '.join(tags) if tags else 'standard build'


def _md_group_details(md_specs_list: list) -> str:
    """Summarize a run of consecutive md tasks into one readable string."""
    parts = []
    i = 0
    while i < len(md_specs_list):
        specs = md_specs_list[i] or {}
        ensemble = specs.get('ensemble', 'NVT')
        if ensemble.lower() == 'minimize':
            parts.append('minimize')
            i += 1
        else:
            # collect consecutive phases of the same ensemble
            phases = []
            ens = ensemble
            while i < len(md_specs_list):
                s = md_specs_list[i] or {}
                e = s.get('ensemble', 'NVT')
                if e.lower() != ens.lower():
                    break
                phases.append(s.get('nsteps', 2000))
                i += 1
            total = sum(phases)
            n = len(phases)
            if n == 1:
                parts.append(f'{ens} ({total:,} steps)')
            else:
                parts.append(f'{ens} ({total:,} steps, {n} phases)')
    return ' → '.join(parts)


def _manipulate_details(specs: dict) -> str:
    mods = list(((specs or {}).get('mods') or {}).keys())
    return ', '.join(f'``{m}``' for m in mods) if mods else '—'


def _terminate_details(specs: dict) -> str:
    specs = specs or {}
    parts = []
    bn = specs.get('basename')
    if bn:
        parts.append(f'basename: ``{bn}``')
    pkg_bn = (specs.get('package') or {}).get('basename')
    if pkg_bn:
        parts.append(f'package: ``{pkg_bn}``')
    return '; '.join(parts) if parts else '—'


def _membrane_details(specs: dict) -> str:
    bilayer = (specs or {}).get('bilayer', {}) or {}
    comp = bilayer.get('composition', {}) or {}
    lipids = set()
    for leaflet in [comp.get('upper_leaflet', []), comp.get('lower_leaflet', [])]:
        for item in (leaflet or []):
            if isinstance(item, dict):
                lipids.add(item.get('name', '?'))
    return ', '.join(sorted(lipids)) if lipids else '—'


def _single_task_details(name: str, specs: dict) -> str:
    if name == 'fetch':
        return _fetch_details(specs)
    if name == 'psfgen':
        return _psfgen_details(specs)
    if name == 'solvate':
        specs = specs or {}
        conc = specs.get('salt_con')
        cation = specs.get('cation')
        anion = specs.get('anion')
        if conc and cation and anion:
            return f'water box + {conc} M {cation}/{anion}'
        elif conc:
            return f'water box + {conc} M salt'
        return 'water box'
    if name == 'ligate':
        steer = (specs or {}).get('steer', {}) or {}
        nsteps = steer.get('nsteps')
        if nsteps:
            return f'ligate chain breaks (steered MD, {nsteps:,} steps)'
        return 'ligate chain breaks'
    if name == 'cleave':
        sites = (specs or {}).get('sites', []) or []
        n = len(sites)
        return f'proteolytic cleavage at {n} site(s)' if n else 'proteolytic cleavage'
    if name == 'pdb2pqr':
        pH = (specs or {}).get('pH')
        return f'protonate at pH {pH}' if pH is not None else 'protonate (pdb2pqr)'
    if name == 'ring_check':
        return 'check for ring-threading defects'
    if name == 'manipulate':
        return _manipulate_details(specs)
    if name == 'terminate':
        return _terminate_details(specs)
    if name == 'continuation':
        psf = (specs or {}).get('psf', '')
        pdb = (specs or {}).get('pdb', '')
        return f'from {psf}' if psf else '—'
    if name == 'validate':
        n = len((specs or {}).get('tests', []))
        return f'{n} test(s)'
    if name == 'mdplot':
        output_dir = (specs or {}).get('output_dir', 'mdplots')
        return f'equilibration time-series plots → ``{output_dir}/``'
    if name == 'make_membrane_system':
        return _membrane_details(specs)
    if name == 'merge':
        return 'merge pipeline systems'
    return '—'


# ---------------------------------------------------------------------------
# Task grouper
# ---------------------------------------------------------------------------

def _task_rows(tasks: list) -> list[tuple[str, str]]:
    """Return list of (task_name, details) rows with md tasks collapsed."""
    rows = []
    i = 0
    while i < len(tasks):
        td = tasks[i]
        name = list(td.keys())[0]
        specs = td[name]
        if name == 'md':
            group = []
            while i < len(tasks):
                t = tasks[i]
                n = list(t.keys())[0]
                if n != 'md':
                    break
                group.append(t['md'] or {})
                i += 1
            rows.append(('md', _md_group_details(group)))
        else:
            rows.append((name, _single_task_details(name, specs)))
            i += 1
    return rows


# ---------------------------------------------------------------------------
# RST generator
# ---------------------------------------------------------------------------

def _table_rst_lines(tasks: list) -> list[str]:
    rows = _task_rows(tasks)
    lines = [
        '.. list-table:: Pipeline task summary',
        '   :header-rows: 1',
        '   :widths: 5 20 75',
        '',
        '   * - Step',
        '     - Task',
        '     - Details',
        '',
    ]
    for step, (name, details) in enumerate(rows, 1):
        lines.append(f'   * - {step}')
        lines.append(f'     - ``{name}``')
        lines.append(f'     - {details}')
        lines.append('')
    return lines


# ---------------------------------------------------------------------------
# Sphinx directive
# ---------------------------------------------------------------------------

class TaskTableDirective(SphinxDirective):
    """Renders a pipeline-summary table from a pestifer YAML config file."""

    required_arguments = 1
    optional_arguments = 0
    has_content = False

    def run(self):
        yaml_rel = self.arguments[0]
        rst_dir = os.path.dirname(
            os.path.join(self.env.srcdir, self.env.docname)
        )
        yaml_abs = os.path.normpath(os.path.join(rst_dir, yaml_rel))

        try:
            with open(yaml_abs) as f:
                config = yaml.safe_load(f)
        except FileNotFoundError:
            msg = self.reporter.error(
                f'task-table: YAML file not found: {yaml_abs}',
                line=self.lineno,
            )
            return [msg]

        tasks = config.get('tasks', [])
        rst_lines = _table_rst_lines(tasks)

        container = nodes.section()
        vl = ViewList(rst_lines, source=yaml_abs)
        self.state.nested_parse(vl, self.content_offset, container)
        return container.children
