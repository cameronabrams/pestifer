#!/usr/bin/env python3
"""Regenerate mdplot doc figures from the sweep's archived NAMD logs.

Each example's own config says what its mdplot task plots, so the spec is read from there rather
than restated here -- a figure regenerated under different settings than the build would produce
is worse than a stale one.  Logs are selected the way mdplot itself selects data: only the runs
after the last task that changed the system, since earlier dynamics describe a different system.
"""
import glob, os, re, subprocess, sys, tarfile, tempfile, shutil, yaml

SWEEP = os.path.expanduser('~/devtests/pestifer/examples-sweep-2026-08-10')
REPO = os.path.expanduser('~/Git/pestifer')
PESTIFER = f'{REPO}/.venv/bin/pestifer'
# tasks that change the system, so dynamics before them are not commensurable
SYSTEM_CHANGING = {'psfgen', 'solvate', 'make_membrane_system', 'merge', 'desolvate', 'pdb2pqr'}
LOG_RE = re.compile(r'^(\d+)-(\d+)-(\d+)_(.*)\.log$')


def mdplot_spec(example):
    cfg = glob.glob(f'{REPO}/pestifer/resources/examples/{example}/inputs/*.yaml')
    if not cfg:
        return None, None
    d = yaml.safe_load(open(cfg[0]))
    tasks = [list(t)[0] for t in d['tasks']]
    spec = next((t['mdplot'] for t in d['tasks'] if list(t)[0] == 'mdplot'), None)
    if spec is None:
        return None, None
    last_change = max([i for i, t in enumerate(tasks) if t in SYSTEM_CHANGING], default=-1)
    return spec, last_change


def extract_logs(example, work):
    tb = glob.glob(f'{SWEEP}/example-{example}/*artifacts.tar.gz')
    if not tb:
        return []
    with tarfile.open(tb[0]) as tf:
        # .xst too: cell dimensions live in the xst trajectory, not the log, so extracting
        # only logs yields an empty a_x/b_y/c_z figure rather than an error
        members = [m for m in tf.getmembers() if m.name.endswith(('.log', '.xst'))]
        tf.extractall(work, members=members)
    return sorted(glob.glob(f'{work}/**/*.log', recursive=True))


def namd_logs_after(logs, last_change):
    """NAMD logs from tasks after the last system-changing one, chronological."""
    out = []
    for p in logs:
        m = LOG_RE.match(os.path.basename(p))
        if not m:
            continue
        controller, task_idx, label = int(m.group(1)), int(m.group(2)), m.group(4)
        # Only the top-level controller.  A nested task such as make_membrane_system runs its own
        # sub-pipeline (calibration patches, the bare quilt) under its own controller index, on
        # systems that are not this one; concatenating those with the main line interleaves
        # unrelated series and produces a meaningless curve rather than an error.
        if controller != 0:
            continue
        if task_idx <= last_change:
            continue
        if not any(k in label for k in ('md-', 'density_equilibrate', 'membrane_equilibrate')):
            continue
        out.append((int(m.group(1)), task_idx, int(m.group(3)), p))
    return [p for *_, p in sorted(out)]


def main(examples):
    for ex in examples:
        spec, last_change = mdplot_spec(ex)
        if spec is None:
            print(f'example-{ex}: no mdplot task in config; skipping')
            continue
        work = tempfile.mkdtemp()
        try:
            logs = namd_logs_after(extract_logs(ex, work), last_change)
            if not logs:
                print(f'example-{ex}: no NAMD logs found; skipping')
                continue
            traces = []
            for t in spec.get('timeseries', []):
                traces.append(','.join(t) if isinstance(t, list) else t)
            basename = spec.get('basename', 'solvated')
            cmd = [PESTIFER, 'mdplot', '--basename', basename, '--logs', *logs]
            # a nested entry means "overlay these on one axes", which the CLI spells
            # --timecoseries; flattening them into --timeseries would make one figure per trace
            scalars = [t for t in spec.get('timeseries', []) if not isinstance(t, list)]
            groups = [t for t in spec.get('timeseries', []) if isinstance(t, list)]
            if scalars:
                cmd += ['--timeseries', *scalars]
            for g in groups:
                cmd += ['--timecoseries', *g]
            if spec.get('profiles'):
                cmd += ['--profiles', *spec['profiles']]
            r = subprocess.run(cmd, cwd=work, capture_output=True, text=True, timeout=3600)
            made = sorted(glob.glob(f'{work}/mdplots/*.png'))
            print(f'example-{ex}: {len(logs)} logs -> {[os.path.basename(m) for m in made]}'
                  + ('' if r.returncode == 0 else f'  [EXIT {r.returncode}]'))
            for m in made:
                print(f'    {m}')
            # leave them in place for the caller to copy
            dest = f'{SWEEP}/example-{ex}/regen'
            os.makedirs(dest, exist_ok=True)
            for m in made:
                shutil.copy(m, dest)
        finally:
            shutil.rmtree(work, ignore_errors=True)


main(sys.argv[1:])
