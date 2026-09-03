# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Offline reconstruction of a *complete* NAMD pressure profile, reciprocal-space term included.

Why this module exists
----------------------

A NAMD run with ``pressureProfile on`` reports "kinetic, bonded and nonbonded (but not reciprocal
space) contributions".  The missing PME reciprocal term is not a small correction -- it can flip
the sign of the integrated profile -- and it cannot be recovered by simply also switching on
``pressureProfileEwald``, because NAMD makes the two mutually exclusive (``SimParameters.C``, in
3.0.3 at line 6699)::

    // if Ewald is on, only calculate Ewald
    if (pressureProfileEwaldOn)
      pressureProfileOn = 0;

So the complete profile needs **two runs, summed frame by frame**.  Running the Ewald half inline
is impractical: ``ComputeEwald`` is registered as an ordinary compute gated only on
``doFullElectrostatics``, never on ``pressureProfileFreq``, so it evaluates every full-electrostatics
step no matter how rarely profiles are sampled -- measured at 14-28x the per-step cost.

Replaying a saved trajectory instead evaluates once per *sampled frame*, which divides that penalty
by the sampling stride: a stage sampled every 100 steps costs ~100x less to reconstruct offline
than to have computed inline.  That is what this module does.

Validation
----------

On a 699-atom TIP3P box the slab-averaged trace/3 of the reconstructed profile agrees with NAMD's
own reported ``PRESSURE`` to ~1 bar, where the real-space half alone is off by ~213 bar.
:func:`validate_against_total_pressure` re-runs that check on any replay.

Frame identity
--------------

Every replayed frame is evaluated at timestep 0, so the ``PRESSUREPROFILE:`` records all carry the
same step number and cannot be matched up by it.  The two passes are aligned **positionally**, and
:func:`combine` refuses to combine passes whose record counts or slab counts differ rather than
silently pairing the wrong frames.
"""
from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess

from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

#: Directives dropped when a production config is turned into a replay config.  Matching is on the
#: first whitespace-delimited token, case-insensitively.
#:
#: Three kinds of thing go here.  *Run control and output* because the replay supplies its own and
#: must not overwrite the trajectory it is reading.  *Thermostats and barostats* because nothing is
#: integrated -- and NAMD special-cases the pressure profile under constant pressure.  *Anything
#: naming the pressure profile* because this module sets all of it.  Everything else -- force
#: field, cutoffs, PME, restraints, colvars, exclusions -- passes through verbatim, because the
#: profile must be computed with the same forces the trajectory was generated with.
DROP_DIRECTIVES = frozenset({
    # run control / restart
    'run', 'minimize', 'numsteps', 'firsttimestep', 'stepspercycle',
    'outputname', 'restartname', 'restartfreq', 'restartsave', 'binaryrestart',
    'reinitvels', 'rescalevels', 'reassignfreq', 'reassigntemp', 'reassignincr', 'reassignhold',
    # trajectory output: the replay must not write over what it is reading
    'dcdfreq', 'dcdfile', 'dcdunitcell', 'veldcdfreq', 'veldcdfile',
    'forcedcdfreq', 'forcedcdfile', 'xstfreq', 'xstfile',
    # initial velocities: replay supplies `temperature`
    'velocities', 'binvelocities', 'temperature',
    # thermostats
    'langevin', 'langevintemp', 'langevindamping', 'langevinhydrogen',
    'langevinfile', 'langevincol', 'loweandersen', 'loweandersentemp',
    'loweandersenrate', 'loweandersencutoff', 'stochrescale', 'stochrescaletemp',
    'stochrescaleperiod', 'stochrescaleheat', 'tcouple', 'tcoupletemp',
    # barostats: nothing is integrated, and the profile is special-cased under constant pressure
    'langevinpiston', 'langevinpistontarget', 'langevinpistonperiod',
    'langevinpistondecay', 'langevinpistontemp',
    'useflexiblecell', 'useconstantarea', 'useconstantratio', 'usegrouppressure',
    'berendsenpressure', 'berendsenpressuretarget', 'berendsenpressurefreq',
    'berendsenpressurecompressibility', 'berendsenpressurerelaxationtime',
    'montecarlopressure', 'montecarlopressurefreq', 'montecarlotemp',
    'surfacetensiontarget', 'strainrate',
    # GPU-resident integration does not support the pressure profile or coorfile replay
    'cudasoaintegrate',
    # this module owns all of these
    'pressureprofile', 'pressureprofilefreq', 'pressureprofileslabs',
    'pressureprofileatomtypes', 'pressureprofileatomtypesfile',
    'pressureprofileatomtypescol', 'pressureprofileewald',
    'pressureprofileewaldx', 'pressureprofileewaldy', 'pressureprofileewaldz',
    # output cadence is set by the replay so every frame reports
    'outputenergies', 'outputtiming', 'outputpressure', 'outputmomenta',
})

#: Tcl commands that would run dynamics or otherwise drive the engine.  A production config often
#: ends in a bare ``run 500000``; a scripted one may loop.  Lines starting with these are dropped
#: along with :data:`DROP_DIRECTIVES` so that only the replay loop drives NAMD.
_TCL_DRIVERS = frozenset({'run', 'minimize', 'move', 'moveallby', 'coorfile', 'reinitatoms'})

_PROFILE_RE = re.compile(r'^PRESSUREPROFILE:\s+(\S+)\s+(.*)$')
_SLAB_THICKNESS_RE = re.compile(r'SLAB THICKNESS:\s*([-\d.eE+]+)')
_FRAMES_RE = re.compile(r'PPREPLAY_FRAMES\s+(\d+)')


@dataclass
class ProfilePass:
    """One replay pass: the per-frame profiles NAMD emitted, plus what it took to get them."""

    mode: str
    """``'real'`` (kinetic + bonded + real-space nonbonded) or ``'ewald'`` (reciprocal only)."""
    profiles: np.ndarray
    """``(nframes, nslabs, 3)`` array of per-slab ``(Pxx, Pyy, Pzz)`` in bar."""
    slab_thickness: float | None = None
    """Angstroms, as NAMD reported it at startup; ``None`` if the log did not say."""
    frames_read: int | None = None
    """Frames the replay loop actually read, from the loop's own counter."""
    log_path: str = ''
    config_path: str = ''
    total_pressure: np.ndarray | None = None
    """NAMD's own ``PRESSURE`` column, one value per frame, when the log carried ENERGY lines."""

    @property
    def nframes(self) -> int:
        return int(self.profiles.shape[0])

    @property
    def nslabs(self) -> int:
        return int(self.profiles.shape[1])


@dataclass
class ReplayResult:
    """The reconstructed profile and the two passes it came from."""

    real: ProfilePass
    ewald: ProfilePass
    total: np.ndarray = field(repr=False)
    """``(nframes, nslabs, 3)``: the elementwise sum, which is the complete profile."""

    @property
    def nframes(self) -> int:
        return int(self.total.shape[0])

    @property
    def nslabs(self) -> int:
        return int(self.total.shape[1])

    def slab_centers(self) -> np.ndarray:
        """Slab-center z coordinates in Angstroms, or slab indices if NAMD gave no thickness."""
        th = self.real.slab_thickness
        n = self.nslabs
        if th is None:
            logger.warning('no slab thickness in the NAMD log; z axis will be slab index')
            return np.arange(n, dtype=float)
        return (np.arange(n, dtype=float) + 0.5) * th

    def mean_profile(self) -> np.ndarray:
        """``(nslabs, 3)``: the complete profile averaged over frames."""
        return self.total.mean(axis=0)

    def sem_profile(self) -> np.ndarray:
        """``(nslabs, 3)``: standard error of the mean, per slab and component."""
        n = self.nframes
        if n < 2:
            return np.zeros_like(self.mean_profile())
        return self.total.std(axis=0, ddof=1) / np.sqrt(n)


def read_config_lines(path: str) -> list[str]:
    with open(path) as f:
        return f.read().splitlines()


def strip_for_replay(lines: list[str]) -> tuple[list[str], list[str]]:
    """Split a production NAMD config into (lines kept, lines dropped).

    Kept lines are emitted verbatim, so the replay sees exactly the force field, cutoffs, PME
    settings and restraints the trajectory was generated with.  Returning the dropped lines rather
    than discarding them silently lets the caller log what it removed -- a config that turns out to
    have driven dynamics from a Tcl loop is worth seeing.
    """
    kept, dropped = [], []
    for raw in lines:
        stripped = raw.strip()
        if not stripped or stripped.startswith('#'):
            kept.append(raw)
            continue
        first = stripped.split()[0].lower()
        if first in DROP_DIRECTIVES or first in _TCL_DRIVERS:
            dropped.append(raw)
        else:
            kept.append(raw)
    return kept, dropped


def build_replay_config(base_lines: list[str], dcd: str, mode: str, *, slabs: int,
                        outname: str, basedir: str = '',
                        ewald_grid: tuple[int, int, int] = (10, 10, 10),
                        stride: int = 1, temperature: float = 300.0) -> str:
    """Render the text of a replay config.

    ``mode`` is ``'real'`` for the ordinary profile or ``'ewald'`` for the reciprocal-space one.
    Both passes get the same slab count; combining passes gridded differently would be meaningless,
    so this is not separately settable.
    """
    if mode not in ('real', 'ewald'):
        raise ValueError(f"mode must be 'real' or 'ewald', not {mode!r}")
    kept, _ = strip_for_replay(base_lines)
    out = []
    if basedir:
        # NAMD chdirs to the config file's own directory (mainfunc.C), which would break every
        # relative path the production config names.  `cwd` puts it back in the directory those
        # paths were written against, so the passthrough lines resolve as they did in production.
        # It must precede them.  Paths this module emits are absolute and so are unaffected.
        out += [f'cwd {os.path.abspath(basedir)}', '']
    out += kept
    out += [
        '',
        '# ---- appended by pestifer pressure-profile-ewald ----',
        '# Nothing is integrated: each frame is read from the trajectory and evaluated in place.',
        f'temperature             {temperature}',
        # NAMD requires outputName even when nothing is integrated.  The production one was
        # dropped so that a replay can never overwrite the run it is reading.
        f'outputName              {outname}',
        'firsttimestep           0',
        'numsteps                0',
        'outputEnergies          1',
        'CUDASOAintegrate        off',
        '',
        f'pressureProfile         on',
        f'pressureProfileSlabs    {slabs}',
        f'pressureProfileFreq     1',
    ]
    if mode == 'ewald':
        out += [
            '# NAMD sets pressureProfileOn = 0 when this is on, so this pass reports the',
            '# reciprocal-space contribution ALONE.  It is summed with the real-space pass.',
            'pressureProfileEwald    on',
            f'pressureProfileEwaldX   {ewald_grid[0]}',
            f'pressureProfileEwaldY   {ewald_grid[1]}',
            f'pressureProfileEwaldZ   {ewald_grid[2]}',
        ]
    skip = ''
    if stride > 1:
        skip = f'\n    for {{set i 1}} {{$i < {stride}}} {{incr i}} {{ coorfile skip }}'
    out += [
        '',
        f'coorfile open dcd {dcd}',
        'set __ppframe 0',
        'while { ![coorfile read] } {',
        '    incr __ppframe',
        '    run 0' + skip,
        '}',
        'coorfile close',
        'print "PPREPLAY_FRAMES $__ppframe"',
        '',
    ]
    return '\n'.join(out) + '\n'


def parse_profiles(log_path: str, mode: str) -> ProfilePass:
    """Read every ``PRESSUREPROFILE:`` record out of a replay log, in order.

    Records are kept in file order deliberately: at ``run 0`` every frame is evaluated at timestep
    0, so the step column is a constant and carries no frame identity.
    """
    rows: list[list[float]] = []
    pressures: list[float] = []
    slab_thickness = None
    frames_read = None
    etitle: list[str] | None = None
    with open(log_path) as f:
        for line in f:
            m = _PROFILE_RE.match(line)
            if m:
                rows.append([float(x) for x in m.group(2).split()])
                continue
            if slab_thickness is None:
                m = _SLAB_THICKNESS_RE.search(line)
                if m:
                    slab_thickness = float(m.group(1))
                    continue
            if frames_read is None:
                m = _FRAMES_RE.search(line)
                if m:
                    frames_read = int(m.group(1))
                    continue
            if line.startswith('ETITLE:'):
                etitle = line.split()[1:]
            elif line.startswith('ENERGY:') and etitle and 'PRESSURE' in etitle:
                pressures.append(float(line.split()[1:][etitle.index('PRESSURE')]))
    if not rows:
        raise RuntimeError(
            f'{log_path}: no PRESSUREPROFILE records. The replay produced no profile -- check the '
            'log for a NAMD error, and that the trajectory atom count matches the PSF.')
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise RuntimeError(f'{log_path}: ragged PRESSUREPROFILE records (widths {sorted(widths)})')
    w = widths.pop()
    if w % 3:
        raise RuntimeError(f'{log_path}: {w} values per record is not a multiple of 3 (xx,yy,zz)')
    arr = np.asarray(rows, dtype=float).reshape(len(rows), w // 3, 3)
    return ProfilePass(mode=mode, profiles=arr, slab_thickness=slab_thickness,
                       frames_read=frames_read, log_path=log_path,
                       total_pressure=np.asarray(pressures) if pressures else None)


def combine(real: ProfilePass, ewald: ProfilePass) -> np.ndarray:
    """Elementwise sum of the two passes, which is the complete profile.

    Refuses mismatched passes.  The frames are aligned by position and nothing in the records
    themselves would reveal a misalignment, so a shape disagreement is the only chance to notice
    that the two passes did not see the same trajectory.
    """
    if real.nframes != ewald.nframes:
        raise RuntimeError(
            f'pass frame counts differ (real {real.nframes}, ewald {ewald.nframes}); the two '
            'passes did not read the same frames, and are aligned by position only')
    if real.nslabs != ewald.nslabs:
        raise RuntimeError(
            f'pass slab counts differ (real {real.nslabs}, ewald {ewald.nslabs}); summing '
            'different slab griddings is meaningless')
    return real.profiles + ewald.profiles


def validate_against_total_pressure(result: ReplayResult) -> dict | None:
    """Compare the slab-averaged trace/3 of each profile with NAMD's own ``PRESSURE``.

    The real-space pass reports the total pressure of a *complete* force evaluation in its ENERGY
    lines even though its profile omits the reciprocal term, so this is an independent check that
    the two-pass sum is right.  Returns ``None`` when the log carried no usable ENERGY lines.
    """
    p = result.real.total_pressure
    if p is None or len(p) < result.nframes:
        return None
    p = np.asarray(p[-result.nframes:], dtype=float)
    real_mean = result.real.profiles.mean(axis=1).sum(axis=1) / 3.0
    total_mean = result.total.mean(axis=1).sum(axis=1) / 3.0
    return {
        'real_only_deviation': float(np.abs(real_mean - p).mean()),
        'reconstructed_deviation': float(np.abs(total_mean - p).mean()),
        'namd_pressure_mean': float(p.mean()),
    }


def run_namd(config_path: str, log_path: str, cwd: str, namd: str = 'namd3',
             nprocs: int = 1) -> None:
    """Run one replay pass, writing NAMD's output to ``log_path``.

    ``cwd`` is the directory the *production* config lives in.  NAMD chdirs to the config file's
    own directory on startup, so the emitted config carries a `cwd` directive that puts it back
    here; starting the process here too keeps the two consistent and makes a relative config path
    on the command line resolve.
    """
    exe = shutil.which(namd)
    if exe is None:
        raise FileNotFoundError(f'NAMD executable {namd!r} not found on PATH')
    cmd = [exe]
    if nprocs > 1:
        cmd.append(f'+p{nprocs}')
    cmd.append(os.path.relpath(os.path.abspath(config_path), os.path.abspath(cwd)))
    logger.info(f'running {" ".join(cmd)} (in {cwd}) > {log_path}')
    with open(log_path, 'w') as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=cwd)
    if proc.returncode != 0:
        raise RuntimeError(
            f'NAMD exited {proc.returncode} on the {os.path.basename(config_path)} pass; '
            f'see {log_path}')


def replay(namd_config: str, dcd: str, *, workdir: str = '.', slabs: int = 20,
           ewald_grid: tuple[int, int, int] = (10, 10, 10), stride: int = 1,
           temperature: float = 300.0, namd: str = 'namd3', nprocs: int = 1,
           prefix: str = 'ppreplay') -> ReplayResult:
    """Run both passes over ``dcd`` and return the reconstructed profile."""
    base = read_config_lines(namd_config)
    kept, dropped = strip_for_replay(base)
    if dropped:
        logger.info(f'dropped {len(dropped)} production directive(s) from {namd_config}: '
                    + ', '.join(sorted({d.strip().split()[0] for d in dropped})))
    os.makedirs(workdir, exist_ok=True)
    # NAMD runs in the production config's own directory so that every relative path it names --
    # PSF, parameters, restraint and colvars files -- resolves the way it did in production.
    # Paths this module supplies are therefore written relative to that directory, not to workdir.
    basedir = os.path.dirname(os.path.abspath(namd_config)) or '.'
    dcd_abs = os.path.abspath(dcd)
    passes = {}
    for mode in ('real', 'ewald'):
        cfg = os.path.join(workdir, f'{prefix}-{mode}.namd')
        log = os.path.join(workdir, f'{prefix}-{mode}.log')
        outname = os.path.abspath(os.path.join(workdir, f'{prefix}-{mode}-out'))
        with open(cfg, 'w') as f:
            f.write(build_replay_config(base, dcd_abs, mode, slabs=slabs, outname=outname,
                                        basedir=basedir, ewald_grid=ewald_grid, stride=stride,
                                        temperature=temperature))
        run_namd(cfg, log, basedir, namd=namd, nprocs=nprocs)
        p = parse_profiles(log, mode)
        p.config_path = cfg
        logger.info(f'{mode} pass: {p.nframes} frames x {p.nslabs} slabs')
        passes[mode] = p
    total = combine(passes['real'], passes['ewald'])
    return ReplayResult(real=passes['real'], ewald=passes['ewald'], total=total)


def write_csv(result: ReplayResult, path: str) -> None:
    """Write the frame-averaged profile: z, then the three components of each contribution."""
    z = result.slab_centers()
    real = result.real.profiles.mean(axis=0)
    ewald = result.ewald.profiles.mean(axis=0)
    tot = result.mean_profile()
    sem = result.sem_profile()
    with open(path, 'w') as f:
        f.write('z,real_xx,real_yy,real_zz,ewald_xx,ewald_yy,ewald_zz,'
                'total_xx,total_yy,total_zz,total_sem_xx,total_sem_yy,total_sem_zz\n')
        for i in range(result.nslabs):
            vals = [z[i], *real[i], *ewald[i], *tot[i], *sem[i]]
            f.write(','.join(f'{v:.6g}' for v in vals) + '\n')


def plot(result: ReplayResult, path: str, *, title: str = '', figsize=(15, 5),
         stamp: str = '') -> None:
    """Three panels -- lateral, normal, isotropic -- each showing both halves and their sum.

    Plotting the two contributions alongside the total is the point rather than a diagnostic
    flourish: it is what shows at a glance that the reciprocal term is not a small correction.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    z = result.slab_centers()
    zlabel = 'z (Å)' if result.real.slab_thickness is not None else 'slab index'

    def components(arr):
        """(lateral, normal, isotropic) for a (nslabs, 3) mean profile."""
        xx, yy, zz = arr[:, 0], arr[:, 1], arr[:, 2]
        return (xx + yy) / 2.0, zz, (xx + yy + zz) / 3.0

    real = components(result.real.profiles.mean(axis=0))
    ewald = components(result.ewald.profiles.mean(axis=0))
    tot = components(result.mean_profile())
    sem = components(result.sem_profile())

    labels = [r'$\frac{1}{2}(P_{xx}+P_{yy})$', r'$P_{zz}$',
              r'$\frac{1}{3}(P_{xx}+P_{yy}+P_{zz})$']
    fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=True)
    for k, ax in enumerate(axes):
        ax.plot(real[k], z, color='#3366aa', lw=1.0, ls='--', label='real space only')
        ax.plot(ewald[k], z, color='#cc8800', lw=1.0, ls=':', label='reciprocal (Ewald) only')
        ax.plot(tot[k], z, color='#cc3333', lw=1.6, label='complete (sum)')
        ax.fill_betweenx(z, tot[k] - sem[k], tot[k] + sem[k], color='#cc3333', alpha=0.2)
        ax.axvline(0, color='k', ls='--', lw=0.8)
        ax.set_xlabel(labels[k] + ' (bar)')
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(zlabel)
    axes[2].legend(fontsize=8, loc='best', framealpha=0.9)
    fig.suptitle(title or f'complete pressure profile ({result.nframes} frames)', fontsize=10)
    fig.tight_layout()
    if stamp:
        fig.text(0.995, 0.005, stamp, ha='right', va='bottom', fontsize=6, color='0.4')
    fig.savefig(path, dpi=120, bbox_inches='tight')
    plt.close(fig)
