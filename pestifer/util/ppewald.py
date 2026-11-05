# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Offline reconstruction of a *complete* NAMD pressure profile, reciprocal-space term included.

Why this module exists
----------------------

A NAMD run with ``pressureProfile on`` reports "kinetic, bonded and nonbonded (but not reciprocal
space) contributions".  The missing PME reciprocal term is not a small correction -- on a DMPC
patch it is about -180 bar of the trace and +33 mN/m of the integrated tension -- and it cannot be
recovered by also switching on ``pressureProfileEwald``, because NAMD makes the two mutually
exclusive (``SimParameters.C``, in 3.0.3 at line 6699)::

    // if Ewald is on, only calculate Ewald
    if (pressureProfileEwaldOn)
      pressureProfileOn = 0;

Running the Ewald profile inline is impractical: ``ComputeEwald`` is registered as an ordinary
compute gated only on ``doFullElectrostatics``, never on ``pressureProfileFreq``, so it evaluates
every full-electrostatics step no matter how rarely profiles are sampled -- measured at 14-28x the
per-step cost.  Replaying the saved trajectory evaluates it once per *sampled frame* instead.

What is summed
--------------

The complete profile is the **production run's own** ``PRESSUREPROFILE`` records plus the
reciprocal profile replayed from its DCD, matched frame by frame **by MD step**.  Two things that
earlier versions did are wrong, and both were measured (2026-11-05, 200-lipid DMPC patch):

- *A real-space profile replayed at* ``run 0`` *is not the run's real-space profile.*  The replay
  draws one velocity set from ``temperature`` and reuses it for every frame, and a single force
  evaluation carries no rigid-bond constraint virial.  It disagreed with the production profile by
  tens to hundreds of bar per component (anisotropically), and a self-check that compared the
  replay with the replay's own ``PRESSURE`` could not see it, since both shared those velocities.
  So no real-space pass is replayed any more.
- *The production run must not use multiple time stepping.*  At ``fullElectFrequency 2`` NAMD's
  ordinary profile already carries part of the reciprocal term, so adding the replay counts it
  twice: 70-95 bar per component, +14.5 mN/m of total tension that is not there.  At
  ``fullElectFrequency 1`` the sum matched NAMD's own pressure tensor to about 1 bar per component.
  :func:`replay` refuses an MTS production config unless told to accept it.

Self-check
----------

When the production run wrote its pressure tensor (``outputPressure`` equal to
``pressureProfileFreq``, so that ``PRESSAVG`` averages the same window as each profile record), the
slab average of the complete profile is compared with ``PRESSAVG`` **per component**.  The trace
alone cannot see an anisotropic error, and an anisotropic error is exactly what a surface tension
is made of.
"""
from __future__ import annotations

import logging
import math
import os
import re
import shutil
import subprocess

from dataclasses import dataclass, field

import numpy as np

from .dcd import read_header, iter_frames

logger = logging.getLogger(__name__)


_PROFILE_RE = re.compile(r'^PRESSUREPROFILE:\s+(\S+)\s+(.*)$')
_PRESSAVG_RE = re.compile(r'^PRESSAVG:\s+(\S+)\s+(.*)$')
_SLAB_THICKNESS_RE = re.compile(r'SLAB THICKNESS:\s*([-\d.eE+]+)')
_FRAMES_RE = re.compile(r'PPREPLAY_FRAMES\s+(\d+)')

#: Target reciprocal-grid spacing, A, for an automatically sized ``pressureProfileEwald`` grid.
#: Measured on example 17's ~115 x 116 x 81 A box (pestifer-sweep, 2026-11-05): the reconstruction
#: missed NAMD's pressure by 182 bar at 10^3 (NAMD's default, ~11 A), 53 bar at 20^3 (~5.8 A),
#: 1.2 bar at 30^3 (~3.8 A) and 0.9 bar at 40^3 (~2.9 A).  3 A is inside the converged range.
EWALD_GRID_SPACING = 3.0

#: Largest per-component disagreement, in bar, between the complete profile's slab average and the
#: production run's PRESSAVG tensor that still counts as closure (or three standard errors of the
#: per-frame difference, if that is larger).  A sound reconstruction agrees to about 1 bar; the
#: double-counting at fullElectFrequency 2 misses by 70-95.
CLOSURE_TOLERANCE_BAR = 10.0


@dataclass
class ProfilePass:
    """Per-frame profiles from one source: the production log (``'real'``) or the replay
    (``'ewald'``)."""

    mode: str
    profiles: np.ndarray
    """``(nframes, nslabs, 3)`` array of per-slab ``(Pxx, Pyy, Pzz)`` in bar."""
    slab_thickness: float | None = None
    """Angstroms, as NAMD reported it at startup; ``None`` if the log did not say."""
    frames_read: int | None = None
    """Frames the replay loop actually read, from the loop's own counter."""
    log_path: str = ''
    config_path: str = ''

    @property
    def nframes(self) -> int:
        return int(self.profiles.shape[0])

    @property
    def nslabs(self) -> int:
        return int(self.profiles.shape[1])


@dataclass
class ReplayResult:
    """The complete profile, the two parts it came from, and what it can be checked against."""

    real: ProfilePass
    ewald: ProfilePass
    total: np.ndarray = field(repr=False)
    """``(nframes, nslabs, 3)``: the elementwise sum, which is the complete profile."""
    steps: np.ndarray | None = None
    """MD step of each frame."""
    tensor: np.ndarray | None = None
    """``(nframes, 3)``: the diagonal of the production run's PRESSAVG at each frame's step, or
    ``None`` when the production log carries no pressure tensor."""
    ewald_grid: tuple[int, int, int] | None = None

    @property
    def nframes(self) -> int:
        return int(self.total.shape[0])

    @property
    def nslabs(self) -> int:
        return int(self.total.shape[1])

    def slab_centers(self) -> np.ndarray:
        """Slab-center coordinates in A, measured from the bottom of the cell (``origin.z - c/2``),
        or slab indices if NAMD gave no thickness."""
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


def config_value(lines: list[str], key: str) -> str | None:
    """The value of the last ``key`` directive in a NAMD config, matched case-insensitively."""
    value = None
    for raw in lines:
        s = raw.split('#', 1)[0].strip()
        if not s:
            continue
        parts = s.split(None, 1)
        if parts[0].lower() == key.lower() and len(parts) > 1:
            value = parts[1].strip()
    return value


def full_elect_frequency(lines: list[str]) -> int:
    """``fullElectFrequency`` as the config sets it (NAMD's default is 1)."""
    v = config_value(lines, 'fullElectFrequency')
    try:
        return int(float(v)) if v is not None else 1
    except ValueError:
        return 1


def auto_ewald_grid(dcd: str, spacing: float = EWALD_GRID_SPACING) -> tuple[int, int, int]:
    """A ``pressureProfileEwald`` grid sized from the first frame's cell: about ``spacing`` A."""
    _, cell, _ = next(iter_frames(dcd, atoms=np.array([0]), frames=[0]))
    if cell is None:
        raise ValueError(f'{dcd}: no unit cell in the DCD (write it with DCDUnitCell yes)')
    return tuple(max(10, int(math.ceil(L / spacing))) for L in cell)


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
    # the replay sets its own working directory -- from this one, if the production config had
    # it -- and NAMD rejects a config that names two
    'cwd',
})

#: Tcl commands that would run dynamics or otherwise drive the engine.  A production config often
#: ends in a bare ``run 500000``; a scripted one may loop.  Lines starting with these are dropped
#: along with :data:`DROP_DIRECTIVES` so that only the replay loop drives NAMD.
_TCL_DRIVERS = frozenset({'run', 'minimize', 'move', 'moveallby', 'coorfile', 'reinitatoms'})


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



def build_replay_config(base_lines: list[str], dcd: str, *, slabs: int, outname: str,
                        ewald_grid: tuple[int, int, int], basedir: str = '', stride: int = 1,
                        temperature: float = 300.0) -> str:
    """Render the text of the reciprocal-space (``pressureProfileEwald``) replay config."""
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
        '# The reciprocal term depends on coordinates and cell only, so the velocities this',
        '# temperature draws do not enter it.',
        f'temperature             {temperature}',
        # NAMD requires outputName even when nothing is integrated.  The production one was
        # dropped so that a replay can never overwrite the run it is reading.
        f'outputName              {outname}',
        'firsttimestep           0',
        'numsteps                0',
        'outputEnergies          1',
        'CUDASOAintegrate        off',
        '',
        # pressureProfile must be on as well, or NAMD writes no PRESSUREPROFILE records at all
        # (measured); with pressureProfileEwald on it then reports the reciprocal term ALONE.
        'pressureProfile         on',
        f'pressureProfileSlabs    {slabs}',
        f'pressureProfileFreq     1',
        '# NAMD sets pressureProfileOn = 0 when this is on, so this pass reports the',
        '# reciprocal-space contribution ALONE.  It is summed with the production profile.',
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


def _profile_array(rows: list[list[float]], log_path: str) -> np.ndarray:
    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise RuntimeError(f'{log_path}: ragged PRESSUREPROFILE records (widths {sorted(widths)})')
    w = widths.pop()
    if w % 3:
        raise RuntimeError(f'{log_path}: {w} values per record is not a multiple of 3 (xx,yy,zz)')
    return np.asarray(rows, dtype=float).reshape(len(rows), w // 3, 3)


def parse_profiles(log_path: str, mode: str = 'ewald') -> ProfilePass:
    """Read every ``PRESSUREPROFILE:`` record out of a replay log, in order.

    Records are kept in file order deliberately: at ``run 0`` every frame is evaluated at timestep
    0, so the step column is a constant and carries no frame identity.
    """
    rows: list[list[float]] = []
    slab_thickness = None
    frames_read = None
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
    if not rows:
        raise RuntimeError(
            f'{log_path}: no PRESSUREPROFILE records. The replay produced no profile -- check the '
            'log for a NAMD error, and that the trajectory atom count matches the PSF.')
    return ProfilePass(mode=mode, profiles=_profile_array(rows, log_path),
                       slab_thickness=slab_thickness, frames_read=frames_read, log_path=log_path)


def parse_production_log(log_path: str):
    """``(profiles, pressavg, slab_thickness)`` from a production NAMD log, keyed by MD step.

    ``profiles[step]`` is ``(nslabs, 3)``; ``pressavg[step]`` is the diagonal ``(xx, yy, zz)`` of
    the ``PRESSAVG`` tensor (present only with ``outputPressure`` set).  A step seen twice -- a log
    holding a restart -- keeps its last record.
    """
    profiles, pressavg = {}, {}
    slab_thickness = None
    with open(log_path, errors='ignore') as f:
        for line in f:
            m = _PROFILE_RE.match(line)
            if m:
                profiles[int(float(m.group(1)))] = [float(x) for x in m.group(2).split()]
                continue
            m = _PRESSAVG_RE.match(line)
            if m:
                t = [float(x) for x in m.group(2).split()]
                if len(t) == 9:
                    pressavg[int(float(m.group(1)))] = np.array([t[0], t[4], t[8]])
                continue
            if slab_thickness is None:
                m = _SLAB_THICKNESS_RE.search(line)
                if m:
                    slab_thickness = float(m.group(1))
    if not profiles:
        raise RuntimeError(
            f'{log_path}: no PRESSUREPROFILE records.  The production run must have had '
            '`pressureProfile on`: its profile is the real-space half of the reconstruction.')
    steps = sorted(profiles)
    arr = _profile_array([profiles[s] for s in steps], log_path)
    return ({s: arr[i] for i, s in enumerate(steps)}, pressavg, slab_thickness)


def align_production(log_path: str, steps: np.ndarray):
    """The production profile and PRESSAVG diagonal at each of ``steps``.

    Returns ``(ProfilePass, tensor)``; ``tensor`` is ``None`` unless every step has a PRESSAVG.
    Raises if any step has no production profile record, naming the first few: a DCD written more
    often than the profile, or a log from a different run, cannot be paired frame by frame.
    """
    profiles, pressavg, th = parse_production_log(log_path)
    missing = [int(s) for s in steps if int(s) not in profiles]
    if missing:
        have = sorted(profiles)
        raise RuntimeError(
            f'{log_path}: {len(missing)} of {len(steps)} trajectory frames have no production '
            f'PRESSUREPROFILE record at their step (first: {missing[:5]}; the log has steps '
            f'{have[0]}..{have[-1]}).  Write the DCD at a multiple of pressureProfileFreq, and '
            'pass the log of the run that wrote this DCD.')
    real = ProfilePass(mode='real', profiles=np.array([profiles[int(s)] for s in steps]),
                       slab_thickness=th, log_path=log_path)
    tensor = None
    if pressavg and all(int(s) in pressavg for s in steps):
        tensor = np.array([pressavg[int(s)] for s in steps])
    return real, tensor


def combine(real: ProfilePass, ewald: ProfilePass) -> np.ndarray:
    """Elementwise sum of the two parts, which is the complete profile.

    Refuses mismatched parts: the replay's records carry no step, so its frame count is the only
    evidence that it read the frames the production records were matched to.
    """
    if real.nframes != ewald.nframes:
        raise RuntimeError(
            f'frame counts differ (production {real.nframes}, ewald replay {ewald.nframes}); the '
            'replay did not read the frames the production records were matched to')
    if real.nslabs != ewald.nslabs:
        raise RuntimeError(
            f'slab counts differ (production {real.nslabs}, ewald replay {ewald.nslabs}); summing '
            'different slab griddings is meaningless')
    return real.profiles + ewald.profiles


def check_against_tensor(result: ReplayResult) -> dict | None:
    """Slab-averaged complete profile minus the production PRESSAVG diagonal, per component.

    Returns ``None`` when the production log had no pressure tensor.  ``deviation`` and
    ``real_only_deviation`` are frame means, ``(xx, yy, zz)`` in bar; ``sem`` is the standard error
    of the per-frame difference (frames treated as independent, so a floor).
    """
    if result.tensor is None:
        return None
    d = result.total.mean(axis=1) - result.tensor
    n = len(d)
    sem = d.std(axis=0, ddof=1) / np.sqrt(n) if n > 1 else np.zeros(3)
    return {
        'deviation': d.mean(axis=0),
        'sem': sem,
        'real_only_deviation': (result.real.profiles.mean(axis=1) - result.tensor).mean(axis=0),
        'nframes': n,
    }


def closure_warnings(check: dict) -> list[str]:
    """Reasons to distrust a reconstruction, from :func:`check_against_tensor`'s result."""
    dev, sem = np.asarray(check['deviation']), np.asarray(check['sem'])
    limit = np.maximum(CLOSURE_TOLERANCE_BAR, 3.0 * sem)
    if np.all(np.abs(dev) <= limit):
        return []
    return [
        'the complete profile does not reproduce the production run\'s own pressure tensor: its '
        'slab average differs from PRESSAVG by '
        + ' / '.join(f'{v:+.1f}' for v in dev)
        + f' bar (xx / yy / zz; a sound reconstruction agrees to about 1 bar, and more than '
        f'{CLOSURE_TOLERANCE_BAR:.0f} is not closure).  Likely causes, most common first: '
        '(1) the Ewald grid is too coarse for this box -- re-run with a finer --ewald-grid and see '
        'whether the gap shrinks; (2) the production config and this replay do not compute the '
        'same forces (the config was edited after the run, or the log is from a different run); '
        '(3) a NAMD build that writes wrong pressure profiles while its energies are right -- '
        'that one misses by thousands of bar, not tens.  Do not use this profile until it closes.']


def replay(namd_config: str, dcd: str, production_log: str, *, workdir: str = '.',
           ewald_grid: tuple[int, int, int] | None = None, stride: int = 1,
           temperature: float = 300.0, namd: str = 'namd3', nprocs: int = 1,
           prefix: str = 'ppreplay', accept_mts: bool = False) -> ReplayResult:
    """Replay the reciprocal term over ``dcd`` and add it to the production profile.

    Everything that can fail without NAMD -- the MTS check, the DCD header, matching every frame
    to a production record -- is done before NAMD starts.
    """
    base = read_config_lines(namd_config)
    fef = full_elect_frequency(base)
    if fef > 1:
        msg = (f'{namd_config} sets fullElectFrequency {fef}.  Under multiple time stepping the '
               'production profile already carries part of the reciprocal term, so adding the '
               'replay counts it twice (measured: 70-95 bar per component, +14.5 mN/m of total '
               'tension on a DMPC patch), and no rescaling repairs it.  Re-run production with '
               'fullElectFrequency 1.')
        if not accept_mts:
            raise RuntimeError(msg + '  (--accept-mts proceeds anyway.)')
        logger.warning(msg + '  Proceeding because MTS was explicitly accepted.')
    hdr = read_header(dcd)
    if not hdr.has_cell:
        raise ValueError(f'{dcd}: no unit cell in the DCD (write it with DCDUnitCell yes)')
    steps = hdr.steps()[::max(1, stride)]
    real, tensor = align_production(production_log, steps)
    _, dropped = strip_for_replay(base)
    if dropped:
        logger.info(f'dropped {len(dropped)} production directive(s) from {namd_config}: '
                    + ', '.join(sorted({d.strip().split()[0] for d in dropped})))
    if ewald_grid is None:
        ewald_grid = auto_ewald_grid(dcd)
        logger.info(f'Ewald grid {ewald_grid[0]} x {ewald_grid[1]} x {ewald_grid[2]}, sized from '
                    f'the first frame\'s cell at ~{EWALD_GRID_SPACING:g} A')
    os.makedirs(workdir, exist_ok=True)
    # NAMD runs in the production config's own directory so that every relative path it names --
    # PSF, parameters, restraint and colvars files -- resolves the way it did in production.
    # A production config that sets its own `cwd` resolves its paths from there instead.
    basedir = os.path.dirname(os.path.abspath(namd_config)) or '.'
    own_cwd = config_value(base, 'cwd')
    if own_cwd:
        basedir = os.path.normpath(os.path.join(basedir, own_cwd))
    cfg = os.path.join(workdir, f'{prefix}-ewald.namd')
    log = os.path.join(workdir, f'{prefix}-ewald.log')
    outname = os.path.abspath(os.path.join(workdir, f'{prefix}-ewald-out'))
    with open(cfg, 'w') as f:
        f.write(build_replay_config(base, os.path.abspath(dcd), slabs=real.nslabs,
                                    outname=outname, basedir=basedir,
                                    ewald_grid=tuple(ewald_grid), stride=stride,
                                    temperature=temperature))
    run_namd(cfg, log, basedir, namd=namd, nprocs=nprocs)
    ewald = parse_profiles(log, 'ewald')
    ewald.config_path = cfg
    logger.info(f'ewald replay: {ewald.nframes} frames x {ewald.nslabs} slabs')
    total = combine(real, ewald)
    return ReplayResult(real=real, ewald=ewald, total=total, steps=np.asarray(steps),
                        tensor=tensor, ewald_grid=tuple(ewald_grid))


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
        ax.plot(real[k], z, color='#3366aa', lw=1.0, ls='--', label='real space (production run)')
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
