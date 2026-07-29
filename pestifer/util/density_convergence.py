# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Density-convergence machinery for the :class:`DensityEquilibrateTask
<pestifer.tasks.density_equilibrate.DensityEquilibrateTask>`.

This module is deliberately free of any NAMD / VMD / task machinery so it can be unit-tested in
isolation.  It provides three separable pieces, mirroring ``docs/design/density-equilibrate.md``:

1. **Density from the ``.xst``.**  NAMD writes the extended-system trajectory one line per frame:
   ``ts  a_x a_y a_z  b_x b_y b_z  c_x c_y c_z  ...``.  The cell volume is the scalar triple product
   of the three cell vectors; combined with the total system mass from the PSF this gives a scalar
   mass-density time series (g/cc).

2. **The stability-bounded chunk length.**  Chunking is *mandatory* for numerical stability: NAMD
   fixes its patch/PME grid at the start of each ``run``, so a chunk must be short enough that the
   shrinking cell does not outrun the patch ``margin`` within it.  :func:`next_chunk_steps` sizes the
   next chunk from the shrink rate just observed, so chunks are short early (fast shrink) and longer
   late (box settled) -- the ladder's shape, derived rather than guessed.  Because NAMD's *effective*
   shrink budget is hard to model exactly (it depends on the patch grid it chose and, on CPU, on a
   ``margin`` that may differ from ours), the predictive sizing is a conservative estimate backed by a
   **reactive safety net**: :func:`is_patch_grid_crash` recognizes the specific NAMD abort so the task
   can roll back and retry shorter, and :func:`parse_patch_grid` reads the grid NAMD actually built
   for observability.

3. **The convergence criterion.**  A *practical fractional-drift* tolerance plus a **precision gate**
   that grows the window until the mean is genuinely resolved above noise, with ``n_consecutive``
   hysteresis.  The gate uses an **autocorrelation-corrected SEM**: NPT cell density is autocorrelated
   over hundreds of steps (an overdamped Langevin piston -- measured integrated autocorrelation time
   ~560 steps for a small TIP3 box), so consecutive ``.xst`` frames are not independent and the honest
   standard error of the mean is ``sigma/sqrt(N_eff)`` with ``N_eff = N/tau_int``.  This also makes the
   gate size-aware for free (``sigma/mean ~ 1/sqrt(N_atoms)``, ``tau`` ~size-independent).  See the
   design doc for the rejected alternatives (fluctuation threshold, ``trend/SEM`` significance) and for
   the fine-sampling (``xstfreq=1``) experiment that measured the autocorrelation.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

import numpy as np

from .densityprofile import AMU_PER_A3_TO_G_PER_CC

logger = logging.getLogger(__name__)

#: Absolute minimum windowed samples before an autocorrelation time (and hence the SEM) is trusted --
#: guards against tau_int being under-estimated from a too-short series (see ``check``).
_MIN_WINDOW_SAMPLES = 24


def _xst_data_rows(path):
    """Yield the numeric fields of each non-comment data line of a NAMD ``.xst`` file."""
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            yield [float(x) for x in s.split()]


def xst_cell_volumes(path):
    """Return ``(timesteps, volumes)`` from a NAMD ``.xst`` file.

    Volume is the scalar triple product ``|a . (b x c)|`` of the three cell vectors, so it is correct
    for any (orthorhombic or triclinic) cell.  Returns two 1-D numpy arrays; empty arrays if the file
    has no data lines yet."""
    ts, vol = [], []
    for v in _xst_data_rows(path):
        a = np.array(v[1:4])
        b = np.array(v[4:7])
        c = np.array(v[7:10])
        ts.append(v[0])
        vol.append(abs(float(np.dot(a, np.cross(b, c)))))
    return np.array(ts, dtype=float), np.array(vol, dtype=float)


def xst_cell_areas(path):
    """Return ``(timesteps, lateral_areas)`` from a NAMD ``.xst`` file.

    The membrane normal is taken along **z**, so the lateral area is ``|a x b|`` of the two in-plane
    cell vectors (``= a_x * b_y`` for an orthorhombic cell).  This is the physically meaningful membrane
    observable (area-per-lipid = area / lipids-per-leaflet).  Returns two 1-D numpy arrays; empty arrays
    if the file has no data lines yet.  Companion to :func:`xst_cell_volumes` for the membrane-aware
    (density + area) equilibration."""
    ts, area = [], []
    for v in _xst_data_rows(path):
        a = np.array(v[1:4])
        b = np.array(v[4:7])
        ts.append(v[0])
        area.append(float(np.linalg.norm(np.cross(a, b))))
    return np.array(ts, dtype=float), np.array(area, dtype=float)


def total_mass_amu(psf_path):
    """Sum the per-atom masses (amu) from the ``!NATOM`` block of an XPLOR/CHARMM PSF."""
    with open(psf_path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines) and '!NATOM' not in lines[i]:
        i += 1
    if i == len(lines):
        raise ValueError(f'{psf_path}: no !NATOM record found; not a PSF?')
    natom = int(lines[i].split()[0])
    m = 0.0
    for line in lines[i + 1:i + 1 + natom]:
        m += float(line.split()[7])
    return m


def integrated_autocorr_time(x):
    """Integrated autocorrelation time ``tau_int`` (in samples) of a 1-D series, via the standard
    ``tau = 1 + 2 * sum_k rho_k`` with the sum truncated at the first non-positive lag (a robust,
    slightly-conservative automatic window).  The effective number of independent samples is then
    ``N / tau_int`` and the standard error of the mean is ``sigma / sqrt(N / tau_int)``.

    Returns ``1.0`` (i.e. samples treated as independent) for a series too short or with zero variance.
    Cost is O(N * tau) in practice because the sum stops at the first zero crossing."""
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 4:
        return 1.0
    x = x - x.mean()
    var = float(np.dot(x, x) / n)
    if var <= 0.0:
        return 1.0
    tau = 1.0
    for k in range(1, n):
        rho = float(np.dot(x[:n - k], x[k:]) / ((n - k) * var))
        if rho <= 0.0:
            break
        tau += 2.0 * rho
    return tau


def total_atoms(psf_path):
    """Return the atom count from the ``!NATOM`` record of an XPLOR/CHARMM PSF (for the size gate)."""
    with open(psf_path) as f:
        for line in f:
            if '!NATOM' in line:
                return int(line.split()[0])
    raise ValueError(f'{psf_path}: no !NATOM record found; not a PSF?')


def volume_to_density(volume_a3, mass_amu):
    """Convert a cell volume (angstrom^3) + total mass (amu) to mass density (g/cc)."""
    return AMU_PER_A3_TO_G_PER_CC * mass_amu / np.asarray(volume_a3, dtype=float)


def xst_max_shrink_rate(path):
    """Largest per-dimension cell-edge *shrink* rate (angstrom/step) across a ``.xst`` chunk.

    Uses the three cell-vector lengths as the "dimensions"; returns ``max(0, (len_first - len_last)))``
    per dimension divided by the elapsed steps.  Growth (a cell that expands) contributes 0.  Returns
    ``0.0`` if the chunk has fewer than two frames or spans zero steps -- callers then fall back to
    the maximum chunk length."""
    rows = list(_xst_data_rows(path))
    if len(rows) < 2:
        return 0.0
    first, last = rows[0], rows[-1]
    dsteps = last[0] - first[0]
    if dsteps <= 0:
        return 0.0

    def lengths(v):
        return (np.linalg.norm(v[1:4]), np.linalg.norm(v[4:7]), np.linalg.norm(v[7:10]))

    l0, l1 = lengths(first), lengths(last)
    shrink = max((a - b) for a, b in zip(l0, l1))  # most-shrunk dimension; negative if all grew
    return max(0.0, shrink) / dsteps


def next_chunk_steps(shrink_rate, margin, shrink_safety, chunk_min, chunk_max):
    """Stability-bounded length (steps) for the next NPT chunk.

    The next chunk is sized so the projected worst-case per-dimension cell shrink stays below
    ``shrink_safety * margin`` (angstrom), given the ``shrink_rate`` (angstrom/step) just observed::

        chunk = clamp(shrink_safety * margin / shrink_rate, chunk_min, chunk_max)

    A zero/near-zero shrink rate (settled box) yields ``chunk_max``; a fast-shrinking box is throttled
    toward ``chunk_min``.  Result is an int in ``[chunk_min, chunk_max]``."""
    if shrink_rate <= 0:
        return int(chunk_max)
    allowed = shrink_safety * margin / shrink_rate
    return int(max(chunk_min, min(chunk_max, allowed)))


def quantize_steps(n, steps_per_cycle, minimum=None):
    """Round ``n`` down to a positive multiple of ``steps_per_cycle``.

    NAMD aborts (``FATAL ERROR: number of steps must be a multiple of stepsPerCycle``) if a ``run``
    count is not a whole number of cycles, so every chunk length must be quantized before use.  Rounds
    *down* (never overshoot a shrink budget) but floors at one cycle -- or, if ``minimum`` is given, at
    ``minimum`` rounded *up* to a whole number of cycles.  Returns an int."""
    spc = max(1, int(steps_per_cycle))
    q = (int(n) // spc) * spc
    if minimum is None:
        floor = spc
    else:
        floor = ((int(minimum) + spc - 1) // spc) * spc
    return max(floor, q)


def is_patch_grid_crash(log_text):
    """True if a NAMD log records the "cell too small for the patch grid" abort.

    This is the specific failure the chunking exists to avoid: the barostat shrank the cell past what
    the patch grid (fixed at ``run`` start) can represent.  NAMD's message is *"Periodic cell has
    become too small for original patch grid! ..."*; we match on the stable ``too small`` +
    ``patch grid`` fragments so a wording change does not silently defeat the reactive retry.  Any
    *other* NAMD failure returns False (the caller re-raises it -- it is not something a shorter chunk
    fixes)."""
    low = (log_text or '').lower()
    return 'too small' in low and 'patch grid' in low


# NAMD startup log fragments (all "Info:" lines).  Kept lenient: PBC on/off wording varies by build.
_PATCH_GRID_RE = re.compile(r'PATCH GRID IS\s+(\d+)[^\d]*?BY\s+(\d+)[^\d]*?BY\s+(\d+)', re.I)
_CELL_BASIS_RE = re.compile(
    r'PERIODIC CELL BASIS\s+([123])\s+([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)', re.I)
_CUTOFF_RE = re.compile(r'\bCUTOFF\s+([\d.eE+-]+)', re.I)
_PAIRLIST_RE = re.compile(r'PAIRLIST DISTANCE\s+([\d.eE+-]+)', re.I)
_MARGIN_RE = re.compile(r'\bMARGIN\s+([\d.eE+-]+)', re.I)


def parse_patch_grid(log_text):
    """Read the patch grid NAMD actually built, for observability (not correctness).

    Returns a dict with whatever could be parsed from the startup log -- ``patches`` (nx, ny, nz),
    ``cell_edges`` (lengths of the three periodic cell-basis vectors), ``cutoff``, ``pairlistdist``,
    ``margin``, and, when both the grid and the cell are known, the minimum current patch dimension
    ``min_patch_dim`` (= min edge / its patch count) and -- *only when it is positive* -- an
    *approximate* ``shrink_headroom`` before that dimension reaches the pairlist distance.  Returns
    ``None`` if no patch grid line is present.

    The headroom is explicitly approximate -- NAMD's exact abort threshold is internal -- and is for
    logging/tuning only; the reactive crash-and-retry, not this number, is what guarantees stability."""
    if not log_text:
        return None
    m = _PATCH_GRID_RE.search(log_text)
    if not m:
        return None
    out = {'patches': (int(m.group(1)), int(m.group(2)), int(m.group(3)))}

    edges = {}
    for bm in _CELL_BASIS_RE.finditer(log_text):
        idx = int(bm.group(1))
        vec = np.array([float(bm.group(2)), float(bm.group(3)), float(bm.group(4))])
        edges[idx] = float(np.linalg.norm(vec))
    if len(edges) == 3:
        out['cell_edges'] = (edges[1], edges[2], edges[3])

    for key, rx in (('cutoff', _CUTOFF_RE), ('pairlistdist', _PAIRLIST_RE), ('margin', _MARGIN_RE)):
        mm = rx.search(log_text)
        if mm:
            out[key] = float(mm.group(1))

    if 'cell_edges' in out:
        patch_dims = [e / n for e, n in zip(out['cell_edges'], out['patches'])]
        out['min_patch_dim'] = min(patch_dims)
        floor = out.get('pairlistdist', out.get('cutoff'))
        # Report headroom only when it is meaningfully positive.  The "shrink until a patch reaches
        # the pairlist distance" heuristic assumes the cutoff-limited multi-patch decomposition NAMD
        # builds on CPU/multi-node.  GPU-resident NAMD tiles the whole system into *few* patches
        # (its single-device scheme), so a small box can already have min_patch_dim < pairlistdist
        # while running perfectly fine -- a negative "headroom" there is nonsense, not a warning.  In
        # that case we omit the number rather than print a misleading one (stability is guaranteed by
        # the reactive crash-and-retry regardless).
        if floor is not None and out['min_patch_dim'] - floor > 0:
            out['shrink_headroom'] = out['min_patch_dim'] - floor
    return out


@dataclass
class ConvergenceParams:
    """Tunables for the density-convergence criterion (see the design doc's Parameters table)."""
    drift_tol: float = 2e-3      #: converged when the *upper confidence bound* on |drift| < this (~0.2%)
    precision_p: float = 3.0     #: precision gate: require SEM/mean < drift_tol / precision_p
    drift_conf: float = 0.0      #: sigmas added to |drift| for its upper bound; 0 = point estimate
                                 #: (the honest precision gate already guards against premature passes,
                                 #: so the bound is redundant; raise it only for extra trend rigor)
    burn_in: int = 2000          #: leading steps discarded before assessing the trend
    window_frac: float = 0.5     #: assess only the trailing this-fraction of post-burn-in samples
    min_steps: int = 4000        #: never declare convergence before this many steps
    n_consecutive: int = 3       #: successive passing checks required before stopping (hysteresis)
    # Honest, autocorrelation-corrected SEM.  The NPT cell density is autocorrelated over ~hundreds of
    # steps (measured integrated autocorr time ~560 steps for a small TIP3 box -- an overdamped Langevin
    # piston, not fast noise), so consecutive .xst frames are NOT independent.  The SEM of the mean is
    # therefore sigma/sqrt(N_eff) with N_eff = N_window / tau_int, not sigma/sqrt(N) -- which is why a
    # naive block-means SEM (blocks treated as independent) is optimistic by ~1.5x.  This also makes the
    # gate size-aware *for free*: sigma/mean ~ 1/sqrt(N_atoms) while tau is ~size-independent, so a large
    # box has a smaller honest SEM at equal duration and clears a fixed gate sooner -- no ad-hoc size
    # factor needed.  We refuse to trust the SEM (and so cannot converge) until the window is at least
    # `autocorr_reliability` correlation times long, since tau_int is not estimable from a short series.
    autocorr_reliability: float = 6.0  #: require window length >= this * tau_int before trusting the SEM


@dataclass
class ConvergenceReport:
    """Verdict + diagnostics from one :meth:`DensityConvergenceMonitor.check` call."""
    converged: bool = False
    passes: int = 0              #: consecutive passing checks accumulated so far
    drift: float | None = None   #: fractional drift over the window (|slope|*span/mean)
    drift_hi: float | None = None  #: upper confidence bound on |drift| (|drift| + drift_conf*SE)
    sem_over_mean: float | None = None
    precision_met: bool = False
    mean_density: float | None = None
    n_window: int = 0            #: samples in the post-burn-in window
    tau_int: float | None = None #: integrated autocorrelation time (in frames) over the window
    n_eff: float | None = None   #: effective independent samples = n_window / tau_int
    last_step: float | None = None
    blowup: bool = False         #: a non-finite density was seen -> caller must abort
    reason: str = ''             #: human-readable one-liner for the report/log


class DensityConvergenceMonitor:
    """Accumulate ``(timestep, density)`` samples and decide, at each chunk boundary, whether the
    box density has stopped drifting.

    Usage::

        mon = DensityConvergenceMonitor(ConvergenceParams())
        for chunk:
            mon.add_samples(times, densities)
            rep = mon.check()
            if rep.blowup: abort
            if rep.converged: break
    """

    def __init__(self, params: ConvergenceParams | None = None):
        self.params = params or ConvergenceParams()
        self._t: list[float] = []
        self._d: list[float] = []
        self._passes = 0

    def precision_gate(self) -> float:
        """The ``SEM/mean`` precision gate, ``drift_tol/precision_p``.  Fixed (not size-scaled): the
        autocorrelation-corrected SEM is already size-aware because ``sigma/mean ~ 1/sqrt(N_atoms)``."""
        return self.params.drift_tol / self.params.precision_p

    def add_samples(self, times, densities):
        """Append one chunk's density time series (parallel sequences of equal length)."""
        times = list(times)
        densities = list(densities)
        if len(times) != len(densities):
            raise ValueError('times and densities must have equal length')
        self._t.extend(float(x) for x in times)
        self._d.extend(float(x) for x in densities)

    def check(self) -> ConvergenceReport:
        """Assess convergence over all accumulated samples.  Advances/resets the hysteresis counter
        and returns a :class:`ConvergenceReport`."""
        p = self.params
        t = np.array(self._t, dtype=float)
        d = np.array(self._d, dtype=float)
        last_step = float(t[-1]) if t.size else None

        # NaN / inf density -> the box blew up (NAMD produced garbage).  Abort, do not converge.
        if d.size and not np.all(np.isfinite(d)):
            self._passes = 0
            return ConvergenceReport(blowup=True, last_step=last_step,
                                     reason='non-finite density (box blowup)')

        # Burn-in: drop samples during the barostat transient.  Measured *relative to the first sample*
        # (the start of this task's NPT), not as an absolute timestep: density_equilibrate continues
        # from a `firsttimestep` that is already well past `burn_in` (e.g. 2400/8000 after the NVT
        # warm-up), so an absolute `t > burn_in` would drop nothing and leave the transient in.  Then
        # assess only the *trailing* `window_frac` of what remains, so the early densification ramp ages
        # out of the window as the run lengthens -- otherwise the ancient transient keeps inflating the
        # drift/SEM forever and a plateaued box never converges (observed: a full-history window on a
        # BPTI box that plateaued by ~30k steps still read drift ~5e-3 at 80k).  The window still
        # *grows* in absolute length as the run grows, so the precision gate keeps tightening.
        origin = float(t[0]) if t.size else 0.0
        keep = t > origin + p.burn_in
        tw, dw = t[keep], d[keep]
        if 0.0 < p.window_frac < 1.0 and tw.size:
            start = int(tw.size * (1.0 - p.window_frac))
            tw, dw = tw[start:], dw[start:]

        # Need a few windowed samples just to fit a line and estimate its variance (the real gate is
        # the autocorrelation-reliability guard below, which needs many more).
        if tw.size < 4:
            self._passes = 0
            return ConvergenceReport(passes=0, n_window=int(tw.size), last_step=last_step,
                                     reason=f'insufficient window ({tw.size} < 4 samples)')

        mean = float(dw.mean())

        # Honest SEM: the .xst frames are autocorrelated (tau_int ~ hundreds of steps for NPT density),
        # so the standard error of the mean is sigma/sqrt(N_eff) with N_eff = N_window/tau_int, not
        # sigma/sqrt(N).  We also refuse to trust the estimate until the window spans several
        # correlation times (tau_int is not estimable from a short series) -- until then the mean is
        # simply not resolved, no matter how flat it looks.
        tau_int = integrated_autocorr_time(dw)
        n_eff = tw.size / tau_int if tau_int > 0 else float(tw.size)
        sigma = float(np.std(dw, ddof=1))
        sem = sigma / np.sqrt(n_eff) if n_eff >= 1.0 else float('inf')
        sem_over_mean = sem / mean if mean else float('inf')
        gate = self.precision_gate()
        # Reliability: the window must span several correlation times AND hold an absolute minimum of
        # samples.  The absolute floor matters because tau_int is *under*-estimated from a short series
        # (you cannot see a long correlation with few samples), which would otherwise let a too-short
        # window pass the tau-relative test spuriously.
        reliable = (tw.size >= _MIN_WINDOW_SAMPLES
                    and tw.size >= p.autocorr_reliability * tau_int)
        precision_met = reliable and (sem_over_mean < gate)

        # Drift: fit a line to the raw windowed density and test whether we can *confidently* say the
        # fractional change over the window is below tolerance.  We use the UPPER confidence bound on
        # |drift| = |slope|*span/mean + drift_conf*SE, not the point estimate, so a noisy slope on a
        # plateau (few independent samples) raises its own uncertainty and correctly fails to certify --
        # instead of the raw point estimate randomly tripping above/below the tolerance and resetting the
        # hysteresis counter.  The slope's standard error is inflated by sqrt(tau_int): OLS assumes
        # independent residuals, but the density is autocorrelated, so the honest slope SE is larger.
        tc = tw - tw.mean()
        Sxx = float(np.dot(tc, tc))
        if Sxx > 0 and tw.size > 2:
            slope = float(np.dot(tc, dw - mean) / Sxx)
            resid = dw - (mean + slope * tc)
            s2 = float(np.dot(resid, resid) / (tw.size - 2))
            se_slope = np.sqrt(s2 / Sxx) * np.sqrt(max(tau_int, 1.0))  # autocorrelation-inflated
        else:
            slope, se_slope = 0.0, float('inf')
        span = float(tw[-1] - tw[0])
        drift = abs(slope) * span / mean if mean else float('inf')
        drift_se = se_slope * span / mean if mean else float('inf')
        drift_hi = drift + p.drift_conf * drift_se       # upper confidence bound on |drift|
        drift_ok = drift_hi < p.drift_tol

        below_min = last_step is None or last_step < p.min_steps
        this_pass = precision_met and drift_ok and not below_min
        self._passes = self._passes + 1 if this_pass else 0
        converged = self._passes >= p.n_consecutive

        if below_min:
            reason = f'below min_steps ({last_step:.0f} < {p.min_steps})'
        elif not reliable:
            reason = (f'window too short for its autocorrelation (N {tw.size}, need '
                      f'>= {_MIN_WINDOW_SAMPLES} and >= {p.autocorr_reliability:g}*tau {tau_int:.0f}); '
                      f'mean not yet resolved')
        elif not precision_met:
            reason = (f'precision gate unmet (SEM/mean {sem_over_mean:.2e} >= {gate:.2e}; '
                      f'tau {tau_int:.0f} fr, N_eff {n_eff:.0f})')
        elif not drift_ok:
            reason = (f'drift not confidently below tol (|drift| {drift:.2e}, '
                      f'{p.drift_conf:g}-sigma upper {drift_hi:.2e} >= {p.drift_tol:.2e})')
        elif converged:
            reason = f'converged: drift confidently < tol ({drift_hi:.2e} < {p.drift_tol:.2e}) x{self._passes}'
        else:
            reason = (f'passing ({self._passes}/{p.n_consecutive}); drift upper {drift_hi:.2e} '
                      f'< tol {p.drift_tol:.2e}')

        return ConvergenceReport(converged=converged, passes=self._passes, drift=drift,
                                 drift_hi=drift_hi, sem_over_mean=sem_over_mean,
                                 precision_met=precision_met, mean_density=mean,
                                 n_window=int(tw.size), tau_int=tau_int, n_eff=n_eff,
                                 last_step=last_step, reason=reason)


@dataclass
class JointConvergenceReport:
    """Verdict from one :meth:`JointConvergence.check` -- the per-observable reports plus the joint
    decision (both observables must be stationary together for ``n_consecutive`` checks)."""
    converged: bool = False
    passes: int = 0
    blowup: bool = False
    reports: dict = field(default_factory=dict)  #: name -> ConvergenceReport
    last_step: float | None = None
    reason: str = ''


class JointConvergence:
    """Track several named observables (e.g. ``density`` and ``area``) and converge only when **all**
    of them are stationary *simultaneously* for ``n_consecutive`` successive checks.

    Each observable gets its own :class:`DensityConvergenceMonitor` (so it may have its own tolerances,
    tau, and window), configured with ``n_consecutive = 1`` so that a monitor's ``check().converged``
    reports whether *this* check passed; this class owns the joint hysteresis counter.  A non-finite
    value in any observable is a blow-up.  Usage mirrors the single monitor::

        jc = JointConvergence({'density': mon_d, 'area': mon_a}, n_consecutive=3)
        for chunk:
            jc.add_samples('density', t, dens); jc.add_samples('area', t, area)
            rep = jc.check()
            if rep.blowup: abort
            if rep.converged: break
    """

    def __init__(self, monitors: dict, n_consecutive: int = 3):
        if not monitors:
            raise ValueError('JointConvergence needs at least one observable monitor')
        self.monitors = monitors
        for m in self.monitors.values():
            m.params.n_consecutive = 1   # this class owns the joint hysteresis
        self.n_consecutive = int(n_consecutive)
        self._passes = 0

    def add_samples(self, name, times, values):
        """Append one chunk's samples for observable ``name``."""
        self.monitors[name].add_samples(times, values)

    def check(self) -> JointConvergenceReport:
        """Assess all observables; advance/reset the joint hysteresis counter."""
        reports = {name: m.check() for name, m in self.monitors.items()}
        last_step = next((r.last_step for r in reports.values() if r.last_step is not None), None)
        if any(r.blowup for r in reports.values()):
            self._passes = 0
            blown = [n for n, r in reports.items() if r.blowup]
            return JointConvergenceReport(blowup=True, reports=reports, last_step=last_step,
                                          reason=f'non-finite {"/".join(blown)} (blowup)')
        joint_pass = all(r.converged for r in reports.values())   # each monitor n_consecutive == 1
        self._passes = self._passes + 1 if joint_pass else 0
        converged = self._passes >= self.n_consecutive
        if converged:
            reason = f'all observables converged for {self._passes} consecutive checks'
        elif joint_pass:
            reason = f'all passing ({self._passes}/{self.n_consecutive})'
        else:
            waiting = [n for n, r in reports.items() if not r.converged]
            reason = f'waiting on {", ".join(waiting)}'
        return JointConvergenceReport(converged=converged, passes=self._passes, reports=reports,
                                      last_step=last_step, reason=reason)
