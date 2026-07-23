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
   late (box settled) -- the ladder's shape, derived rather than guessed.

3. **The convergence criterion.**  A *practical fractional-drift* tolerance made size-independent by
   a **precision gate** (grow the window/duration until the mean is resolved above noise) plus
   ``n_consecutive`` hysteresis.  See the design doc for why the rejected alternatives (fluctuation
   threshold, ``trend/SEM`` significance) do not generalize.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from .densityprofile import AMU_PER_A3_TO_G_PER_CC

logger = logging.getLogger(__name__)


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


@dataclass
class ConvergenceParams:
    """Tunables for the density-convergence criterion (see the design doc's Parameters table)."""
    drift_tol: float = 1e-3      #: converged when fractional drift over the window < this (~0.1%)
    precision_p: float = 3.0     #: precision gate: require SEM/mean < drift_tol / precision_p
    n_blocks: int = 6            #: blocks the trailing window is averaged into
    burn_in: int = 2000          #: leading steps discarded before assessing the trend
    min_steps: int = 4000        #: never declare convergence before this many steps
    n_consecutive: int = 3       #: successive passing checks required before stopping (hysteresis)


@dataclass
class ConvergenceReport:
    """Verdict + diagnostics from one :meth:`DensityConvergenceMonitor.check` call."""
    converged: bool = False
    passes: int = 0              #: consecutive passing checks accumulated so far
    drift: float | None = None   #: fractional drift over the window (|slope|*span/mean)
    sem_over_mean: float | None = None
    precision_met: bool = False
    mean_density: float | None = None
    n_window: int = 0            #: samples in the post-burn-in window
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

        # Burn-in: keep only samples strictly after the barostat transient.
        keep = t > p.burn_in
        tw, dw = t[keep], d[keep]

        # Need enough post-burn-in samples to populate the blocks.
        if tw.size < p.n_blocks:
            self._passes = 0
            return ConvergenceReport(passes=0, n_window=int(tw.size), last_step=last_step,
                                     reason=f'insufficient window ({tw.size} < {p.n_blocks} samples)')

        # Contiguous block means (density and time).
        idx = np.array_split(np.arange(tw.size), p.n_blocks)
        block_t = np.array([tw[ix].mean() for ix in idx])
        block_d = np.array([dw[ix].mean() for ix in idx])

        mean = float(dw.mean())
        sem = float(np.std(block_d, ddof=1) / np.sqrt(p.n_blocks))
        sem_over_mean = sem / mean if mean else float('inf')
        precision_met = sem_over_mean < (p.drift_tol / p.precision_p)

        # Fractional drift from the line fit to the block means.
        slope = float(np.polyfit(block_t, block_d, 1)[0])
        span = float(block_t[-1] - block_t[0])
        drift = abs(slope) * span / mean if mean else float('inf')

        below_min = last_step is None or last_step < p.min_steps
        this_pass = precision_met and (drift < p.drift_tol) and not below_min
        self._passes = self._passes + 1 if this_pass else 0
        converged = self._passes >= p.n_consecutive

        if below_min:
            reason = f'below min_steps ({last_step:.0f} < {p.min_steps})'
        elif not precision_met:
            reason = (f'precision gate unmet (SEM/mean {sem_over_mean:.2e} '
                      f'>= {p.drift_tol / p.precision_p:.2e})')
        elif drift >= p.drift_tol:
            reason = f'drift {drift:.2e} >= tol {p.drift_tol:.2e}'
        elif converged:
            reason = f'converged: drift {drift:.2e} < tol for {self._passes} consecutive checks'
        else:
            reason = f'passing ({self._passes}/{p.n_consecutive}); drift {drift:.2e} < tol {p.drift_tol:.2e}'

        return ConvergenceReport(converged=converged, passes=self._passes, drift=drift,
                                 sem_over_mean=sem_over_mean, precision_met=precision_met,
                                 mean_density=mean, n_window=int(tw.size), last_step=last_step,
                                 reason=reason)
