# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`DensityEquilibrateTask` -- a self-terminating replacement for the
hand-written ladder of progressively longer NPT runs that used to end every solvated build.

The task runs NPT as a series of **stability-bounded short restarts** (each a NAMD ``run`` that
recomputes the patch/PME grid before the shrinking cell outruns the patch ``margin``) and **stops
itself** when the box density has converged, or when a hard step ceiling is hit.  The chunk loop,
crash-and-retry, and chunk sizing live in :class:`~pestifer.tasks.equilibrate_base.ChunkedEquilibrateTask`;
this class supplies the density-specific hooks (read density from the ``.xst`` cell volume + PSF mass,
track it with the autocorrelation-corrected :class:`DensityConvergenceMonitor`, and write the density
report + plot).  See ``docs/design/density-equilibrate.md`` and :mod:`pestifer.util.density_convergence`.
"""
import logging
import os

from .equilibrate_base import ChunkedEquilibrateTask, _fmt
from ..core.artifacts import DataFileArtifact, PNGImageFileArtifact
from ..util.density_convergence import (
    ConvergenceParams,
    DensityConvergenceMonitor,
    total_atoms,
    total_mass_amu,
    volume_to_density,
    xst_cell_volumes,
)

logger = logging.getLogger(__name__)


class DensityEquilibrateTask(ChunkedEquilibrateTask):
    """Convergence-driven NPT box-density equilibration (self-terminating NPT ladder)."""

    _yaml_header = 'density_equilibrate'
    _ensemble = 'NPT'
    _provision_message = 'ensemble: NPT (density-convergence)'

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, MD_OUTPUT
        return TaskContract(requires=(STATE,), provides=(STATE, MD_OUTPUT))

    # --- ChunkedEquilibrateTask hooks --------------------------------------------------------------
    def _setup(self, state, min_steps):
        specs = self.specs
        self._mass_amu = total_mass_amu(state.psf.path)
        n_atoms = total_atoms(state.psf.path)
        self._monitor = DensityConvergenceMonitor(ConvergenceParams(
            drift_tol=float(specs['drift_tol']),
            precision_p=float(specs['precision_p']),
            drift_conf=float(specs['drift_conf']),
            burn_in=int(specs['burn_in']),
            window_frac=float(specs['window_frac']),
            min_steps=min_steps,
            n_consecutive=int(specs['n_consecutive']),
            autocorr_reliability=float(specs['autocorr_reliability']),
        ))
        self._all_t, self._all_d = [], []   # full density-vs-time series for the plot
        logger.info(f'{self.taskname}: {n_atoms} atoms; precision gate SEM/mean < '
                    f'{self._monitor.precision_gate():.2e} (autocorrelation-corrected)')

    def _ingest_chunk(self, xst, total_steps, this_chunk, n_chunk):
        ts, vols = xst_cell_volumes(xst)
        if ts.size == 0:
            logger.warning(f'{self.taskname}: chunk {n_chunk} .xst had no frames; '
                           f'continuing (raise xstfreq resolution if this persists)')
        else:
            dens = volume_to_density(vols, self._mass_amu)
            self._monitor.add_samples(ts, dens)
            self._all_t.extend(ts.tolist())
            self._all_d.extend(dens.tolist())
        report = self._monitor.check()
        self._rows.append((n_chunk, total_steps, this_chunk, report))
        # signed drift (negative = de-densifying) so the trend direction is visible; the
        # convergence gate keys off the magnitude (report.drift).
        logger.info(f'{self.taskname}: chunk {n_chunk} @ step {total_steps} -- '
                    f'rho={_fmt(report.mean_density)} g/cc, drift={_fmt(report.signed_drift)}, '
                    f'SEM/mean={_fmt(report.sem_over_mean)} -- {report.reason}')
        return report.blowup, report.converged

    def _converged_stop_reason(self, total_steps):
        r = self._rows[-1][3]
        return (f'CONVERGED at step {total_steps}: fractional drift {_fmt(r.signed_drift)} '
                f'(|.| < tol {self.specs["drift_tol"]:g}) for {r.passes} consecutive checks')

    def _ceiling_stop_reason(self, total_steps, max_steps):
        last = self._rows[-1][3] if self._rows else None
        resid = _fmt(last.signed_drift) if last else 'n/a'
        gate = 'met' if (last and last.precision_met) else 'UNMET'
        return (f'CEILING: reached max_steps ({max_steps}) without convergence; '
                f'residual drift {resid}, precision gate {gate} -- system may not have settled')

    def _blowup_stop_reason(self, total_steps):
        return 'ABORTED: box blowup (non-finite density)'

    def _blowup_error(self, total_steps):
        return (f'{self.taskname}: box blew up (non-finite density) at '
                f'step {total_steps}; see convergence report')

    def _write_outputs(self, stop_reason):
        self._write_report(self._rows, self._mass_amu, stop_reason=stop_reason)
        self._write_plot(self._all_t, self._all_d, self._rows, stop_reason=stop_reason)

    # --- report + plot -----------------------------------------------------------------------------
    def _write_report(self, rows, mass_amu, stop_reason):
        """Write a per-chunk convergence report (``<basename>-density.dat``) and register it."""
        fn = f'{self.basename}-density.dat'
        with open(fn, 'w') as f:
            f.write(f'# density_equilibrate convergence report -- {self.taskname}\n')
            f.write(f'# total system mass: {mass_amu:.1f} amu\n')
            f.write(f'# stop: {stop_reason}\n')
            f.write('# chunk  step  nsteps  rho[g/cc]  drift  drift_hi  SEM/mean  tau[fr]  N_eff  '
                    'precision  passes  reason\n')
            for n, step, nsteps, r in rows:
                f.write(f'{n:6d}  {step:8d}  {nsteps:6d}  {_fmt(r.mean_density):>9}  '
                        f'{_fmt(r.signed_drift):>9}  {_fmt(r.drift_hi):>9}  {_fmt(r.sem_over_mean):>9}  '
                        f'{_fmt(r.tau_int):>7}  {_fmt(r.n_eff):>6}  '
                        f'{"yes" if r.precision_met else "no":>9}  {r.passes:6d}  {r.reason}\n')
        self.register(f'{self.basename}-density', key='density_report', artifact_type=DataFileArtifact)
        logger.info(f'{self.taskname}: convergence report -> {fn}')

    def _write_plot(self, times, densities, rows, stop_reason=''):
        """Write a density-vs-time PNG (``<basename>-density.png``) and register it.

        Shows the whole NPT densification from the task's first sample (not padded back to timestep 0,
        where no data exists -- the minimize/NVT warm-up are separate tasks), shades the burn-in
        transient that is discarded and the trailing window that is actually assessed, overlays the
        per-chunk window-mean density, and marks the convergence step when reached."""
        if not times:
            return
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover - plotting is best-effort
            logger.warning(f'{self.taskname}: could not import matplotlib for density plot ({e})')
            return
        import numpy as np
        t = np.asarray(times, dtype=float)
        d = np.asarray(densities, dtype=float)
        t0, t1 = float(t.min()), float(t.max())

        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.plot(t, d, lw=0.8, color='#3366aa', label='density (per frame)')

        # Burn-in region: discarded transient at the *start of this task's NPT* (relative to the first
        # sample), matching the monitor's `t > origin + burn_in`.
        burn = float(self.specs['burn_in'])
        window_frac = float(self.specs['window_frac'])
        if burn > 0 and t0 + burn < t1:
            ax.axvspan(t0, t0 + burn, color='0.85', label=f'burn-in ({burn:g} steps, discarded)')

        # Trailing window actually assessed at the final check: the last `window_frac` of the
        # post-burn-in samples.  Shade it so the discarded-ramp vs. assessed-plateau split is visible.
        post = t[t > t0 + burn]
        if window_frac and 0.0 < window_frac < 1.0 and post.size:
            win_start = float(post[int(post.size * (1.0 - window_frac))])
            ax.axvspan(win_start, t1, color='#ddeecc', alpha=0.6,
                       label=f'final window (trailing {window_frac:g})')

        means = [(step, r.mean_density) for (_n, step, _ns, r) in rows if r.mean_density is not None]
        if means:
            mx, my = zip(*means)
            ax.plot(mx, my, 'o-', color='#cc3333', ms=3, lw=1, label='window mean')
            ax.axhline(my[-1], ls='--', color='#cc3333', lw=0.8, alpha=0.7,
                       label=f'final rho = {my[-1]:.4f} g/cc')

        # Convergence marker.
        if stop_reason.startswith('CONVERGED'):
            conv_step = rows[-1][1] if rows else t1
            ax.axvline(conv_step, ls='-', color='#118811', lw=1.2, label=f'converged @ {conv_step:.0f}')

        ax.set_xlim(t0, t1)
        ax.set_xlabel('timestep')
        ax.set_ylabel('density (g/cc)')
        ax.set_title(f'{self.taskname}\nbox density vs. time', fontsize=9)
        ax.legend(fontsize=7, loc='lower right', framealpha=0.9)
        fig.tight_layout()
        # keep convergence plots out of the run-dir root, in the mdplots/ subdir the mdplot task
        # already uses (kept, not swept into the artifacts tarball, so they stay easy to eyeball)
        output_dir = self.specs.get('output_dir', 'mdplots')
        os.makedirs(output_dir, exist_ok=True)
        fn = os.path.join(output_dir, f'{self.basename}-density.png')
        fig.savefig(fn, dpi=120)
        plt.close(fig)
        self.register(fn, key='density_plot', artifact_type=PNGImageFileArtifact, keep=True)
        logger.info(f'{self.taskname}: density plot -> {fn}')
