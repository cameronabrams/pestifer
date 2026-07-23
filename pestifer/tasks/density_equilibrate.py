# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`DensityEquilibrateTask` -- a self-terminating replacement for the
hand-written ladder of progressively longer NPT runs that used to end every solvated build.

The task runs NPT as a series of **stability-bounded short restarts** (each a NAMD ``run`` that
recomputes the patch/PME grid before the shrinking cell outruns the patch ``margin``) and **stops
itself** when the box density has converged, or when a hard step ceiling is hit.  See
``docs/design/density-equilibrate.md`` and :mod:`pestifer.util.density_convergence` for the design
and the (deliberately NAMD-free, unit-tested) convergence machinery.

It automates both jobs the ladder did by hand: (1) staging the restarts -- chunk length is a
*stability* constraint, sized adaptively from the observed shrink rate (short early, longer late);
(2) deciding when enough is enough -- a practical fractional-drift tolerance made size-independent by
a precision gate, with ``n_consecutive`` hysteresis.
"""
import logging
import os

from .mdtask import MDTask
from ..core.artifacts import DataFileArtifact, PNGImageFileArtifact
from ..core.errors import PestiferBuildError
from ..util.density_convergence import (
    ConvergenceParams,
    DensityConvergenceMonitor,
    next_chunk_steps,
    total_mass_amu,
    volume_to_density,
    xst_cell_volumes,
    xst_max_shrink_rate,
)

logger = logging.getLogger(__name__)


class DensityEquilibrateTask(MDTask):
    """Convergence-driven NPT box-density equilibration (self-terminating NPT ladder)."""

    _yaml_header = 'density_equilibrate'

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, MD_OUTPUT
        return TaskContract(requires=(STATE,), provides=(STATE, MD_OUTPUT))

    def provision(self, packet: dict):
        super().provision(packet)
        # This task *is* NPT by definition; message it as such regardless of any stray spec.
        self.extra_message = 'ensemble: NPT (density-convergence)'

    def do(self):
        """Run stability-bounded NPT chunks until the box density converges or ``max_steps`` is
        reached; register the final state and a convergence report + density-vs-time plot."""
        specs = self.specs
        specs['ensemble'] = 'NPT'  # force -- this task is defined as NPT
        margin = float(specs['margin'])
        shrink_safety = float(specs['shrink_safety'])
        chunk_min = int(specs['chunk_min'])
        chunk_max = int(specs['chunk_max'])
        max_steps = int(specs['max_steps'])
        min_steps = int(specs['min_steps'])

        state = self.get_current_artifact('state')
        if not state or not state.psf:
            raise PestiferBuildError(f'{self.taskname}: no state PSF to compute density from')
        mass_amu = total_mass_amu(state.psf.path)

        monitor = DensityConvergenceMonitor(ConvergenceParams(
            drift_tol=float(specs['drift_tol']),
            precision_p=float(specs['precision_p']),
            n_blocks=int(specs['n_blocks']),
            burn_in=int(specs['burn_in']),
            min_steps=min_steps,
            n_consecutive=int(specs['n_consecutive']),
        ))

        rows = []             # per-chunk diagnostics for the convergence report
        all_t, all_d = [], []  # full density-vs-time series for the plot
        total_steps = self.get_current_artifact_data('firsttimestep') or 0
        chunk = chunk_min      # conservative first chunk: the box shrinks fastest now
        stop_reason = None
        n_chunk = 0

        while total_steps < max_steps:
            n_chunk += 1
            # Do not overshoot the ceiling.
            this_chunk = min(chunk, max_steps - total_steps)
            specs['nsteps'] = this_chunk
            logger.info(f'{self.taskname}: NPT chunk {n_chunk} -- {this_chunk} steps '
                        f'(firsttimestep {total_steps})')
            rc = self.namdrun()
            if rc != 0:
                raise PestiferBuildError(f'{self.taskname}: NAMD chunk {n_chunk} failed (rc={rc})')
            total_steps += this_chunk

            xst = f'{self.basename}.xst'
            if not os.path.exists(xst):
                raise PestiferBuildError(f'{self.taskname}: chunk {n_chunk} produced no .xst ({xst})')
            ts, vols = xst_cell_volumes(xst)
            if ts.size == 0:
                logger.warning(f'{self.taskname}: chunk {n_chunk} .xst had no frames; '
                               f'continuing (raise xstfreq resolution if this persists)')
            else:
                dens = volume_to_density(vols, mass_amu)
                monitor.add_samples(ts, dens)
                all_t.extend(ts.tolist())
                all_d.extend(dens.tolist())

            report = monitor.check()
            rows.append((n_chunk, total_steps, this_chunk, report))
            logger.info(f'{self.taskname}: chunk {n_chunk} @ step {total_steps} -- '
                        f'rho={_fmt(report.mean_density)} g/cc, drift={_fmt(report.drift)}, '
                        f'SEM/mean={_fmt(report.sem_over_mean)} -- {report.reason}')

            if report.blowup:
                self._write_report(rows, mass_amu, stop_reason='ABORTED: box blowup (non-finite density)')
                raise PestiferBuildError(f'{self.taskname}: box blew up (non-finite density) at '
                                         f'step {total_steps}; see convergence report')

            if report.converged:
                stop_reason = (f'CONVERGED at step {total_steps}: fractional drift {_fmt(report.drift)} '
                               f'< tol {specs["drift_tol"]:g} for {report.passes} consecutive checks')
                logger.info(f'{self.taskname}: {stop_reason}')
                break

            # STABILITY: size the next chunk from the shrink rate just observed so the projected
            # per-dimension cell shrink stays below shrink_safety * margin within the next run.
            rate = xst_max_shrink_rate(xst)
            chunk = next_chunk_steps(rate, margin, shrink_safety, chunk_min, chunk_max)
            logger.debug(f'{self.taskname}: observed shrink rate {rate:.3e} A/step -> next chunk {chunk} steps')
        else:
            last = rows[-1][3] if rows else None
            resid = _fmt(last.drift) if last else 'n/a'
            gate = 'met' if (last and last.precision_met) else 'UNMET'
            stop_reason = (f'CEILING: reached max_steps ({max_steps}) without convergence; '
                           f'residual drift {resid}, precision gate {gate} -- system may not have settled')
            logger.warning(f'{self.taskname}: {stop_reason}')

        self._write_report(rows, mass_amu, stop_reason=stop_reason)
        self._write_plot(all_t, all_d, rows)
        return 0

    def _write_report(self, rows, mass_amu, stop_reason):
        """Write a per-chunk convergence report (``<basename>-density.dat``) and register it."""
        fn = f'{self.basename}-density.dat'
        with open(fn, 'w') as f:
            f.write(f'# density_equilibrate convergence report -- {self.taskname}\n')
            f.write(f'# total system mass: {mass_amu:.1f} amu\n')
            f.write(f'# stop: {stop_reason}\n')
            f.write('# chunk  step  nsteps  rho[g/cc]  drift  SEM/mean  precision  passes  reason\n')
            for n, step, nsteps, r in rows:
                f.write(f'{n:6d}  {step:8d}  {nsteps:6d}  {_fmt(r.mean_density):>9}  '
                        f'{_fmt(r.drift):>9}  {_fmt(r.sem_over_mean):>9}  '
                        f'{"yes" if r.precision_met else "no":>9}  {r.passes:6d}  {r.reason}\n')
        self.register(f'{self.basename}-density', key='density_report', artifact_type=DataFileArtifact)
        logger.info(f'{self.taskname}: convergence report -> {fn}')

    def _write_plot(self, times, densities, rows):
        """Write a density-vs-time PNG (``<basename>-density.png``) and register it."""
        if not times:
            return
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover - plotting is best-effort
            logger.warning(f'{self.taskname}: could not import matplotlib for density plot ({e})')
            return
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(times, densities, lw=0.8, color='#3366aa')
        # burn-in marker and per-chunk running mean
        burn = float(self.specs['burn_in'])
        if burn > 0 and burn < max(times):
            ax.axvline(burn, ls=':', color='0.6', label=f'burn-in ({burn:g})')
        means = [(step, r.mean_density) for (_n, step, _ns, r) in rows if r.mean_density is not None]
        if means:
            mx, my = zip(*means)
            ax.plot(mx, my, 'o-', color='#cc3333', ms=3, lw=1, label='window mean')
        ax.set_xlabel('timestep')
        ax.set_ylabel('density (g/cc)')
        ax.set_title(f'{self.taskname}: box density vs. time')
        ax.legend(fontsize=8, loc='best')
        fig.tight_layout()
        fn = f'{self.basename}-density.png'
        fig.savefig(fn, dpi=120)
        plt.close(fig)
        self.register(fn, key='density_plot', artifact_type=PNGImageFileArtifact, keep=True)
        logger.info(f'{self.taskname}: density plot -> {fn}')


def _fmt(x):
    """Compact fixed/scientific formatting that tolerates ``None``."""
    if x is None:
        return 'n/a'
    ax = abs(x)
    if ax != 0 and (ax < 1e-3 or ax >= 1e4):
        return f'{x:.2e}'
    return f'{x:.4f}'
