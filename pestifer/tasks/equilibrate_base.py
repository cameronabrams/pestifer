# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Shared base for the self-terminating, chunked equilibration tasks
(:class:`~pestifer.tasks.density_equilibrate.DensityEquilibrateTask` and
:class:`~pestifer.tasks.membrane_equilibrate.MembraneEquilibrateTask`).

:class:`ChunkedEquilibrateTask` owns the parts that are identical regardless of *what* observable is
being converged: the stability-bounded chunk loop, the reactive patch-grid crash-and-retry, the
shrink-rate-driven next-chunk sizing, and the ``max_steps`` ceiling.  Everything observable-specific --
which ensemble to run, what to read from the ``.xst``, how to decide convergence, and what report/plot
to write -- is delegated to hooks a subclass implements.  See ``docs/design/density-equilibrate.md``
and ``docs/design/membrane-equilibrate.md``.
"""
import logging
import os

from .mdtask import MDTask
from ..core.errors import PestiferBuildError
from ..util.density_convergence import (
    is_patch_grid_crash,
    next_chunk_steps,
    parse_patch_grid,
    quantize_steps,
    xst_max_shrink_rate,
)

logger = logging.getLogger(__name__)

# Hard floor for crash-recovery chunk shrinking: if even this many steps outruns the patch grid, the
# box is at the grid limit and no shorter chunk helps -- the fix is a larger `margin`, so we stop
# retrying and surface that.  Independent of `chunk_min` (the nominal first-chunk length).
_MIN_RETRY_CHUNK = 20


def _read_text(path):
    """Read a text file, returning '' if it is missing/unreadable (e.g. a crash before any output)."""
    try:
        with open(path, errors='replace') as f:
            return f.read()
    except OSError:
        return ''


def _fmt(x):
    """Compact fixed/scientific formatting that tolerates ``None``."""
    if x is None:
        return 'n/a'
    ax = abs(x)
    if ax != 0 and (ax < 1e-3 or ax >= 1e4):
        return f'{x:.2e}'
    return f'{x:.4f}'


class ChunkedEquilibrateTask(MDTask):
    """Run an ensemble in stability-bounded chunks until a subclass-defined observable set converges.

    Subclass contract (all observable-specific behavior):

    * ``_ensemble`` / ``_provision_message`` -- class attributes.
    * ``_setup(state, min_steps)`` -- build the convergence tracker + per-observable context, log the
      gate; may raise.  Initialize ``self._rows`` accumulation here or rely on the base (it sets
      ``self._rows = []``).
    * ``_ingest_chunk(xst, total_steps, this_chunk, n_chunk) -> (blowup, converged)`` -- read the
      chunk's observables from ``xst``, update the tracker, accumulate plot series, append a row to
      ``self._rows``, and log the per-chunk line.
    * ``_converged_stop_reason(total_steps)`` / ``_ceiling_stop_reason(total_steps, max_steps)`` /
      ``_blowup_stop_reason(total_steps)`` / ``_blowup_error(total_steps)`` -- human-readable strings.
    * ``_write_outputs(stop_reason)`` -- write the report + plot artifacts.
    """

    _ensemble = 'NPT'
    _provision_message = 'ensemble: NPT (convergence-driven)'

    def provision(self, packet: dict):
        super().provision(packet)
        self.extra_message = self._provision_message

    def do(self):
        specs = self.specs
        specs['ensemble'] = self._ensemble  # force -- the ensemble is part of the task's definition
        margin = float(specs['margin'])
        shrink_safety = float(specs['shrink_safety'])
        chunk_min = int(specs['chunk_min'])
        chunk_max = int(specs['chunk_max'])
        max_steps = int(specs['max_steps'])
        min_steps = int(specs['min_steps'])
        chunk_growth = float(specs['chunk_growth'])
        max_shrink_retries = int(specs['max_shrink_retries'])
        # NAMD requires every `run` to be a whole number of cycles; quantize all chunk lengths to it.
        steps_per_cycle = int(self.namd_global_config.get('solvated', {}).get('stepspercycle', 10))

        state = self.get_current_artifact('state')
        if not state or not state.psf:
            raise PestiferBuildError(f'{self.taskname}: no state PSF to compute observables from')
        self._rows = []
        self._setup(state, min_steps)

        start_step = self.get_current_artifact_data('firsttimestep') or 0
        total_steps = start_step
        # max_steps is a per-run BUDGET (steps THIS equilibration may run), not an absolute step ceiling.
        # A membrane_equilibrate whose firsttimestep is already large -- e.g. an asymmetric build's second
        # calibration patch, running after the first on the same accumulating subcontroller step counter
        # -- must still get its full budget, not be starved by the inherited offset (which had it ceiling
        # after a few thousand steps and report an un-relaxed area). Absolute ~ budget when firsttimestep
        # is small (a fresh post-embed run), so this changes nothing there.
        step_ceiling = start_step + max_steps
        chunk = chunk_min      # conservative first chunk: the box shrinks fastest now
        stop_reason = None
        n_chunk = 0
        shrink_retries = 0     # consecutive patch-grid retries at the current boundary
        grid_logged = False

        while total_steps < step_ceiling:
            # Do not overshoot the ceiling; quantize to a whole number of NAMD cycles.  If less than one
            # cycle remains in the budget, we cannot run a valid chunk -- stop here.
            if step_ceiling - total_steps < steps_per_cycle:
                break
            this_chunk = quantize_steps(min(chunk, step_ceiling - total_steps), steps_per_cycle,
                                        minimum=_MIN_RETRY_CHUNK)
            specs['nsteps'] = this_chunk
            attempt = f' (retry {shrink_retries}/{max_shrink_retries})' if shrink_retries else ''
            logger.info(f'{self.taskname}: {self._ensemble} chunk {n_chunk + 1} -- {this_chunk} steps '
                        f'(firsttimestep {total_steps}){attempt}')
            try:
                self.namdrun()
            except PestiferBuildError:
                # Reactive safety net.  If NAMD died because the shrinking cell outran its
                # (start-of-run) patch grid, roll back and retry with a shorter chunk: a crashed namdrun
                # registers no new `state`, so `get_current_artifact('state')` still points at the
                # *previous* completed chunk -- the safe, pre-threshold rollback point -- and the next
                # namdrun() continues from it.  We simply run fewer steps so we stop before the cell
                # shrinks past the grid.  NAMD's own abort is the just-in-time trigger; the threshold
                # need not be modeled exactly.  Any other failure is re-raised unchanged.
                log_text = _read_text(f'{self.basename}.log')
                if (is_patch_grid_crash(log_text) and shrink_retries < max_shrink_retries
                        and this_chunk > _MIN_RETRY_CHUNK):
                    shrink_retries += 1
                    chunk = max(_MIN_RETRY_CHUNK, this_chunk // 2)
                    logger.warning(
                        f'{self.taskname}: {self._ensemble} chunk {n_chunk + 1} outran the patch grid '
                        f'at {this_chunk} steps; rolling back to the last completed state and retrying '
                        f'with {chunk} steps (retry {shrink_retries}/{max_shrink_retries})')
                    continue
                if is_patch_grid_crash(log_text):
                    raise PestiferBuildError(
                        f'{self.taskname}: NAMD keeps outrunning the patch grid even at '
                        f'{this_chunk} steps after {shrink_retries} retries. The freshly solvated '
                        f'cell is shrinking faster than the patch margin allows; increase `margin` '
                        f'(currently {margin:g} A) in this task or the `namd` config, or pre-shrink '
                        f'the box.')
                raise
            shrink_retries = 0
            n_chunk += 1
            total_steps += this_chunk

            xst = f'{self.basename}.xst'
            if not os.path.exists(xst):
                raise PestiferBuildError(f'{self.taskname}: chunk {n_chunk} produced no .xst ({xst})')
            if not grid_logged:
                self._log_patch_grid(f'{self.basename}.log')
                grid_logged = True

            blowup, converged = self._ingest_chunk(xst, total_steps, this_chunk, n_chunk)
            if blowup:
                self._write_outputs(self._blowup_stop_reason(total_steps))
                raise PestiferBuildError(self._blowup_error(total_steps))
            if converged:
                stop_reason = self._converged_stop_reason(total_steps)
                logger.info(f'{self.taskname}: {stop_reason}')
                break

            # STABILITY: size the next chunk from the shrink rate just observed so the projected
            # per-dimension cell shrink stays below shrink_safety * margin within the next run.  The
            # floor is _MIN_RETRY_CHUNK, not chunk_min, so a size learned by crash-recovery is not yanked
            # back up.  Growth is capped at chunk_growth x the chunk we just ran (AIMD-style: halve on a
            # crash, re-grow gently) so we do not leap straight back to a size that just crashed --
            # important on CPU, where the true margin can be far below our estimate.
            rate = xst_max_shrink_rate(xst)
            rate_sized = next_chunk_steps(rate, margin, shrink_safety, _MIN_RETRY_CHUNK, chunk_max)
            grow_cap = max(_MIN_RETRY_CHUNK, int(this_chunk * chunk_growth))
            chunk = min(rate_sized, grow_cap, chunk_max)
            logger.debug(f'{self.taskname}: shrink rate {rate:.3e} A/step -> rate-sized {rate_sized}, '
                         f'grow-cap {grow_cap} -> next chunk {chunk} steps')

        # No convergence/blowup break -> we exhausted the budget (either total_steps >= step_ceiling or
        # too few steps remained for a whole cycle).
        if stop_reason is None:
            stop_reason = self._ceiling_stop_reason(total_steps, max_steps)
            logger.warning(f'{self.taskname}: {stop_reason}')

        # expose the outcome so a caller (e.g. make_membrane_system's calibration) can trust this
        # self-terminating equilibration's own verdict instead of re-deriving convergence downstream
        self.converged = stop_reason.startswith('CONVERGED')
        self._write_outputs(stop_reason)
        return 0

    def _log_patch_grid(self, log_path):
        """Report the patch grid NAMD actually built (observability only; correctness is guaranteed by
        the crash-and-retry, not by this estimate)."""
        info = parse_patch_grid(_read_text(log_path))
        if not info:
            return
        px, py, pz = info['patches']
        msg = f'{self.taskname}: NAMD patch grid {px}x{py}x{pz}'
        if 'cell_edges' in info:
            ex, ey, ez = info['cell_edges']
            msg += f'; cell {ex:.1f}x{ey:.1f}x{ez:.1f} A'
        if 'min_patch_dim' in info:
            msg += f'; min patch dim {info["min_patch_dim"]:.1f} A'
        if 'shrink_headroom' in info:
            msg += (f'; ~{info["shrink_headroom"]:.1f} A shrink headroom to pairlist '
                    f'(approx; crash-and-retry is the real guard)')
        logger.info(msg)

    # --- hooks a subclass must implement -------------------------------------------------------
    def _setup(self, state, min_steps):
        raise NotImplementedError

    def _ingest_chunk(self, xst, total_steps, this_chunk, n_chunk):
        raise NotImplementedError

    def _converged_stop_reason(self, total_steps):
        raise NotImplementedError

    def _ceiling_stop_reason(self, total_steps, max_steps):
        raise NotImplementedError

    def _blowup_stop_reason(self, total_steps):
        raise NotImplementedError

    def _blowup_error(self, total_steps):
        raise NotImplementedError

    def _write_outputs(self, stop_reason):
        raise NotImplementedError
