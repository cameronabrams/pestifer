# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of :class:`MembraneEquilibrateTask` -- a self-terminating NPgT equilibration for membrane
systems, the anisotropic sibling of :class:`~pestifer.tasks.density_equilibrate.DensityEquilibrateTask`.

It runs NPgT (``useFlexibleCell``, tensionless) in the same stability-bounded chunks and stops when
**both** the box density **and** the membrane lateral area (``A = a_x·b_y``) have converged.  Density +
area fully pin the anisotropic cell (``V = A·c_z``), so the box thickness follows.  The chunk loop,
crash-and-retry, and sizing are inherited from
:class:`~pestifer.tasks.equilibrate_base.ChunkedEquilibrateTask`; this class supplies the membrane
hooks and tracks the observables with :class:`JointConvergence`.  See
``docs/design/membrane-equilibrate.md``.

**Two-stage protocol (default, `two_stage: true`).**  A freshly gridded bilayer is built at the correct
lateral area (~SAPL) but is under-dense in *z* (the leaflet z-reservation over-sizes the box height).
A single tensionless piston asked to fix that density *and* find the equilibrium area at once takes the
excess volume out of the lateral dimensions too, co-compacting the area down to a metastable gel-like
APL.  So we decouple the two: **stage 1** settles the density at *constant lateral area*
(``useConstantArea``) -- excess volume leaves from z only, the area stays at build; **stage 2** then
relaxes the area under tensionless NPgT (``useConstantRatio``) from the already-correct density, a small
move that no longer collapses.  Stage 1 gates on density alone (area is pinned); stage 2 gates on both.
Set ``two_stage: false`` for the legacy single-stage behavior (e.g. a post-embed stage already near
density equilibrium).

Area tolerances currently reuse the density convergence specs (area is a soft mode that may want its own
values once its autocorrelation is measured -- see the design doc); optional `area_*` specs override
per-observable.  Area-per-lipid (APL) is reported when `lipids_per_leaflet` is supplied.
"""
import logging

from .equilibrate_base import ChunkedEquilibrateTask, _fmt
from ..core.artifacts import DataFileArtifact, PNGImageFileArtifact
from ..util.density_convergence import (
    ConvergenceParams,
    DensityConvergenceMonitor,
    JointConvergence,
    membrane_leaflet_geometry,
    total_atoms,
    total_mass_amu,
    volume_to_density,
    xst_cell_areas,
    xst_cell_volumes,
)

logger = logging.getLogger(__name__)


class MembraneEquilibrateTask(ChunkedEquilibrateTask):
    """Convergence-driven NPgT membrane equilibration on box density + lateral area."""

    _yaml_header = 'membrane_equilibrate'
    _ensemble = 'NPgT'
    _provision_message = 'ensemble: NPgT (density+area convergence)'

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, MD_OUTPUT
        return TaskContract(requires=(STATE,), provides=(STATE, MD_OUTPUT))

    def _params_for(self, prefix, min_steps):
        """Build ConvergenceParams from the shared convergence specs, letting ``<prefix>_<key>`` specs
        (e.g. ``area_drift_tol``) override per-observable.  Area defaults to the density values."""
        specs = self.specs

        def g(key, cast):
            v = specs.get(f'{prefix}_{key}')
            return cast(v) if v is not None else cast(specs[key])
        return ConvergenceParams(
            drift_tol=g('drift_tol', float),
            precision_p=g('precision_p', float),
            drift_conf=g('drift_conf', float),
            burn_in=g('burn_in', int),
            window_frac=g('window_frac', float),
            min_steps=min_steps,
            n_consecutive=int(specs['n_consecutive']),
            autocorr_reliability=g('autocorr_reliability', float),
        )

    # --- barostat mode per stage -------------------------------------------------------------------
    def _set_barostat_mode(self, phase):
        """Point ``other_parameters`` at the barostat mode for a stage, overriding the ``membrane``
        NAMD block (``useFlexibleCell`` stays on; we only flip the lateral constraint):

        * **phase 1 -- constant area** (``useConstantArea yes``, ``useConstantRatio no``): x,y are
          frozen and only z relaxes, so the under-dense box removes its excess volume *from z alone*
          and the lateral area (already correct at build, ~SAPL) is not dragged down while density settles.
        * **phase 2 -- tensionless** (``useConstantArea no``, ``useConstantRatio yes``): the current
          NPgT area relaxation; density is already at target, so only the small lateral move remains.

        A user-supplied ``other_parameters`` is preserved and merged under these keys (lowercase, to
        override the ``membrane`` block's ``useconstantratio``)."""
        mode = ({'useconstantarea': 'yes', 'useconstantratio': 'no'} if phase == 1
                else {'useconstantarea': 'no', 'useconstantratio': 'yes'})
        self.specs['other_parameters'] = {**self._user_other, **mode}

    # --- ChunkedEquilibrateTask hooks --------------------------------------------------------------
    def _setup(self, state, min_steps):
        specs = self.specs
        self._mass_amu = total_mass_amu(state.psf.path)
        n_atoms = total_atoms(state.psf.path)
        # Measure per-leaflet lipid counts + protein cross-sections from the embedded frame so APL is
        # reported protein-corrected and per-leaflet (see membrane_leaflet_geometry).  A supplied
        # `lipids_per_leaflet` spec, if any, is only a fallback for when the geometry can't be measured.
        self._geom = None
        coor = getattr(state, 'pdb', None) or getattr(state, 'coor', None)
        try:
            if coor is None:
                raise ValueError('no coordinate file in state')
            self._geom = membrane_leaflet_geometry(state.psf.path, coor.path)
            logger.info(f'{self.taskname}: leaflet geometry -- lipids/leaflet '
                        f'{self._geom.n_lower}(lower)/{self._geom.n_upper}(upper); protein footprint '
                        f'{self._geom.a_prot_lower:.0f}/{self._geom.a_prot_upper:.0f} A^2 '
                        f'(APL reported protein-corrected, per leaflet)')
        except Exception as e:
            logger.warning(f'{self.taskname}: could not measure leaflet geometry ({e}); '
                           f'APL will use the lipids_per_leaflet spec if given')
        self._lipids_per_leaflet = int(specs.get('lipids_per_leaflet') or 0)
        self._mon_d = DensityConvergenceMonitor(self._params_for('density', min_steps))
        # Area is a soft, slow mode; a larger floor (area_min_steps) prevents premature convergence on
        # a momentary flat spot during the initial area relaxation (see the area_* schema defaults).
        self._min_steps = min_steps
        self._area_min = max(min_steps, int(specs.get('area_min_steps') or 0) or min_steps)
        self._mon_a = DensityConvergenceMonitor(self._params_for('area', self._area_min))
        self._n_consecutive = int(specs['n_consecutive'])
        self._all_t, self._all_d, self._all_a = [], [], []   # series for the two-panel plot
        self._phase2_start_step = None                       # set at the stage-1 -> stage-2 handoff

        # Two-stage protocol (default on): stage 1 settles density at *constant lateral area* so the
        # under-dense box loses its excess volume from z only; stage 2 then relaxes the lateral area at
        # constant (already-correct) density under tensionless NPgT.  A single tensionless piston asked
        # to do both at once co-compacts the area down with the shrinking box (metastable gel-like APL);
        # decoupling the two removes that failure mode.  See docs/design/membrane-equilibrate.md.
        self._user_other = dict(specs.get('other_parameters', {}))
        self._two_stage = bool(specs.get('two_stage', True))
        self._phase = 1 if self._two_stage else 2
        self._set_barostat_mode(self._phase)
        if self._two_stage:
            # stage-1 gates on DENSITY alone (area is pinned, so it is not a meaningful observable yet)
            self._conv = JointConvergence({'density': self._mon_d}, n_consecutive=self._n_consecutive)
            logger.info(f'{self.taskname}: {n_atoms} atoms; two-stage NPgT -- stage 1 settles density '
                        f'at constant lateral area (gate < {self._mon_d.precision_gate():.2e}), then '
                        f'stage 2 relaxes area tensionless (gate < {self._mon_a.precision_gate():.2e}); '
                        f'autocorrelation-corrected')
        else:
            self._conv = JointConvergence({'density': self._mon_d, 'area': self._mon_a},
                                          n_consecutive=self._n_consecutive)
            logger.info(f'{self.taskname}: {n_atoms} atoms; single-stage NPgT density+area convergence '
                        f'(density gate < {self._mon_d.precision_gate():.2e}, area gate < '
                        f'{self._mon_a.precision_gate():.2e}; autocorrelation-corrected)')

    def _ingest_chunk(self, xst, total_steps, this_chunk, n_chunk):
        ts, vols = xst_cell_volumes(xst)
        _, areas = xst_cell_areas(xst)
        last_area = None
        if ts.size == 0:
            logger.warning(f'{self.taskname}: chunk {n_chunk} .xst had no frames; '
                           f'continuing (raise xstfreq resolution if this persists)')
        else:
            dens = volume_to_density(vols, self._mass_amu)
            last_area = float(areas.mean())
            self._conv.add_samples('density', ts, dens)
            if self._phase == 2:                     # area is only a live observable once it can move
                self._conv.add_samples('area', ts, areas)
            self._all_t.extend(ts.tolist())
            self._all_d.extend(dens.tolist())
            self._all_a.extend(areas.tolist())
        jr = self._conv.check()
        self._rows.append((n_chunk, total_steps, this_chunk, jr, self._phase))
        rd, ra = jr.reports.get('density'), jr.reports.get('area')
        # in stage 1 the area monitor is absent; report the measured (pinned) area/APL for legibility,
        # with drift shown as n/a since the area is not yet relaxing
        area = ra.mean_density if ra else last_area
        apl = self._apl_mean(area)
        # display the *signed* drift (negative = shrinking, positive = growing) so the direction of
        # the trend is legible; the convergence gate still keys off the magnitude (rd/ra.drift).
        stage = f'[stage {self._phase}] ' if self._two_stage else ''
        logger.info(f'{self.taskname}: {stage}chunk {n_chunk} @ step {total_steps} -- '
                    f'rho={_fmt(rd.mean_density if rd else None)} g/cc '
                    f'(drift {_fmt(rd.signed_drift if rd else None)}), '
                    f'area={_fmt(area)} A^2 (drift {_fmt(ra.signed_drift if ra else None)}'
                    f'{"" if apl is None else f", APL~{_fmt(apl)}"}) -- {jr.reason}')

        # stage 1 -> stage 2 handoff: density has settled at constant area; hand off to tensionless
        # area relaxation.  Reset both monitors so stage-2 convergence is judged on stage-2 samples only,
        # and require area_min more steps of *actual* relaxation before the area may certify.
        if self._phase == 1 and jr.converged:
            self._phase = 2
            self._phase2_start_step = total_steps
            self._set_barostat_mode(2)
            self._mon_d.reset()
            self._mon_a.reset()
            self._mon_d.params.min_steps = total_steps + self._min_steps
            self._mon_a.params.min_steps = total_steps + self._area_min
            self._conv = JointConvergence({'density': self._mon_d, 'area': self._mon_a},
                                          n_consecutive=self._n_consecutive)
            logger.info(f'{self.taskname}: stage 1 complete at step {total_steps} (density settled at '
                        f'constant area, APL~{_fmt(apl)}); switching to tensionless area relaxation')
            return jr.blowup, False
        return jr.blowup, jr.converged

    def _converged_stop_reason(self, total_steps):
        return (f'CONVERGED at step {total_steps}: density + lateral area both stationary for '
                f'{self._rows[-1][3].passes} consecutive checks')

    def _ceiling_stop_reason(self, total_steps, max_steps):
        jr = self._rows[-1][3] if self._rows else None
        bits = []
        if jr:
            for name in ('density', 'area'):
                r = jr.reports.get(name)
                if r:
                    bits.append(f'{name} drift {_fmt(r.signed_drift)} (precision '
                                f'{"met" if r.precision_met else "UNMET"})')
        return (f'CEILING: reached max_steps ({max_steps}) without joint convergence; '
                f'{"; ".join(bits) or "n/a"} -- system may not have settled')

    def _blowup_stop_reason(self, total_steps):
        return 'ABORTED: box blowup (non-finite density or area)'

    def _blowup_error(self, total_steps):
        return (f'{self.taskname}: box blew up (non-finite density or area) at '
                f'step {total_steps}; see convergence report')

    def _write_outputs(self, stop_reason):
        self._write_report(stop_reason)
        self._write_plot(stop_reason)

    # --- report + plot -----------------------------------------------------------------------------
    def _apl_pair(self, area):
        """Return ``(apl_lower, apl_upper)`` in A^2/lipid -- protein-corrected and per-leaflet when the
        leaflet geometry was measured; otherwise the naive ``area/lipids_per_leaflet`` (spec) for both;
        ``(None, None)`` if neither is available."""
        if area is None:
            return None, None
        if self._geom is not None:
            return self._geom.apl_lower(area), self._geom.apl_upper(area)
        if self._lipids_per_leaflet:
            v = area / self._lipids_per_leaflet
            return v, v
        return None, None

    def _apl_mean(self, area):
        lo, up = self._apl_pair(area)
        if lo is not None and up is not None:
            return 0.5 * (lo + up)
        return lo if lo is not None else up

    def _write_report(self, stop_reason):
        """Write a per-chunk two-observable convergence report (``<basename>-membrane.dat``)."""
        fn = f'{self.basename}-membrane.dat'
        with open(fn, 'w') as f:
            f.write(f'# membrane_equilibrate convergence report -- {self.taskname}\n')
            if self._geom is not None:
                g = self._geom
                f.write(f'# total system mass: {self._mass_amu:.1f} amu; lipids/leaflet: '
                        f'{g.n_lower}(lower)/{g.n_upper}(upper); protein xy footprint: '
                        f'{g.a_prot_lower:.0f}/{g.a_prot_upper:.0f} A^2 (lower/upper)\n')
                f.write('# APL is protein-corrected per leaflet: (area - protein_footprint)/n_lipids\n')
            else:
                f.write(f'# total system mass: {self._mass_amu:.1f} amu; '
                        f'lipids/leaflet: {self._lipids_per_leaflet or "n/a"} (geometry unmeasured; '
                        f'APL = area/lipids_per_leaflet, NOT protein-corrected)\n')
            if self._two_stage:
                s2 = ('n/a (still in stage 1)' if self._phase2_start_step is None
                      else f'step {self._phase2_start_step}')
                f.write(f'# protocol: two-stage (stage 1 = const-area density settle, '
                        f'stage 2 = tensionless area relaxation; stage 2 began at {s2})\n')
            else:
                f.write('# protocol: single-stage tensionless NPgT (density+area jointly)\n')
            f.write(f'# stop: {stop_reason}\n')
            f.write('# stage  chunk  step  nsteps  rho[g/cc]  rho_drift  rho_SEM/m  area[A^2]  '
                    'APL_lo[A^2]  APL_up[A^2]  area_drift  area_SEM/m  passes  reason\n')
            for n, step, nsteps, jr, ph in self._rows:
                rd, ra = jr.reports.get('density'), jr.reports.get('area')
                area = ra.mean_density if ra else None
                apl_lo, apl_up = self._apl_pair(area)
                f.write(f'{ph:6d}  {n:6d}  {step:8d}  {nsteps:6d}  '
                        f'{_fmt(rd.mean_density if rd else None):>9}  '
                        f'{_fmt(rd.signed_drift if rd else None):>9}  {_fmt(rd.sem_over_mean if rd else None):>9}  '
                        f'{_fmt(area):>9}  {_fmt(apl_lo):>11}  {_fmt(apl_up):>11}  '
                        f'{_fmt(ra.signed_drift if ra else None):>10}  '
                        f'{_fmt(ra.sem_over_mean if ra else None):>10}  '
                        f'{jr.passes:6d}  {jr.reason}\n')
        self.register(f'{self.basename}-membrane', key='membrane_report', artifact_type=DataFileArtifact)
        logger.info(f'{self.taskname}: convergence report -> {fn}')

    def _write_plot(self, stop_reason=''):
        """Write a two-panel density-and-area-vs-time PNG (``<basename>-membrane.png``)."""
        if not self._all_t:
            return
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except Exception as e:  # pragma: no cover - plotting is best-effort
            logger.warning(f'{self.taskname}: could not import matplotlib for membrane plot ({e})')
            return
        import numpy as np
        t = np.asarray(self._all_t, dtype=float)
        d = np.asarray(self._all_d, dtype=float)
        a = np.asarray(self._all_a, dtype=float)
        t0, t1 = float(t.min()), float(t.max())
        burn = float(self.specs.get('burn_in', 0))
        window_frac = float(self.specs.get('window_frac', 0.5))
        conv_step = self._rows[-1][1] if (stop_reason.startswith('CONVERGED') and self._rows) else None

        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        panels = [(axes[0], d, 'density (g/cc)', '#3366aa'),
                  (axes[1], a, 'lateral area (A^2)', '#227744')]
        for ax, y, ylabel, color in panels:
            ax.plot(t, y, lw=0.8, color=color)
            if burn > 0 and t0 + burn < t1:
                ax.axvspan(t0, t0 + burn, color='0.85')
            post = t[t > t0 + burn]
            if 0.0 < window_frac < 1.0 and post.size:
                ws = float(post[int(post.size * (1.0 - window_frac))])
                ax.axvspan(ws, t1, color='#ddeecc', alpha=0.6)
            if conv_step is not None:
                ax.axvline(conv_step, ls='-', color='#118811', lw=1.2)
            # mark the stage 1 -> stage 2 (const-area -> tensionless) handoff
            if self._phase2_start_step is not None:
                ax.axvline(self._phase2_start_step, ls='--', color='#aa6633', lw=1.0)
            ax.set_ylabel(ylabel)
            ax.set_xlim(t0, t1)
        axes[0].set_title(f'{self.taskname}: density + lateral area vs. time'
                          + (f' — converged @ {conv_step:.0f}' if conv_step is not None else ''))
        # Secondary APL axis: mean protein-corrected APL is affine in area, so a linear twin axis with
        # limits mapped through _apl_mean is exact.  (Falls back to naive area/lipids if geometry is
        # unmeasured but a spec count was given.)
        lo0, lo1 = axes[1].get_ylim()
        apl0, apl1 = self._apl_mean(lo0), self._apl_mean(lo1)
        if apl0 is not None and apl1 is not None:
            ax2 = axes[1].twinx()
            ax2.set_ylim(apl0, apl1)
            ax2.set_ylabel('mean APL (A^2/lipid, protein-corr.)'
                           if self._geom is not None else 'APL (A^2/lipid)')
        axes[1].set_xlabel('timestep')
        fig.tight_layout()
        fn = f'{self.basename}-membrane.png'
        fig.savefig(fn, dpi=120)
        plt.close(fig)
        self.register(fn, key='membrane_plot', artifact_type=PNGImageFileArtifact, keep=True)
        logger.info(f'{self.taskname}: membrane plot -> {fn}')
