"""Unit tests for pestifer.util.ppewald: the production profile plus a replayed reciprocal term.

The engine is deliberately testable without NAMD: config synthesis, log parsing, step matching,
the MTS refusal and the per-component closure check all run on strings and arrays.  The
`needs_tools` tests at the end drive real NAMD end to end.
"""
import os
import textwrap

import numpy as np
import pytest

from pestifer.util.ppewald import (
    CLOSURE_TOLERANCE_BAR,
    DROP_DIRECTIVES,
    ProfilePass,
    ReplayResult,
    align_production,
    build_replay_config,
    check_against_tensor,
    closure_warnings,
    combine,
    full_elect_frequency,
    parse_profiles,
    replay,
    strip_for_replay,
    write_csv,
)

PROD = textwrap.dedent("""\
    structure               sys.psf
    coordinates             sys.pdb
    parameters              par_all36.prm
    paraTypeCharmm          on
    cutoff                  12.0
    PME                     on
    PMEGridSizeX            64
    # a comment
    langevin                on
    langevinTemp            310
    langevinPiston          on
    langevinPistonTarget    1.01325
    useFlexibleCell         yes
    temperature             310
    outputName              prod
    dcdfreq                 500
    DCDUnitCell             yes
    restartfreq             5000
    xstfreq                 500
    CUDASOAintegrate        on
    pressureProfile         on
    pressureProfileFreq     100
    pressureProfileSlabs    20
    colvars                 on
    colvarsConfig           restraints.in
    run 128000
    """).splitlines()


class TestStripForReplay:
    """What survives into a replay config, and what must not."""

    def test_force_field_lines_pass_through_verbatim(self):
        kept, _ = strip_for_replay(PROD)
        for line in ('structure               sys.psf',
                     'parameters              par_all36.prm',
                     'cutoff                  12.0',
                     'PME                     on',
                     'PMEGridSizeX            64'):
            assert line in kept

    def test_restraints_survive_because_they_enter_the_virial(self):
        # colvars forces contribute to the pressure, so dropping them would change the profile.
        kept, _ = strip_for_replay(PROD)
        assert any(k.startswith('colvars ') for k in kept)
        assert any('colvarsConfig' in k for k in kept)

    @pytest.mark.parametrize('directive', [
        'langevin', 'langevinPiston', 'useFlexibleCell',   # nothing is integrated
        'outputName', 'dcdfreq', 'restartfreq', 'xstfreq',  # must not clobber the run
        'CUDASOAintegrate',                                 # unsupported for this
        'pressureProfile', 'pressureProfileFreq', 'pressureProfileSlabs',
        'temperature',
    ])
    def test_production_directives_are_dropped(self, directive):
        kept, dropped = strip_for_replay(PROD)
        assert any(d.strip().split()[0] == directive for d in dropped), \
            f'{directive} should have been dropped'
        assert not any(k.strip().split()[:1] == [directive] for k in kept)

    def test_the_trailing_run_command_is_dropped(self):
        # A production config ends by running dynamics; only the replay loop may drive NAMD.
        _, dropped = strip_for_replay(PROD)
        assert 'run 128000' in [d.strip() for d in dropped]

    def test_comments_and_blank_lines_survive(self):
        kept, _ = strip_for_replay(PROD)
        assert '# a comment' in kept

    def test_matching_is_case_insensitive(self):
        _, dropped = strip_for_replay(['LANGEVINPISTON on', 'DCDfreq 100'])
        assert len(dropped) == 2

    def test_drop_list_is_lowercase(self):
        # matching lowercases the token, so an uppercase entry here would never fire
        assert all(d == d.lower() for d in DROP_DIRECTIVES)


def _cfg(**kw):
    args = dict(slabs=20, outname='out', ewald_grid=(10, 10, 10))
    args.update(kw)
    return build_replay_config(PROD, kw.pop('dcd', 'traj.dcd'), **args)


class TestBuildReplayConfig:
    def test_the_replay_is_the_reciprocal_term_only(self):
        cfg = _cfg(ewald_grid=(12, 13, 14))
        assert 'pressureProfileEwald    on' in cfg
        for axis, n in zip('XYZ', (12, 13, 14)):
            assert f'pressureProfileEwald{axis}   {n}' in cfg

    def test_pressureprofile_is_on_too_or_namd_writes_nothing(self):
        # measured 2026-11-05: with pressureProfileEwald on but pressureProfile absent, NAMD ran
        # every frame and wrote no PRESSUREPROFILE record at all
        assert 'pressureProfile         on' in _cfg()

    def test_numsteps_zero_and_outputname_present(self):
        # NAMD requires outputName even when nothing is integrated, and numsteps must be 0 or the
        # replay integrates instead of replaying.
        cfg = _cfg(outname='myout')
        assert 'numsteps                0' in cfg
        assert 'outputName              myout' in cfg

    def test_cwd_precedes_the_passthrough_lines(self):
        # NAMD chdirs to the config file's directory, so `cwd` must come before any relative path
        # the production config named, or those paths resolve against the wrong directory.
        cfg = _cfg(basedir='/run/dir')
        assert cfg.index('cwd /run/dir') < cfg.index('structure               sys.psf')

    def test_no_basedir_emits_no_cwd(self):
        assert 'cwd ' not in _cfg()

    def test_replay_loop_reads_the_named_trajectory(self):
        cfg = build_replay_config(PROD, '/abs/traj.dcd', slabs=20, outname='out',
                                  ewald_grid=(10, 10, 10))
        assert 'coorfile open dcd /abs/traj.dcd' in cfg
        assert 'run 0' in cfg
        assert 'coorfile close' in cfg

    def test_stride_skips_between_reads(self):
        cfg = _cfg(stride=5)
        assert 'coorfile skip' in cfg
        assert '$i < 5' in cfg

    def test_stride_one_emits_no_skip(self):
        assert 'coorfile skip' not in _cfg(stride=1)

    def test_slab_count_is_the_one_given(self):
        assert 'pressureProfileSlabs    17' in _cfg(slabs=17)


class TestFullElectFrequency:
    def test_default_is_one(self):
        assert full_elect_frequency(PROD) == 1

    def test_last_setting_wins_case_insensitively(self):
        assert full_elect_frequency(['fullElectFrequency 1', 'FULLELECTFREQUENCY 2 # mts']) == 2


def _log(records, thickness=2.0, frames=None):
    """A minimal replay log carrying the lines the parser looks for."""
    out = ['Info: PRESSURE PROFILE CALCULATIONS ACTIVE',
           f'Info:       SLAB THICKNESS: {thickness}']
    for rec in records:
        out.append('PRESSUREPROFILE: 0 ' + ' '.join(str(v) for v in rec))
    if frames is not None:
        out.append(f'TCL: PPREPLAY_FRAMES {frames}')
    return '\n'.join(out) + '\n'


def _prodlog(steps, profiles, pressavg=None, thickness=2.0):
    """A minimal production log: a step-numbered profile record, and optionally PRESSAVG, per step."""
    out = [f'Info:       SLAB THICKNESS: {thickness}']
    for i, st in enumerate(steps):
        out.append(f'PRESSUREPROFILE: {st} ' + ' '.join(str(v) for v in np.ravel(profiles[i])))
        if pressavg is not None:
            xx, yy, zz = pressavg[i]
            out.append(f'PRESSAVG: {st} {xx} 0 0 0 {yy} 0 0 0 {zz}')
    return '\n'.join(out) + '\n'


class TestParseProfiles:
    def test_shape_is_frames_by_slabs_by_three(self, tmp_path):
        p = tmp_path / 'r.log'
        p.write_text(_log([[1, 2, 3, 4, 5, 6], [7, 8, 9, 10, 11, 12]], frames=2))
        got = parse_profiles(str(p), 'real')
        assert got.profiles.shape == (2, 2, 3)
        assert got.nframes == 2 and got.nslabs == 2
        assert got.frames_read == 2
        assert got.slab_thickness == 2.0

    def test_records_are_kept_in_file_order(self, tmp_path):
        # every frame is evaluated at timestep 0, so the step column carries no frame identity and
        # order is the only thing aligning the replay with the frames it read
        p = tmp_path / 'r.log'
        p.write_text(_log([[1, 1, 1], [2, 2, 2], [3, 3, 3]]))
        got = parse_profiles(str(p), 'real')
        assert [row[0][0] for row in got.profiles] == [1.0, 2.0, 3.0]

    def test_no_records_raises_rather_than_returning_empty(self, tmp_path):
        p = tmp_path / 'r.log'
        p.write_text('Info: nothing here\nFATAL ERROR: something\n')
        with pytest.raises(RuntimeError, match='no PRESSUREPROFILE records'):
            parse_profiles(str(p), 'real')

    def test_ragged_records_raise(self, tmp_path):
        p = tmp_path / 'r.log'
        p.write_text(_log([[1, 2, 3], [1, 2, 3, 4, 5, 6]]))
        with pytest.raises(RuntimeError, match='ragged'):
            parse_profiles(str(p), 'real')

    def test_width_not_a_multiple_of_three_raises(self, tmp_path):
        p = tmp_path / 'r.log'
        p.write_text(_log([[1, 2, 3, 4]]))
        with pytest.raises(RuntimeError, match='not a multiple of 3'):
            parse_profiles(str(p), 'real')

    def test_missing_slab_thickness_is_tolerated(self, tmp_path):
        p = tmp_path / 'r.log'
        p.write_text('PRESSUREPROFILE: 0 1 2 3\n')
        assert parse_profiles(str(p), 'real').slab_thickness is None


def _pass(mode, arr, thickness=2.0):
    return ProfilePass(mode=mode, profiles=np.asarray(arr, dtype=float), slab_thickness=thickness)


class TestAlignProduction:
    """The production profile is matched to trajectory frames by MD step, never by position."""

    def test_records_are_matched_by_step(self, tmp_path):
        steps = [0, 100, 200, 300]
        prof = [np.full((2, 3), s) for s in steps]
        p = tmp_path / 'prod.log'
        p.write_text(_prodlog(steps, prof))
        # a DCD's first frame is firsttimestep + dcdfreq, so step 0 is never a frame
        real, tensor = align_production(str(p), np.array([100, 300]))
        assert [r[0, 0] for r in real.profiles] == [100.0, 300.0]
        assert tensor is None

    def test_a_frame_with_no_profile_record_raises(self, tmp_path):
        p = tmp_path / 'prod.log'
        p.write_text(_prodlog([0, 100], [np.zeros((2, 3))] * 2))
        with pytest.raises(RuntimeError, match='no production'):
            align_production(str(p), np.array([100, 150]))

    def test_pressavg_diagonal_is_captured(self, tmp_path):
        p = tmp_path / 'prod.log'
        p.write_text(_prodlog([100, 200], [np.zeros((1, 3))] * 2,
                              pressavg=[(1, 2, 3), (4, 5, 6)]))
        _, tensor = align_production(str(p), np.array([100, 200]))
        assert np.allclose(tensor, [[1, 2, 3], [4, 5, 6]])

    def test_a_log_without_profiles_says_what_it_needs(self, tmp_path):
        p = tmp_path / 'prod.log'
        p.write_text('ENERGY: 0 1 2 3\n')
        with pytest.raises(RuntimeError, match='pressureProfile on'):
            align_production(str(p), np.array([100]))


class TestCombine:
    def test_sum_is_elementwise(self):
        a = _pass('real', [[[1, 2, 3], [4, 5, 6]]])
        b = _pass('ewald', [[[10, 20, 30], [40, 50, 60]]])
        assert np.allclose(combine(a, b), [[[11, 22, 33], [44, 55, 66]]])

    def test_mismatched_frame_counts_raise(self):
        # the replay's records carry no step, so its frame count is the one chance to catch a
        # replay that did not read the frames the production records were matched to
        a = _pass('real', np.zeros((3, 2, 3)))
        b = _pass('ewald', np.zeros((2, 2, 3)))
        with pytest.raises(RuntimeError, match='frame counts differ'):
            combine(a, b)

    def test_mismatched_slab_counts_raise(self):
        a = _pass('real', np.zeros((2, 4, 3)))
        b = _pass('ewald', np.zeros((2, 5, 3)))
        with pytest.raises(RuntimeError, match='slab counts differ'):
            combine(a, b)


def _result(real, ewald, thickness=2.0, tensor=None):
    r, e = _pass('real', real, thickness), _pass('ewald', ewald, thickness)
    return ReplayResult(real=r, ewald=e, total=combine(r, e),
                        tensor=None if tensor is None else np.asarray(tensor, dtype=float))


class TestReplayResult:
    def test_slab_centers_are_offset_by_half_a_slab(self):
        res = _result(np.zeros((1, 3, 3)), np.zeros((1, 3, 3)), thickness=10.0)
        assert np.allclose(res.slab_centers(), [5.0, 15.0, 25.0])

    def test_slab_centers_fall_back_to_index_without_thickness(self):
        res = _result(np.zeros((1, 3, 3)), np.zeros((1, 3, 3)), thickness=None)
        assert np.allclose(res.slab_centers(), [0, 1, 2])

    def test_sem_is_zero_for_a_single_frame(self):
        res = _result(np.ones((1, 2, 3)), np.ones((1, 2, 3)))
        assert np.allclose(res.sem_profile(), 0.0)

    def test_mean_profile_averages_over_frames(self):
        real = np.array([[[0, 0, 0]], [[2, 2, 2]]], dtype=float)
        res = _result(real, np.zeros_like(real))
        assert np.allclose(res.mean_profile(), [[1, 1, 1]])

    def test_csv_has_one_row_per_slab(self, tmp_path):
        res = _result(np.ones((2, 4, 3)), np.zeros((2, 4, 3)))
        out = tmp_path / 'p.csv'
        write_csv(res, str(out))
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 5          # header + 4 slabs
        assert lines[0].startswith('z,real_xx')


class TestClosureCheck:
    """Per component, against the production run's own PRESSAVG -- not the trace.

    The shape that motivated it (2026-11-05): at fullElectFrequency 2 the sum missed PRESSAVG by
    -86 / -94 / -71 bar, an anisotropic error worth +14.5 mN/m, which is a surface tension's own
    currency.  A check on the trace alone weighs that error by its mean and sees nothing of its
    anisotropy; this one names each component.
    """

    def test_a_closed_profile_passes(self):
        real = np.full((4, 2, 3), 100.0)
        ewald = np.full((4, 2, 3), -98.0)
        res = _result(real, ewald, tensor=np.full((4, 3), 2.5))
        check = check_against_tensor(res)
        assert np.allclose(check['deviation'], -0.5)
        assert np.allclose(check['real_only_deviation'], 97.5)
        assert closure_warnings(check) == []

    def test_an_anisotropic_miss_is_flagged_even_with_a_small_trace_error(self):
        real = np.zeros((4, 1, 3))
        ewald = np.tile(np.array([[[+40.0, +40.0, -80.0]]]), (4, 1, 1))  # trace error 0
        res = _result(real, ewald, tensor=np.zeros((4, 3)))
        warnings = closure_warnings(check_against_tensor(res))
        assert len(warnings) == 1 and 'Ewald grid' in warnings[0]

    def test_the_tolerance_is_what_the_message_names(self):
        assert f'{CLOSURE_TOLERANCE_BAR:.0f}' in closure_warnings(
            {'deviation': np.array([0, 0, 1e3]), 'sem': np.zeros(3)})[0]

    def test_no_tensor_means_no_check_not_a_pass(self):
        res = _result(np.zeros((2, 1, 3)), np.zeros((2, 1, 3)))
        assert check_against_tensor(res) is None


class TestMultipleTimeSteppingIsRefused:
    def test_replay_refuses_before_running_namd(self, tmp_path, monkeypatch):
        from pestifer.util import ppewald
        cfg = tmp_path / 'prod.namd'
        cfg.write_text('\n'.join(PROD + ['fullElectFrequency 2']))
        ran = []
        monkeypatch.setattr(ppewald, 'run_namd', lambda *a, **k: ran.append(a))
        with pytest.raises(RuntimeError, match='counts it twice'):
            replay(str(cfg), str(tmp_path / 'x.dcd'), str(tmp_path / 'x.log'))
        assert ran == []


class TestSubcommandWarnings:
    """The subcommand reports the closure check and warns when there is none to report."""

    def _run(self, tmp_path, monkeypatch, caplog, res):
        import argparse
        import logging
        from pestifer.subcommands import pressure_profile_ewald as sub
        monkeypatch.setattr(sub, 'replay', lambda *a, **k: res)
        for name in ('prod.namd', 'prod.log', 'prod.dcd'):
            (tmp_path / name).write_text('')
        args = argparse.Namespace(
            namd_config=str(tmp_path / 'prod.namd'), log=str(tmp_path / 'prod.log'),
            dcd=str(tmp_path / 'prod.dcd'), workdir=str(tmp_path), slabs=None, ewald_grid=None,
            stride=1, temperature=300, namd='namd3', np=1, prefix='t', title=None,
            figsize=[6, 2], accept_mts=False)
        with caplog.at_level(logging.WARNING):
            sub.PressureProfileEwaldSubcommand.func(args)
        return [rec.getMessage() for rec in caplog.records if rec.levelno >= logging.WARNING]

    def test_a_profile_that_misses_the_tensor_is_flagged(self, tmp_path, monkeypatch, caplog):
        res = _result(np.full((3, 1, 3), -8200.0), np.full((3, 1, 3), 200.0),
                      tensor=np.full((3, 3), 150.0))
        warnings = self._run(tmp_path, monkeypatch, caplog, res)
        assert any('does not reproduce' in w for w in warnings), warnings

    def test_a_sound_profile_is_not_flagged(self, tmp_path, monkeypatch, caplog):
        res = _result(np.full((3, 1, 3), 100.0), np.full((3, 1, 3), 50.5),
                      tensor=np.full((3, 3), 150.0))
        assert self._run(tmp_path, monkeypatch, caplog, res) == []

    def test_no_tensor_is_said_out_loud(self, tmp_path, monkeypatch, caplog):
        res = _result(np.full((3, 1, 3), 100.0), np.full((3, 1, 3), 50.5))
        warnings = self._run(tmp_path, monkeypatch, caplog, res)
        assert any('NOT checked' in w for w in warnings), warnings


@pytest.mark.needs_tools
class TestReplayEndToEnd:
    """Drive real NAMD: a production run with a profile and a pressure tensor, then the replay.

    The trajectory is generated here rather than committed: ``*.dcd`` is ignored repo-wide (they
    are build output), and this test needs NAMD anyway, so making its own input costs a couple of
    seconds and keeps the fixture to text files.
    """

    DATA = os.path.join(os.path.dirname(__file__), 'test_ppewald')

    def _make_run(self, tmp_path):
        """Run the fixture water box for 200 steps with profile, tensor and DCD every 20."""
        from pestifer.util.ppewald import run_namd
        cfg = tmp_path / 'gen.namd'
        body = [ln for ln in open(os.path.join(self.DATA, 'prod.namd')).read().splitlines()
                if not ln.strip().lower().startswith(('run', 'outputname'))]
        cfg.write_text('\n'.join([f'cwd {self.DATA}', *body,
                                  f'outputName {tmp_path / "gen"}',
                                  'outputPressure 20', 'outputEnergies 20', 'run 200', '']))
        run_namd(str(cfg), str(tmp_path / 'gen.log'), self.DATA, nprocs=2)
        dcd = tmp_path / 'gen.dcd'
        assert dcd.exists(), f'NAMD wrote no trajectory; see {tmp_path / "gen.log"}'
        return cfg, tmp_path / 'gen.log', dcd

    def test_the_complete_profile_closes_on_the_pressure_tensor(self, tmp_path):
        cfg, log, dcd = self._make_run(tmp_path)
        res = replay(str(cfg), str(dcd), str(log), workdir=str(tmp_path), nprocs=2)
        assert res.nframes == 10 and res.nslabs == 10
        assert list(res.steps) == list(range(20, 201, 20))
        check = check_against_tensor(res)
        assert check is not None
        # the real-space half misses by hundreds of bar on every component; the sum closes
        assert np.all(np.abs(check['real_only_deviation']) > 100)
        assert closure_warnings(check) == [], check

    def test_outputs_are_written(self, tmp_path):
        from pestifer.util.ppewald import plot
        cfg, log, dcd = self._make_run(tmp_path)
        res = replay(str(cfg), str(dcd), str(log), workdir=str(tmp_path), nprocs=2)
        csv, png = tmp_path / 'p.csv', tmp_path / 'p.png'
        write_csv(res, str(csv))
        plot(res, str(png))
        assert len(csv.read_text().strip().splitlines()) == 11   # header + 10 slabs
        assert png.stat().st_size > 1000
