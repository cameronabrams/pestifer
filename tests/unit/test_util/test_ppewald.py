"""Unit tests for pestifer.util.ppewald (the NAMD-free half of the two-pass replay).

The engine is deliberately testable without NAMD: config synthesis, log parsing, the guard against
combining mismatched passes, and the validation arithmetic all run on strings and arrays.  One
`needs_tools` test at the end drives real NAMD end to end.
"""
import os
import textwrap

import numpy as np
import pytest

from pestifer.util.ppewald import (
    DROP_DIRECTIVES,
    ProfilePass,
    ReplayResult,
    build_replay_config,
    combine,
    parse_profiles,
    strip_for_replay,
    validate_against_total_pressure,
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


class TestBuildReplayConfig:
    def test_real_pass_has_no_ewald_block(self):
        cfg = build_replay_config(PROD, 'traj.dcd', 'real', slabs=20, outname='out')
        assert 'pressureProfileEwald' not in cfg
        assert 'pressureProfile         on' in cfg

    def test_ewald_pass_turns_ewald_on(self):
        cfg = build_replay_config(PROD, 'traj.dcd', 'ewald', slabs=20, outname='out',
                                  ewald_grid=(12, 13, 14))
        assert 'pressureProfileEwald    on' in cfg
        for axis, n in zip('XYZ', (12, 13, 14)):
            assert f'pressureProfileEwald{axis}   {n}' in cfg

    def test_numsteps_zero_and_outputname_present(self):
        # NAMD requires outputName even when nothing is integrated, and numsteps must be 0 or the
        # replay integrates instead of replaying.
        cfg = build_replay_config(PROD, 'traj.dcd', 'real', slabs=20, outname='myout')
        assert 'numsteps                0' in cfg
        assert 'outputName              myout' in cfg

    def test_cwd_precedes_the_passthrough_lines(self):
        # NAMD chdirs to the config file's directory, so `cwd` must come before any relative path
        # the production config named, or those paths resolve against the wrong directory.
        cfg = build_replay_config(PROD, 'traj.dcd', 'real', slabs=20, outname='out',
                                  basedir='/run/dir')
        assert cfg.index('cwd /run/dir') < cfg.index('structure               sys.psf')

    def test_no_basedir_emits_no_cwd(self):
        cfg = build_replay_config(PROD, 'traj.dcd', 'real', slabs=20, outname='out')
        assert 'cwd ' not in cfg

    def test_replay_loop_reads_the_named_trajectory(self):
        cfg = build_replay_config(PROD, '/abs/traj.dcd', 'real', slabs=20, outname='out')
        assert 'coorfile open dcd /abs/traj.dcd' in cfg
        assert 'run 0' in cfg
        assert 'coorfile close' in cfg

    def test_stride_skips_between_reads(self):
        cfg = build_replay_config(PROD, 't.dcd', 'real', slabs=20, outname='out', stride=5)
        assert 'coorfile skip' in cfg
        assert '$i < 5' in cfg

    def test_stride_one_emits_no_skip(self):
        cfg = build_replay_config(PROD, 't.dcd', 'real', slabs=20, outname='out', stride=1)
        assert 'coorfile skip' not in cfg

    def test_both_passes_get_the_same_slab_count(self):
        # summing profiles gridded differently would be meaningless
        a = build_replay_config(PROD, 't.dcd', 'real', slabs=17, outname='o')
        b = build_replay_config(PROD, 't.dcd', 'ewald', slabs=17, outname='o')
        assert 'pressureProfileSlabs    17' in a
        assert 'pressureProfileSlabs    17' in b

    def test_rejects_an_unknown_mode(self):
        with pytest.raises(ValueError, match="'real' or 'ewald'"):
            build_replay_config(PROD, 't.dcd', 'reciprocal', slabs=20, outname='o')


def _log(records, thickness=2.0, pressures=None, frames=None):
    """A minimal NAMD log carrying the lines the parser looks for."""
    out = ['Info: PRESSURE PROFILE CALCULATIONS ACTIVE',
           f'Info:       SLAB THICKNESS: {thickness}']
    if pressures is not None:
        out.append('ETITLE:      TS           BOND            TEMP       PRESSURE')
        for i, (rec, p) in enumerate(zip(records, pressures)):
            out.append(f'ENERGY:       0           0.0           300.0       {p}')
            out.append('PRESSUREPROFILE: 0 ' + ' '.join(str(v) for v in rec))
    else:
        for rec in records:
            out.append('PRESSUREPROFILE: 0 ' + ' '.join(str(v) for v in rec))
    if frames is not None:
        out.append(f'TCL: PPREPLAY_FRAMES {frames}')
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
        # order is the only thing aligning the two passes
        p = tmp_path / 'r.log'
        p.write_text(_log([[1, 1, 1], [2, 2, 2], [3, 3, 3]]))
        got = parse_profiles(str(p), 'real')
        assert [row[0][0] for row in got.profiles] == [1.0, 2.0, 3.0]

    def test_energy_pressure_column_is_captured(self, tmp_path):
        p = tmp_path / 'r.log'
        p.write_text(_log([[1, 1, 1], [2, 2, 2]], pressures=[-100.5, 250.25]))
        got = parse_profiles(str(p), 'real')
        assert np.allclose(got.total_pressure, [-100.5, 250.25])

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


def _pass(mode, arr, thickness=2.0, pressures=None):
    return ProfilePass(mode=mode, profiles=np.asarray(arr, dtype=float),
                       slab_thickness=thickness,
                       total_pressure=None if pressures is None else np.asarray(pressures))


class TestCombine:
    def test_sum_is_elementwise(self):
        a = _pass('real', [[[1, 2, 3], [4, 5, 6]]])
        b = _pass('ewald', [[[10, 20, 30], [40, 50, 60]]])
        assert np.allclose(combine(a, b), [[[11, 22, 33], [44, 55, 66]]])

    def test_mismatched_frame_counts_raise(self):
        # positional alignment is the only alignment there is, so this is the one chance to catch
        # two passes that did not read the same frames
        a = _pass('real', np.zeros((3, 2, 3)))
        b = _pass('ewald', np.zeros((2, 2, 3)))
        with pytest.raises(RuntimeError, match='frame counts differ'):
            combine(a, b)

    def test_mismatched_slab_counts_raise(self):
        a = _pass('real', np.zeros((2, 4, 3)))
        b = _pass('ewald', np.zeros((2, 5, 3)))
        with pytest.raises(RuntimeError, match='slab counts differ'):
            combine(a, b)


class TestReplayResult:
    def _result(self, real, ewald, thickness=2.0, pressures=None):
        r, e = _pass('real', real, thickness, pressures), _pass('ewald', ewald, thickness)
        return ReplayResult(real=r, ewald=e, total=combine(r, e))

    def test_slab_centers_are_offset_by_half_a_slab(self):
        res = self._result(np.zeros((1, 3, 3)), np.zeros((1, 3, 3)), thickness=10.0)
        assert np.allclose(res.slab_centers(), [5.0, 15.0, 25.0])

    def test_slab_centers_fall_back_to_index_without_thickness(self):
        res = self._result(np.zeros((1, 3, 3)), np.zeros((1, 3, 3)), thickness=None)
        assert np.allclose(res.slab_centers(), [0, 1, 2])

    def test_sem_is_zero_for_a_single_frame(self):
        res = self._result(np.ones((1, 2, 3)), np.ones((1, 2, 3)))
        assert np.allclose(res.sem_profile(), 0.0)

    def test_mean_profile_averages_over_frames(self):
        real = np.array([[[0, 0, 0]], [[2, 2, 2]]], dtype=float)
        res = self._result(real, np.zeros_like(real))
        assert np.allclose(res.mean_profile(), [[1, 1, 1]])

    def test_validation_prefers_the_reconstruction_when_the_sum_is_right(self):
        # one slab, so the slab average is the value itself; NAMD's pressure equals real+ewald
        real = np.array([[[100.0, 100.0, 100.0]]])
        ewald = np.array([[[50.0, 50.0, 50.0]]])
        res = self._result(real, ewald, pressures=[150.0])
        check = validate_against_total_pressure(res)
        assert check['reconstructed_deviation'] == pytest.approx(0.0)
        assert check['real_only_deviation'] == pytest.approx(50.0)

    def test_validation_returns_none_without_energy_lines(self):
        res = self._result(np.zeros((1, 1, 3)), np.zeros((1, 1, 3)))
        assert validate_against_total_pressure(res) is None

    def test_csv_has_one_row_per_slab(self, tmp_path):
        res = self._result(np.ones((2, 4, 3)), np.zeros((2, 4, 3)))
        out = tmp_path / 'p.csv'
        write_csv(res, str(out))
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 5          # header + 4 slabs
        assert lines[0].startswith('z,real_xx')


@pytest.mark.needs_tools
class TestReplayEndToEnd:
    """Drive real NAMD over a real trajectory and check the reconstruction against NAMD itself.

    The trajectory is generated here rather than committed: ``*.dcd`` is ignored repo-wide (they
    are build output), and this test needs NAMD anyway, so making its own input costs a couple of
    seconds and keeps the fixture to text files.
    """

    DATA = os.path.join(os.path.dirname(__file__), 'test_ppewald')

    def _make_trajectory(self, tmp_path):
        """Run 60 steps of the fixture water box, writing a 3-frame DCD into tmp_path."""
        from pestifer.util.ppewald import run_namd
        cfg = tmp_path / 'gen.namd'
        body = [ln for ln in open(os.path.join(self.DATA, 'prod.namd')).read().splitlines()
                if not ln.strip().lower().startswith(('run', 'outputname', 'pressureprofile'))]
        cfg.write_text('\n'.join([f'cwd {self.DATA}', *body,
                                  f'outputName {tmp_path / "gen"}',
                                  'run 60', '']))
        run_namd(str(cfg), str(tmp_path / 'gen.log'), self.DATA, nprocs=2)
        dcd = tmp_path / 'gen.dcd'
        assert dcd.exists(), f'NAMD wrote no trajectory; see {tmp_path / "gen.log"}'
        return dcd

    def test_reconstruction_beats_real_space_alone(self, tmp_path):
        from pestifer.util.ppewald import replay
        dcd = self._make_trajectory(tmp_path)
        res = replay(os.path.join(self.DATA, 'prod.namd'), str(dcd),
                     workdir=str(tmp_path), slabs=10, nprocs=2)
        assert res.nframes > 0
        assert res.nslabs == 10
        # the two passes must be genuinely different quantities, not the same one twice
        assert not np.allclose(res.real.profiles, res.ewald.profiles)
        check = validate_against_total_pressure(res)
        assert check is not None
        # the whole point: the sum agrees with NAMD's own pressure and the real-space half does not
        assert check['reconstructed_deviation'] < check['real_only_deviation'] / 10

    def test_outputs_are_written(self, tmp_path):
        from pestifer.util.ppewald import plot, replay
        dcd = self._make_trajectory(tmp_path)
        res = replay(os.path.join(self.DATA, 'prod.namd'), str(dcd),
                     workdir=str(tmp_path), slabs=10, nprocs=2)
        csv, png = tmp_path / 'p.csv', tmp_path / 'p.png'
        write_csv(res, str(csv))
        plot(res, str(png))
        assert len(csv.read_text().strip().splitlines()) == 11   # header + 10 slabs
        assert png.stat().st_size > 1000
