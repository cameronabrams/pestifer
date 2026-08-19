# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for NAMD restart-configuration generation.

``make-namd-restart`` reads a finished (or crashed) NAMD run and writes a config that resumes it.
It is the clearest example of pestifer operating on a run it did not manage -- the user brings
their own log, config and checkpoint files -- which makes it the piece of the NAMD layer most
obviously public API, and it was the least-tested module in that group (13%).

Nothing here runs NAMD.  ``NAMDConfig`` is a text parser and the restart logic is a decision about
which timestep and which checkpoint files to point at; both are exercised with synthesized logs
and configs.  The log fixtures are minimal but real in form: the parser dispatches on
``line.startswith(key)``, so a restart line carries no ``Info:`` prefix, which is the sort of
detail a hand-written fixture gets wrong silently.
"""

import argparse
import logging
import os
import tempfile
import unittest

from pestifer.util.namdrestart import NAMDConfig, make_namd_restart_subcommand


CONFIG = """\
# a NAMD configuration
set temperature 310
set outputbase run00
structure my.psf
coordinates my.pdb
bincoordinates old.coor
binvelocities old.vel
extendedsystem old.xsc
outputname $outputbase
firsttimestep 0
temperature $temperature
run 1000
"""


def _log(*, output='run00', running_for=1000, restarts=(100, 200), finished=True):
    """A NAMD log carrying exactly the fields the restart logic reads."""
    lines = [f'Info: OUTPUT FILENAME        {output}',
             f'TCL: Running for {running_for} steps']
    lines += [f'WRITING COORDINATES TO RESTART FILE AT STEP {t}' for t in restarts]
    if finished:
        lines.append('WallClock: 40.19  CPUTime: 38.63  Memory: 15.17 MB')
    return '\n'.join(lines) + '\n'


class _Sandbox(unittest.TestCase):
    """Each test runs in its own directory: the subcommand reads and writes relative paths."""

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(os.chdir, self._cwd)

    def write(self, name, text):
        with open(name, 'w') as f:
            f.write(text)
        return name

    def checkpoints(self, base):
        """The .coor/.vel/.xsc a restart points at; absent unless a test says otherwise."""
        for ext in ('coor', 'vel', 'xsc'):
            self.write(f'{base}.{ext}', 'binary-ish')

    def args(self, **kw):
        d = dict(namd_log='run00.log', config='run00.namd', new_base='run01', run=0, slurm=None)
        d.update(kw)
        return argparse.Namespace(**d)


class TestNAMDConfigParsing(_Sandbox):

    def cfg(self, text=CONFIG):
        return NAMDConfig(self.write('c.namd', text))

    def test_a_missing_file_is_reported_not_guessed(self):
        with self.assertRaises(FileNotFoundError):
            NAMDConfig('nope.namd')

    def test_line_kinds_are_classified(self):
        kinds = {l['ltype'] for l in self.cfg().lines}
        self.assertEqual(kinds, {'comment', 'varassign', 'command', 'blank'})

    def test_variable_assignments_are_collected(self):
        self.assertEqual(self.cfg().varsdefined,
                         {'temperature': '310', 'outputbase': 'run00'})

    def test_only_variables_actually_referenced_are_cited(self):
        cited = self.cfg().varscited
        self.assertCountEqual(cited, ['outputbase', 'temperature'])

    def test_a_banner_comment_is_not_mistaken_for_a_comment(self):
        c = self.cfg('###### banner ######\n# real comment\nrun 10\n')
        self.assertEqual([l['ltype'] for l in c.lines if l['ltype'] != 'blank'],
                         ['banner', 'comment', 'command'])

    def test_a_variable_set_twice_is_flagged(self):
        with self.assertLogs('pestifer.util.namdrestart', level='WARNING') as cm:
            self.cfg('set a 1\nset a 2\nrun 5\n')
        self.assertIn('set twice', ''.join(cm.output))

    def test_a_referenced_but_undefined_variable_is_flagged(self):
        # a silent miss here would produce a restart config NAMD rejects
        with self.assertLogs('pestifer.util.namdrestart', level='WARNING') as cm:
            self.cfg('outputname $nosuchvar\nrun 5\n')
        self.assertIn('nosuchvar', ''.join(cm.output))

    def test_braced_variable_references_are_understood(self):
        c = self.cfg('set base run00\noutputname ${base}\nrun 5\n')
        self.assertIn('base', c.varscited)


class TestNAMDConfigEditing(_Sandbox):

    def cfg(self, text=CONFIG):
        return NAMDConfig(self.write('c.namd', text))

    @staticmethod
    def _command(c, name):
        return next(l for l in c.lines
                    if l['ltype'] == 'command' and l['commandname'] == name)

    def test_replacing_a_literal_argument(self):
        c = self.cfg()
        c.replace_command('run', ['500'])
        self.assertEqual(self._command(c, 'run')['commandargs'], ['500'])

    def test_a_variable_argument_is_kept_and_its_value_rewritten(self):
        """``outputname $outputbase`` must stay a reference; the *variable* changes.

        Overwriting the reference with a literal would silently orphan every other use of that
        variable in the file.
        """
        c = self.cfg()
        c.replace_command('outputname', ['run01'])
        self.assertEqual(self._command(c, 'outputname')['commandargs'], ['$outputbase'])
        self.assertEqual(c.varsdefined['outputbase'], 'run01')

    def test_command_lookup_is_case_insensitive(self):
        c = self.cfg()
        c.replace_command('RUN', ['7'])
        self.assertEqual(self._command(c, 'run')['commandargs'], ['7'])

    def test_back_resolving_an_unknown_variable_is_reported(self):
        c = self.cfg()
        with self.assertLogs('pestifer.util.namdrestart', level='WARNING') as cm:
            c.var_backresolve('nosuchvar', '1')
        self.assertIn('Cannot back-resolve', ''.join(cm.output))

    def test_written_config_round_trips_through_the_parser(self):
        c = self.cfg()
        c.replace_command('run', ['500'])
        c.write('out.namd')
        again = NAMDConfig('out.namd')
        self.assertEqual(self._command(again, 'run')['commandargs'], ['500'])
        self.assertEqual(again.varsdefined['temperature'], '310')


class TestMakeRestart(_Sandbox):
    """The decision the subcommand exists to make: resume from where, for how long, using which
    checkpoint files."""

    def setup_run(self, **logkw):
        self.write('run00.namd', CONFIG)
        self.write('run00.log', _log(**logkw))

    def test_a_crashed_run_resumes_from_the_restart_checkpoint(self):
        self.setup_run(finished=False)
        self.checkpoints('run00.restart')
        make_namd_restart_subcommand(self.args(run=0))
        out = NAMDConfig('run01.namd')
        cmd = {l['commandname']: l['commandargs'] for l in out.lines if l['ltype'] == 'command'}
        self.assertEqual(cmd['bincoordinates'], ['run00.restart.coor'])
        self.assertEqual(cmd['firsttimestep'], ['200'], 'resume from the last restart written')
        self.assertEqual(cmd['run'], ['800'], '1000 requested - 200 done = 800 remaining')

    def test_a_finished_run_extends_using_the_plain_checkpoint(self):
        # a successful run has final .coor/.vel/.xsc, not .restart.* -- and the new step count is
        # what the user asked for, not a remainder
        self.setup_run(finished=True)
        self.checkpoints('run00')
        make_namd_restart_subcommand(self.args(run=250))
        out = NAMDConfig('run01.namd')
        cmd = {l['commandname']: l['commandargs'] for l in out.lines if l['ltype'] == 'command'}
        self.assertEqual(cmd['bincoordinates'], ['run00.coor'])
        self.assertEqual(cmd['run'], ['250'])

    def test_extending_a_finished_run_by_zero_steps_does_nothing(self):
        self.setup_run(finished=True)
        self.checkpoints('run00')
        with self.assertLogs('pestifer.util.namdrestart', level='WARNING') as cm:
            make_namd_restart_subcommand(self.args(run=0))
        self.assertIn('did not request any new time steps', ''.join(cm.output))
        self.assertFalse(os.path.exists('run01.namd'), 'no config should be written')

    def test_a_log_with_no_restart_checkpoint_is_refused(self):
        self.setup_run(restarts=(), finished=False)
        self.checkpoints('run00.restart')
        with self.assertRaises(ValueError) as cm:
            make_namd_restart_subcommand(self.args())
        self.assertIn('restart time step', str(cm.exception))

    def test_a_log_with_no_step_target_is_refused(self):
        self.write('run00.namd', CONFIG)
        self.write('run00.log',
                   'Info: OUTPUT FILENAME        run00\n'
                   'WRITING COORDINATES TO RESTART FILE AT STEP 100\n')
        self.checkpoints('run00.restart')
        with self.assertRaises(ValueError) as cm:
            make_namd_restart_subcommand(self.args())
        self.assertIn('Running for', str(cm.exception))

    def test_a_missing_checkpoint_file_names_itself(self):
        """Refusing here is the point: a config pointing at an absent .coor fails inside NAMD,
        minutes later, with a worse message."""
        self.setup_run(finished=False)
        # deliberately provide only two of the three
        self.write('run00.restart.coor', 'x')
        self.write('run00.restart.vel', 'x')
        with self.assertRaises(FileNotFoundError) as cm:
            make_namd_restart_subcommand(self.args())
        self.assertIn('run00.restart.xsc', str(cm.exception))

    def test_a_log_without_an_output_filename_is_reported(self):
        self.write('run00.namd', CONFIG)
        self.write('run00.log', 'TCL: Running for 1000 steps\n'
                                'WRITING COORDINATES TO RESTART FILE AT STEP 100\n')
        with self.assertLogs('pestifer.util.namdrestart', level='ERROR') as cm:
            with self.assertRaises(FileNotFoundError):
                make_namd_restart_subcommand(self.args())
        self.assertIn('No output filename', ''.join(cm.output))


class TestVariableBasedConfigs(_Sandbox):
    """Regression: a config using the standard ``set outputbase`` idiom silently kept the old
    basename, so the restart would have overwritten the very files it was resuming from.

    ``replace_command`` passes the argument through verbatim -- ``$outputbase`` -- while
    ``varsdefined`` is keyed on bare names, so the back-resolve never matched and the replacement
    was a no-op that only logged.
    """

    VAR_CONFIG = ('set outputbase run00\n'
                  'structure my.psf\n'
                  'bincoordinates old.coor\n'
                  'binvelocities old.vel\n'
                  'extendedsystem old.xsc\n'
                  'outputname $outputbase\n'
                  'firsttimestep 0\n'
                  'run 1000\n')

    def setup_run(self):
        self.write('run00.namd', self.VAR_CONFIG)
        self.write('run00.log', _log(finished=False))
        self.checkpoints('run00.restart')

    def test_the_new_basename_reaches_the_variable(self):
        self.setup_run()
        make_namd_restart_subcommand(self.args())
        out = NAMDConfig('run01.namd')
        self.assertEqual(out.varsdefined['outputbase'], 'run01')

    def test_the_restart_does_not_write_over_the_run_it_resumes(self):
        """The consequence that makes this worth a test: output going to run00 would clobber the
        checkpoint files the restart is reading."""
        self.setup_run()
        make_namd_restart_subcommand(self.args())
        text = open('run01.namd').read()
        self.assertNotIn('set outputbase run00', text)

    def test_a_braced_reference_is_handled_too(self):
        self.write('run00.namd', self.VAR_CONFIG.replace('$outputbase', '${outputbase}'))
        self.write('run00.log', _log(finished=False))
        self.checkpoints('run00.restart')
        make_namd_restart_subcommand(self.args())
        self.assertEqual(NAMDConfig('run01.namd').varsdefined['outputbase'], 'run01')


class TestSlurmScriptRotation(_Sandbox):
    """With ``--slurm`` the batch script is rewritten in place, so the previous one is preserved
    under a numbered name rather than overwritten."""

    def setup_run(self):
        self.write('run00.namd', CONFIG)
        self.write('run00.log', _log(finished=False))
        self.checkpoints('run00.restart')
        self.write('job.sh', '#!/bin/bash\nBASENAME=run00\nnamd3 $BASENAME.namd\n')

    def test_the_previous_script_is_kept_and_the_basename_updated(self):
        self.setup_run()
        make_namd_restart_subcommand(self.args(slurm='job.sh'))
        self.assertTrue(os.path.exists('%job.sh%-1'), 'the original script was not preserved')
        self.assertIn('run01', open('job.sh').read())

    def test_successive_restarts_do_not_clobber_the_backup(self):
        self.setup_run()
        make_namd_restart_subcommand(self.args(slurm='job.sh'))
        make_namd_restart_subcommand(self.args(new_base='run02', slurm='job.sh'))
        self.assertTrue(os.path.exists('%job.sh%-1'))
        self.assertTrue(os.path.exists('%job.sh%-2'))


if __name__ == '__main__':
    unittest.main()
