import unittest
import os
from pathlib import Path
from pestifer.core.controller import Controller
from pestifer.core.config import Config 
from pestifer.tasks.continuation import ContinuationTask
from pestifer.tasks.terminate import TerminateTask
from pestifer.molecule.molecule import Molecule
from pestifer.core.artifacts import StateArtifacts, YAMLFileArtifact
from pestifer.core.errors import PestiferBuildError
class TestContinuationTask(unittest.TestCase):

    def setUp(self):
        self.controller=Controller().configure(Config().configure_new(), terminate=False)    
        input_dir = Path('../fixtures/continuation_inputs')
        # copy inputs to working directory
        psf = input_dir / 'my_6pti.psf'
        pdb = input_dir / 'my_6pti.pdb'
        xsc = input_dir / 'my_6pti.xsc'
        coor = input_dir / 'my_6pti.coor'
        vel = input_dir / 'my_6pti.vel'
        # copy to cwd
        dest_psf = Path('my_6pti.psf')
        dest_pdb = Path('my_6pti.pdb')
        dest_xsc = Path('my_6pti.xsc')
        dest_coor = Path('my_6pti.coor')
        dest_vel = Path('my_6pti.vel')
        dest_chainmap = Path('chainmap.yaml')
        if dest_chainmap.exists():
            dest_chainmap.unlink()
        dest_tarball = Path('my_system.tar.gz')
        if dest_tarball.exists():
            dest_tarball.unlink()
        if dest_psf.exists():
            dest_psf.unlink()
        if dest_pdb.exists():
            dest_pdb.unlink()
        if dest_xsc.exists():
            dest_xsc.unlink()
        if dest_coor.exists():
            dest_coor.unlink()
        if dest_vel.exists():
            dest_vel.unlink()
        os.symlink(psf.resolve(), dest_psf)
        os.symlink(pdb.resolve(), dest_pdb)
        os.symlink(xsc.resolve(), dest_xsc)
        os.symlink(coor.resolve(), dest_coor)
        os.symlink(vel.resolve(), dest_vel)

    def test_continuation_task(self):
        task_list = [{'continuation': dict(pdb='my_6pti.pdb', psf='my_6pti.psf', xsc='my_6pti.xsc', coor='my_6pti.coor',vel='my_6pti.vel')}]
        self.controller.reconfigure_tasks(task_list)
        self.assertEqual(len(self.controller.tasks), 1)
        self.assertIsInstance(self.controller.tasks[0], ContinuationTask)
        self.controller.do_tasks()
        task = self.controller.tasks[0]
        state: StateArtifacts = task.get_current_artifact('state')
        psf = state.psf.name
        self.assertEqual(psf, 'my_6pti.psf')
        # continuation registers the STATE fileset only; the full Molecule ingest is deferred
        # (MD-family tasks never need it), so no base_molecule is produced here.
        base_molecule = task.get_current_artifact_data('base_molecule')
        self.assertIsNone(base_molecule)

    def test_continuation_hex_serials_no_choke(self):
        # A ready-to-run system large enough to carry hexadecimal PDB atom serials must not make
        # the continuation task choke: because continuation no longer ingests a full Molecule,
        # the hex-serial PDB is never parsed on this path.  Build a hex-serial variant of the
        # fixture and confirm the task completes cleanly with no molecule ingested.
        hex_pdb = Path('my_6pti_hex.pdb')
        has_hex = False
        with open('my_6pti.pdb') as src, open(hex_pdb, 'w') as dst:
            n = 0
            for line in src:
                if line.startswith(('ATOM', 'HETATM')):
                    n += 1
                    s = n + 90000                     # push serials across the 99999 -> hex line
                    enc = format(s, 'X') if s > 99999 else str(s)
                    if any(ch in 'ABCDEF' for ch in enc):
                        has_hex = True
                    line = f'{line[:6]}{enc:>5}{line[11:]}'
                dst.write(line)
        self.assertTrue(has_hex, 'test setup failed to produce hexadecimal serials')
        task_list = [{'continuation': dict(pdb=str(hex_pdb), psf='my_6pti.psf',
                                            coor='my_6pti.coor', xsc='my_6pti.xsc', vel='my_6pti.vel')}]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        task = self.controller.tasks[0]
        self.assertEqual(task.result, 0)
        self.assertIsNone(task.get_current_artifact_data('base_molecule'))
        # And when the molecule IS genuinely needed later (the lazy path), parsing the hex-serial
        # PDB must succeed rather than choke -- this is what would happen for a hex system that
        # also requests, e.g., a terminate chainmap.
        bm = task.ensure_base_molecule()
        self.assertIsInstance(bm, Molecule)
        self.assertGreater(bm.num_atoms(), 0)
        hex_pdb.unlink()

    def test_verify_parameters_rejects_foreign_atom_type(self):
        # Give one atom a bogus atom type absent from the build's CHARMM release; the default-on
        # force-field-consistency preflight must hard-error (naming the term) rather than let the
        # system through to fail cryptically in a later NAMD run.  Only the PSF text is changed --
        # the coordinate files still match atom-for-atom, so the check reads the PSF alone.
        lines = Path('my_6pti.psf').read_text().splitlines()
        seen_natom = False
        out = []
        patched = False
        for line in lines:
            if '!NATOM' in line:
                seen_natom = True
                out.append(line)
                continue
            toks = line.split()
            if seen_natom and not patched and len(toks) == 9 and toks[0].isdigit():
                toks[5] = 'QQ9'          # a type no CHARMM release defines
                line = ' '.join(toks)
                patched = True
            out.append(line)
        self.assertTrue(patched, 'test setup failed to patch an atom type')
        bad_psf = Path('my_6pti_badtype.psf')
        bad_psf.write_text('\n'.join(out) + '\n')
        task_list = [{'continuation': dict(psf=str(bad_psf), pdb='my_6pti.pdb',
                                           coor='my_6pti.coor', xsc='my_6pti.xsc', vel='my_6pti.vel')}]
        self.controller.reconfigure_tasks(task_list)
        with self.assertRaises(PestiferBuildError) as ctx:
            self.controller.do_tasks()
        self.assertIn('QQ9', str(ctx.exception))
        bad_psf.unlink()

    def test_verify_parameters_false_skips_check(self):
        # The escape hatch pestifer's own internal continuations use: with verify_parameters=False
        # the same foreign-type PSF passes through untouched (no parse, no check).
        lines = Path('my_6pti.psf').read_text().splitlines()
        seen_natom = False
        out = []
        patched = False
        for line in lines:
            if '!NATOM' in line:
                seen_natom = True
                out.append(line)
                continue
            toks = line.split()
            if seen_natom and not patched and len(toks) == 9 and toks[0].isdigit():
                toks[5] = 'QQ9'
                line = ' '.join(toks)
                patched = True
            out.append(line)
        bad_psf = Path('my_6pti_badtype2.psf')
        bad_psf.write_text('\n'.join(out) + '\n')
        task_list = [{'continuation': dict(psf=str(bad_psf), pdb='my_6pti.pdb', coor='my_6pti.coor',
                                           xsc='my_6pti.xsc', vel='my_6pti.vel', verify_parameters=False)}]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        self.assertEqual(self.controller.tasks[0].result, 0)
        bad_psf.unlink()

    def test_termination_task(self):
        task_list = [
            {'continuation':
              dict(pdb='my_6pti.pdb', psf='my_6pti.psf', xsc='my_6pti.xsc', coor='my_6pti.coor', vel='my_6pti.vel')}, 
            {'terminate':dict(cleanup=False,chainmapfile='chainmap.yaml',package=dict(basename='my_system'))}]
        self.controller.reconfigure_tasks(task_list)
        self.assertEqual(len(self.controller.tasks), 2)
        self.assertIsInstance(self.controller.tasks[0], ContinuationTask)
        self.assertIsInstance(self.controller.tasks[1], TerminateTask)
        self.controller.do_tasks()
        task = self.controller.tasks[0]
        state = self.controller.pipeline.get_state_artifact(produced_by=task)
        # state: StateArtifacts = task.get_current_artifact('state')
        psf = state.psf
        self.assertEqual(psf.name, 'my_6pti.psf')
        pdb = state.pdb
        self.assertEqual(pdb.name, 'my_6pti.pdb')
        vel = state.vel
        self.assertEqual(vel.name, 'my_6pti.vel')
        xsc = state.xsc
        self.assertEqual(xsc.name, 'my_6pti.xsc')
        coor = state.coor
        self.assertEqual(coor.name, 'my_6pti.coor')
        self.assertTrue(psf.exists())
        self.assertTrue(pdb.exists())
        self.assertTrue(vel.exists())
        self.assertTrue(xsc.exists())
        self.assertTrue(coor.exists())
        psf.path.unlink()
        pdb.path.unlink()
        vel.path.unlink()
        xsc.path.unlink()
        coor.path.unlink()
        base_molecule: Molecule = task.get_current_artifact_data('base_molecule')
        self.assertIsInstance(base_molecule, Molecule)
        chainmapfile: YAMLFileArtifact = task.get_current_artifact('chainmapfile')
        self.assertIsInstance(chainmapfile, YAMLFileArtifact)
        self.assertTrue(Path('my_system.tar.gz').exists())
        self.assertTrue(chainmapfile.exists())
        Path('my_system.tar.gz').unlink()
        chainmapfile.path.unlink()

