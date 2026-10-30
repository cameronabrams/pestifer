# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""``md: stereo_restraints`` -- holding L alpha carbons L and trans peptide bonds trans through an md
stage with NAMD extraBonds.

The oracle for what the lines say is the PSF itself (each restrained quartet must name the right
atoms of one residue or one peptide bond), and for whether they reach the simulation it is NAMD's own
log: a file pestifer writes but NAMD never reads is the failure this guards against.
"""
import glob
import os
import re
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

from pestifer.psfutil.stereo import StereoChecker
from pestifer.tasks.mdtask import MDTask

HERE = Path(__file__).resolve().parent
PSF = HERE / 'fixtures' / 'continuation_inputs' / 'my_6pti.psf'
PDB = HERE / 'fixtures' / 'continuation_inputs' / 'my_6pti.pdb'
BPTI = HERE / 'fixtures' / 'bpti_build'


def parse(lines):
    out = []
    for ln in lines:
        f = ln.split()
        out.append((f[0], [int(x) for x in f[1:5]], float(f[5]), float(f[6])))
    return out


class TestRestraintLines(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.checker = StereoChecker(str(PSF))
        cls.atoms = cls.checker.topol.atoms.data

    def setUp(self):
        self.xyz = self.checker.load_coords(str(PDB))

    def test_every_center_and_trans_bond_of_a_real_build_is_held(self):
        rows = parse(self.checker.restraint_lines(self.xyz, 50.0))
        chiral = [r for r in rows if r[3] == -33.0]
        omega = [r for r in rows if r[3] == 180.0]
        ncenters = sum(1 for t in self.checker.centers if t[0] == 'D-CA')
        self.assertEqual(len(chiral), ncenters)
        self.assertEqual(len(omega), len(self.checker.peptides), 'BPTI has no cis bond')
        self.assertTrue(all(r[0] == 'dihedral' and r[2] == 50.0 for r in rows))

    def test_each_quartet_names_the_right_atoms(self):
        for kind, idx, _k, ref in parse(self.checker.restraint_lines(self.xyz)):
            a = [self.atoms[i] for i in idx]
            names = [x.atomname for x in a]
            if ref == -33.0:
                self.assertEqual(names, ['N', 'C', 'CB', 'CA'])
                self.assertEqual(len({(x.segname, x.resid.resid) for x in a}), 1, 'one residue')
            else:
                self.assertEqual(names, ['CA', 'C', 'N', 'CA'])
                self.assertEqual(len({x.segname for x in a}), 1)
                self.assertEqual(a[0].resid.resid, a[1].resid.resid)
                self.assertEqual(a[2].resid.resid, a[3].resid.resid)
                self.assertNotEqual(a[1].resid.resid, a[2].resid.resid, 'a peptide bond spans two residues')

    def test_a_center_already_d_is_left_alone(self):
        """NEGATIVE CONTROL: in a mirror image every center is D and no chirality restraint may be
        written (pulling one through the plane is not holding it), while omega, mirror-symmetric,
        is still held."""
        self.xyz[:, 0] *= -1
        rows = parse(self.checker.restraint_lines(self.xyz))
        self.assertEqual([r for r in rows if r[3] == -33.0], [])
        self.assertEqual(len(rows), len(self.checker.peptides))


class TestTheFile(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.d)
        self.task = MDTask.__new__(MDTask)
        self.task.basename = '00-02-000_md-minimize'
        self.task.taskname = 'md'
        self.registered = []
        self.task.register = lambda data, key, artifact_type=None: self.registered.append((data, key))
        self.state = NS(psf=NS(name=str(PSF)), pdb=NS(name=str(PDB)))

    def test_merged_after_an_ssrestraints_file(self):
        Path('ss.extrabonds').write_text('bond 0 1 10.0 1.5\n')
        name = self.task.make_stereo_restraints_file(self.state, 100.0, prior='ss.extrabonds')
        self.assertEqual(name, '00-02-000_md-minimize-stereo.extrabonds')
        lines = Path(name).read_text().splitlines()
        self.assertEqual(lines[0], 'bond 0 1 10.0 1.5', 'NAMD reads one file, so the prior comes first')
        self.assertTrue(all(ln.startswith('dihedral ') for ln in lines[1:]))
        self.assertEqual(self.registered, [('00-02-000_md-minimize-stereo', 'stereo_extrabonds')])

    def test_no_protein_writes_nothing(self):
        import pestifer.psfutil.stereo as stereo
        orig = stereo.StereoChecker.restraint_lines
        stereo.StereoChecker.restraint_lines = lambda self, xyz, k=100.0: []
        self.addCleanup(setattr, stereo.StereoChecker, 'restraint_lines', orig)
        self.assertIsNone(self.task.make_stereo_restraints_file(self.state, 100.0))
        self.assertEqual(glob.glob('*.extrabonds'), [])


@pytest.mark.needs_tools
class TestNAMDReadsThem(unittest.TestCase):
    """A short vacuum minimize of BPTI, with and without the option: NAMD must report reading
    exactly the dihedrals pestifer wrote, and none without it."""

    def setUp(self):
        self.d = tempfile.mkdtemp(dir=os.path.expanduser('~/devtests/pestifer')
                                  if os.path.isdir(os.path.expanduser('~/devtests/pestifer')) else None)
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.d)
        for f in ('00-01-00_psfgen-build.psf', '00-01-00_psfgen-build.pdb'):
            shutil.copy(BPTI / f, f)

    def _run(self, stereo):
        from pestifer.core.config import Config
        from pestifer.core.controller import Controller
        c = Controller().configure(Config().configure_new(), terminate=False)
        md = {'ensemble': 'minimize', 'nsteps': 20}
        if stereo:
            md['stereo_restraints'] = {'enabled': True, 'k': 100.0}
        c.reconfigure_tasks([{'continuation': {'psf': '00-01-00_psfgen-build.psf',
                                               'pdb': '00-01-00_psfgen-build.pdb'}},
                             {'md': md}])
        c.do_tasks()
        logs = [p for p in glob.glob('*md*.log') if 'vmd' not in p.lower()]
        self.assertTrue(logs, 'no NAMD log')
        text = Path(sorted(logs)[-1]).read_text()
        m = re.search(r'READ\s+(\d+)\s+EXTRA DIHEDRALS', text)
        return int(m.group(1)) if m else 0, glob.glob('*-stereo.extrabonds')

    def test_namd_reads_every_line(self):
        n, files = self._run(stereo=True)
        self.assertEqual(len(files), 1)
        written = sum(1 for ln in Path(files[0]).read_text().splitlines() if ln.startswith('dihedral'))
        self.assertGreater(written, 50)
        self.assertEqual(n, written)

    def test_off_by_default(self):
        n, files = self._run(stereo=False)
        self.assertEqual((n, files), (0, []))


if __name__ == '__main__':
    unittest.main()
