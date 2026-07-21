# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for PsfgenTask behavior when following a continuation task.

The slow test here verifies that DISU patch REMARK records written by a first
psfgen build (vacuum BPTI, mirroring example 1) survive into the PSF produced
by a second psfgen run that loads the first build via a continuation task and
applies a single non-cysteine mutation.

Requires VMD; run with --runslow.
"""
import shutil
import unittest
from pathlib import Path

import pytest

from pestifer.core.artifacts import StateArtifacts
from pestifer.core.config import Config
from pestifer.core.controller import Controller

pytestmark = pytest.mark.needs_tools

# The three disulfide bonds present in 6pti (chain A): C5-C55, C14-C38, C30-C51
_BPTI_DISU_REMARKS = {
    'patch DISU A:5 A:55',
    'patch DISU A:14 A:38',
    'patch DISU A:30 A:51',
}


def _psf_patch_remarks(psf_path: str) -> set[str]:
    """Return normalised patch remark strings from the header of a psfgen PSF."""
    remarks = set()
    with open(psf_path) as fh:
        for line in fh:
            line = line.strip()
            if '!NATOM' in line:
                break
            if line.startswith('REMARKS patch'):
                remarks.add(' '.join(line[len('REMARKS '):].split()))
    return remarks


class TestPsfgenContinuation(unittest.TestCase):

    @pytest.mark.slow
    def test_disu_preserved_after_continuation_psfgen(self):
        """
        Stage 1: fetch 6pti and build a vacuum BPTI PSF/PDB (no mods, no
        solvation — mirrors the psfgen step of example 1).
        Stage 2: load stage-1 output via a continuation task, then re-run
        psfgen with a single non-cysteine mutation (THR11 → ALA).
        The three DISU patch REMARKs must survive into the stage-2 PSF.
        """
        before = set(Path('.').iterdir())
        success = False
        try:
            # ── Stage 1: fetch 6pti + psfgen (vacuum, no mods) ───────────────
            c1 = Controller().configure(Config().configure_new(), terminate=False)
            c1.reconfigure_tasks([
                {'fetch': {'sourceID': '6pti'}},
                {'psfgen': {}},
            ])
            c1.do_tasks()

            state1: StateArtifacts = c1.tasks[-1].get_current_artifact('state')
            stage1_psf = state1.psf.name
            stage1_pdb = state1.pdb.name

            # Sanity-check: stage-1 PSF must carry all three DISU remarks
            stage1_remarks = _psf_patch_remarks(stage1_psf)
            for r in _BPTI_DISU_REMARKS:
                self.assertIn(r, stage1_remarks,
                              f'stage 1 PSF is missing expected remark: {r!r}')

            # ── Stage 2: continuation → psfgen with one non-cysteine mutation ─
            # Use controller index 1 so basenames do not collide with stage 1.
            c2 = Controller().configure(
                Config().configure_new(), index=1, terminate=False)
            c2.reconfigure_tasks([
                {'continuation': {'psf': stage1_psf, 'pdb': stage1_pdb}},
                {'psfgen': {'mods': {'mutations': ['A:T11A']}}},
            ])
            c2.do_tasks()

            state2: StateArtifacts = c2.tasks[-1].get_current_artifact('state')
            stage2_psf = state2.psf.name

            # Primary assertion: all three DISU remarks must survive stage 2
            stage2_remarks = _psf_patch_remarks(stage2_psf)
            for r in _BPTI_DISU_REMARKS:
                self.assertIn(r, stage2_remarks,
                              f'stage 2 PSF lost DISU remark after '
                              f'continuation + psfgen: {r!r}')

            success = True
        finally:
            if success:
                after = set(Path('.').iterdir())
                for p in after - before:
                    if p.is_file():
                        p.unlink()
                    elif p.is_dir():
                        shutil.rmtree(p)


class TestAssertCoordinatesSet(unittest.TestCase):
    """PsfgenTask._assert_coordinates_set hard-errors when psfgen leaves atoms at the origin
    (e.g. a monatomic ion whose PDB resname collides with a multi-atom CHARMM/CGenFF residue)."""

    def _task_for_pdb(self, pdb_path):
        from pestifer.tasks.psfgen import PsfgenTask
        from types import SimpleNamespace
        t = PsfgenTask.__new__(PsfgenTask)
        state = SimpleNamespace(pdb=SimpleNamespace(name=pdb_path))
        t.get_current_artifact = lambda key: state
        return t

    def _write_pdb(self, atoms):
        # standard PDB columns: x starts at index 30 (col 31), matching what the check reads
        import tempfile, os
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as f:
            for i, (resname, chain, resid, name, xyz) in enumerate(atoms, start=1):
                x, y, z = xyz
                f.write(f'ATOM  {i:>5} {name:<4} {resname:>3} {chain}{resid:>4}    '
                        f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      A    C\n')
        self.addCleanup(os.remove, path)
        return path

    def test_no_origin_atoms_passes(self):
        pdb = self._write_pdb([
            ('ALA', 'A', 1, 'CA', (10.0, 11.0, 12.0)),
            ('CAL', 'B', 401, 'CAL', (5.0, 6.0, 7.0)),
        ])
        self._task_for_pdb(pdb)._assert_coordinates_set()  # must not raise

    def test_origin_atoms_hard_error_names_residue_and_fix(self):
        from pestifer.core.errors import PestiferBuildError
        pdb = self._write_pdb([
            ('ALA', 'A', 1, 'CA', (10.0, 11.0, 12.0)),
            ('CA', 'A', 401, 'C3', (0.0, 0.0, 0.0)),
            ('CA', 'A', 401, 'O3', (0.0, 0.0, 0.0)),
        ])
        with self.assertRaises(PestiferBuildError) as ctx:
            self._task_for_pdb(pdb)._assert_coordinates_set()
        msg = str(ctx.exception)
        self.assertIn('origin', msg)
        self.assertIn('CA A401', msg)        # names the offending residue
        self.assertIn('CA CAL', msg)         # shows the alias fix


class TestCheckIonAliases(unittest.TestCase):
    """PsfgenTask._check_ion_aliases hard-errors early (before psfgen) when a monatomic ion's PDB
    resname needs a CHARMM alias that isn't configured, but never fires when it is aliased, when
    the residue is multi-atom (a genuine organic residue of the same name), or for non-ions."""

    def _task(self, residues, residue_aliases=()):
        from pestifer.tasks.psfgen import PsfgenTask
        from types import SimpleNamespace
        t = PsfgenTask.__new__(PsfgenTask)
        au = SimpleNamespace(residues=SimpleNamespace(data=residues))
        t.base_molecule = SimpleNamespace(asymmetric_unit=au)
        t.scripters = {'psfgen': SimpleNamespace(
            psfgen_config={'aliases': {'residue': list(residue_aliases)}})}
        return t

    def _res(self, resname, natoms=1, chainID='B', resid=401):
        from types import SimpleNamespace
        return SimpleNamespace(resname=resname, atoms=[0] * natoms,
                               chainID=chainID, resid=SimpleNamespace(resid=resid))

    def test_monatomic_ion_unaliased_hard_errors(self):
        from pestifer.core.errors import PestiferBuildError
        t = self._task([self._res('ALA', natoms=5), self._res('CA', natoms=1, resid=401)])
        with self.assertRaises(PestiferBuildError) as ctx:
            t._check_ion_aliases()
        msg = str(ctx.exception)
        self.assertIn('CA (calcium) at B401', msg)
        self.assertIn('"CA CAL"', msg)
        self.assertIn('"CAL CA CAL"', msg)

    def test_aliased_ion_passes(self):
        t = self._task([self._res('CA', natoms=1)], residue_aliases=['CA CAL'])
        t._check_ion_aliases()  # must not raise

    def test_multiatom_same_name_residue_passes(self):
        # a genuine 68-atom CGenFF residue named CA is not a bare ion -> no error
        self._task([self._res('CA', natoms=68)])._check_ion_aliases()

    def test_non_ion_passes(self):
        self._task([self._res('HOH', natoms=1), self._res('LIG', natoms=20)])._check_ion_aliases()

    def test_multiple_ions_all_reported(self):
        from pestifer.core.errors import PestiferBuildError
        t = self._task([self._res('CA', natoms=1, resid=401), self._res('ZN', natoms=1, resid=500)])
        with self.assertRaises(PestiferBuildError) as ctx:
            t._check_ion_aliases()
        msg = str(ctx.exception)
        self.assertIn('CA (calcium)', msg)
        self.assertIn('ZN (zinc)', msg)
        self.assertIn('"ZN ZN2"', msg)
