# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
import os
from pestifer import resources
from pestifer.charmmff.charmmffcontent import CHARMMFFContent

import logging
logger = logging.getLogger(__name__)
from pathlib import Path

class TestCharmmffContent(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Pick the newest CHARMMFF version directory using the same
        # semantic month/year parse that ResourceManager uses in production
        # (mtime ordering is unreliable: directories can be re-extracted
        # out of release order, which silently picks the older release).
        logger.debug("Setting up TestCharmmffContent class...")
        from pestifer.core.resourcemanager import ResourceManager
        rm = ResourceManager()
        charmmff_path = rm.charmmff_version_dirs()[-1]
        cls.C = CHARMMFFContent(charmmff_path, force_rebuild=True)
        logger.debug("tarfilename: " + cls.C.tarfilename)
        logger.debug(f"Selected CHARMMFF version: {charmmff_path.name}")
        logger.debug("Done setting up TestCharmmffContent class...")

    def test_charmmffcontent_full_provisioning(self):
        self.C.deprovision()
        self.assertEqual(len(self.C.residues), 0)
        self.assertEqual(len(self.C.patches), 0)
        self.C.provision()
        self.assertEqual(len(self.C.residues), 3875)
        self.assertEqual(len(self.C.patches), 792)
        # the base repo ships 'lipid' and 'solvent'; a user may also have on-demand-cache
        # collections registered from ~/.pestifer/pdbrepository/<release>/ (e.g. a solvent box
        # built on the fly), so require at least the two base collections rather than exactly two
        self.assertGreaterEqual(len(self.C.pdbrepository.collections), 2)
        self.assertIn('lipid', self.C.pdbrepository.collections)
        self.assertIn('solvent', self.C.pdbrepository.collections)
        self.assertEqual(len(self.C.pdbrepository.collections['lipid'].info), 219)
        self.assertEqual(len(self.C.pdbrepository.collections['solvent'].info), 15)

    def test_charmmffcontent_restricted_provisioning(self):
        self.C.deprovision()
        self.assertEqual(len(self.C.residues), 0)
        self.assertEqual(len(self.C.patches), 0)
        self.C.provision(resnames=['ALA', 'TIP3', 'PSM', 'CHL1', 'DOPC', 'NNEU'])
        self.assertEqual(len(self.C.residues), 5)
        self.assertEqual(len(self.C.patches), 1)
        self.assertTrue('ALA' in self.C.residues)
        self.assertFalse('LYS' in self.C.residues)
        self.assertTrue('TIP3' in self.C.residues)
        self.assertTrue('TIP3' in self.C.pdbrepository)
        self.assertTrue('PSM' in self.C.residues)
        self.assertTrue('PSM' in self.C.pdbrepository)
        self.assertTrue('CHL1' in self.C.residues)
        self.assertTrue('DOPC' in self.C.residues)
        self.assertTrue('NNEU' in self.C.patches)

    def test_charmmffcontent_initialized(self):
        self.C.clean_local_charmmff_files()
        self.assertTrue(self.C != None)
        self.assertTrue(self.C.filenamemap != {})
        basenames = [k for k in self.C.filenamemap['top'].keys()]
        basenames.extend([k for k in self.C.filenamemap['toppar'].keys()])
        basenames.extend([k for k in self.C.filenamemap['par'].keys()])
        # feb26 ships both CGenFF v5 and v4.6 topology/parameter files; pestifer loads only v5
        # (v4.6 is a strict subset), so the two v4.6 files are excluded from the loaded set.
        # The count includes pestifer's custom/ additions (moreions.str, the degenerate-torsion
        # fills par, etc.).
        self.assertEqual(len(basenames), 57)
        self.assertEqual(len(set(basenames)), len(basenames))  # check for duplicates
        self.assertTrue(len(self.C.streams) > 0)
        self.assertEqual(self.C.streams.sort(), ['prot', 'carb', 'na', 'lipid'].sort())
        # self.assertIn('pestifer.top', self.C.custom_files)
        self.C.copy_charmmfile_local('par_all36m_prot.prm')
        self.assertTrue(os.path.exists('par_all36m_prot.prm'))
        self.C.copy_charmmfile_local('top_all36_prot.rtf')
        self.assertTrue(os.path.exists('top_all36_prot.rtf'))
        self.C.copy_charmmfile_local('toppar_water_ions.str')
        self.assertTrue(os.path.exists('toppar_water_ions.str'))
        with open('toppar_water_ions.str','r') as f:
            contents=f.read()
            self.assertTrue('! commented out by pestifer' in contents)
        self.C.copy_charmmfile_local('toppar_all36_carb_glycopeptide.str')
        self.assertTrue(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        with open('toppar_all36_carb_glycopeptide.str','r') as f:
            contents=f.read()
            self.assertFalse('! commented out by pestifer' in contents)
        self.C.copy_charmmfile_local('toppar_all36_moreions.str')
        self.assertTrue(os.path.exists('toppar_all36_moreions.str'))
        self.C.clean_local_charmmff_files()
        self.assertFalse(os.path.exists('par_all36m_prot.prm'))
        self.assertFalse(os.path.exists('top_all36_prot.rtf'))
        self.assertFalse(os.path.exists('toppar_water_ions.str'))
        self.assertFalse(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        self.assertFalse(os.path.exists('toppar_all36_moreions.str'))

    def test_charmmffcontent_pdbrepository_initialized(self):
        self.C.provision()
        self.assertTrue(self.C.pdbrepository!=None)
        self.assertTrue(self.C.pdbrepository.collections!=None)
        self.assertTrue(len(self.C.pdbrepository.collections)>0)
        self.assertTrue('lipid' in self.C.pdbrepository.collections)
        self.assertTrue('solvent' in self.C.pdbrepository.collections)
        self.assertTrue('TIP3' in self.C.pdbrepository)
        self.assertTrue('PSM' in self.C.pdbrepository)

    def test_charmmffcontent_pdbrepository_checkout(self):
        self.C.provision()
        c=self.C.pdbrepository.checkout('PSM')
        self.assertTrue(c!=None)
        self.assertEqual(c.get_charge(),0.0)
        params=c.get_parameters()
        self.assertIn('toppar_all36_lipid_sphingo.str',params)
        self.assertEqual(len(c.info['conformers']),10)
        self.assertAlmostEqual(c.info['conformers'][0]['head-tail-length'],26.57,places=2)
        self.assertAlmostEqual(c.info['conformers'][0]['max-internal-length'],30.54,places=2)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('PSM-00.pdb'))
        os.remove('PSM-00.pdb')

    def test_charmffcontent_get_topfile_of_patchname(self):
        """Test that the topfile of a patchname is returned correctly."""
        patchname = 'NNEU'
        topfile = self.C.resi_to_topfile_map[patchname]
        self.assertEqual(topfile, 'top_all36_prot.rtf')
        patchname = 'TYRO'
        topfile = self.C.resi_to_topfile_map[patchname]
        self.assertEqual(topfile, 'pestifer.top')

    def test_charmmffcontent_get_charge(self):
        self.C.provision()
        resname='POT'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,1.0)
        resname='CLA'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,-1.0)
        resname='PO4'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,-3.0)

    def test_charmmffcontent_get_resi(self):
        self.C.provision()
        self.assertEqual(len(self.C.residues),3875)
        self.assertEqual(len(self.C.patches),792)
        self.assertTrue('ALA' in self.C.residues)
        self.assertTrue('TIP3' in self.C.residues)
        self.assertTrue('FAKE' not in self.C.residues)
        self.assertTrue('NNEU' in self.C.patches)
        self.assertTrue('ASPP' in self.C.patches)
        self.assertTrue('GLUP' in self.C.patches)
        resname='ALA'
        ala=self.C.get_resi(resname)
        self.assertTrue(ala!=None)
        self.assertTrue(ala.metadata['streamID']=='prot')
        self.assertTrue(ala.metadata['substreamID']=='')
        self.assertAlmostEqual(ala.mass,71.0794,places=3)
        resname='TOCL1'
        TOCL1=self.C.get_resi(resname)
        self.assertTrue(TOCL1!=None)
        self.assertTrue(TOCL1.metadata['streamID']=='lipid')
        self.assertTrue(TOCL1.metadata['substreamID']=='cardiolipin')
        self.assertAlmostEqual(TOCL1.mass, 1457.02, places=2)

    def test_charmffcontent_detect_structure_CHL1(self):
        self.C.provision()
        resname='CHL1'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cholesterol')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['C3'])
        self.assertTrue(tails==['C27'])

    def test_charmffcontent_detect_structure_CHM1(self):
        self.C.provision()
        resname='CHM1'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cholesterol')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['C1'])
        self.assertTrue(tails==['C6'])

    def test_charmffcontent_detect_structure_DPPC(self):
        self.C.provision()
        resname='DPPC'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        shortest_paths=resi.annotation['shortest_paths']
        dist1=shortest_paths[heads[0]][tails[0]]
        dist2=shortest_paths[heads[0]][tails[1]]
        # logger.debug(f'dist1 {dist1} dist2 {dist2}')
        self.assertEqual(dist1,25)
        self.assertEqual(dist2,26)

    def test_charmffcontent_detect_structure_SDS(self):
        self.C.provision()
        resname='SDS'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='detergent')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['S'])
        self.assertTrue(tails==['C12'])

    def test_charmffcontent_detect_structure_TOCL1(self):
        self.C.provision()
        resname='TOCL1'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cardiolipin')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertEqual(len(tails),4)
        self.assertTrue(heads==['C2'])
        self.assertTrue(tails==['CA18', 'CB18', 'CC18', 'CD18'])

class TestAddCustomDirectory(unittest.TestCase):
    """
    Exercise :meth:`CHARMMFFContent.add_custom_directory` in isolation.

    These tests build a stub :class:`CHARMMFFContent` via ``__new__`` so we
    don't pay the full tarball-load cost just to verify directory-scanning
    semantics: tilde expansion, duplicate-basename handling, and duplicate
    ``RESI``/``PRES`` handling. All of those code paths only touch the
    file/RESI registry, not the residue/patch objects.
    """
    import tempfile, shutil

    def _make_stub(self):
        c = CHARMMFFContent.__new__(CHARMMFFContent)
        c.filenamemap = {'top': {}, 'par': {}, 'toppar': {}}
        c.all_topology_files = {}
        c.all_parameter_files = {}
        c.resi_to_topfile_map = {}
        c.user_custom_resnames = set()
        c.custom_files = []
        c.massdict = {}
        return c

    @staticmethod
    def _write_str(directory, basename, resname, comment=""):
        import os
        path = os.path.join(directory, basename)
        with open(path, "w") as f:
            f.write(f"* {comment}\nRESI {resname}      0.000\n")
        return path

    def test_tilde_expansion_in_searchpath(self):
        """``~/foo`` and ``$VAR/foo`` should be expanded before the
        :func:`os.path.isdir` check."""
        import os
        testdir = os.path.expanduser("~/.pestifer_unit_test_tilde")
        os.makedirs(testdir, exist_ok=True)
        try:
            self._write_str(testdir, "tt1.str", "TTT1")
            c = self._make_stub()
            # Should NOT raise NotADirectoryError on the literal '~/'.
            c.add_custom_directory("~/.pestifer_unit_test_tilde")
            self.assertIn("TTT1", c.resi_to_topfile_map)
        finally:
            import shutil
            shutil.rmtree(testdir)

    def test_duplicate_basename_warns_and_overrides(self):
        """Two searchpath entries with a same-basename file: the later one
        wins, and a warning is emitted (regression on the old AssertionError
        behavior)."""
        import os, tempfile
        d1 = tempfile.mkdtemp(prefix="ccdir_a_")
        d2 = tempfile.mkdtemp(prefix="ccdir_b_")
        try:
            path_a = self._write_str(d1, "DUP.str", "RA")
            path_b = self._write_str(d2, "DUP.str", "RB")
            c = self._make_stub()
            c.add_custom_directory(d1)
            with self.assertLogs(
                "pestifer.charmmff.charmmffcontent", level="WARNING"
            ) as caplog:
                c.add_custom_directory(d2)
            self.assertTrue(any("overrides earlier registration" in m
                                for m in caplog.output))
            # Last-loaded wins on the file registration
            self.assertEqual(c.filenamemap["toppar"]["DUP.str"], path_b)
        finally:
            import shutil
            shutil.rmtree(d1); shutil.rmtree(d2)

    def test_duplicate_resi_across_files_warns_and_overrides(self):
        """A ``RESI`` defined in two differently-named files triggers the
        RESI-override warning and last-wins. Same-basename clashes are
        handled by the basename warning and intentionally stay quiet here
        (one warning per real conflict)."""
        import os, tempfile
        d1 = tempfile.mkdtemp(prefix="ccdir_resi_a_")
        d2 = tempfile.mkdtemp(prefix="ccdir_resi_b_")
        try:
            self._write_str(d1, "from_a.str", "SAME")
            self._write_str(d2, "from_b.str", "SAME")
            c = self._make_stub()
            c.add_custom_directory(d1)
            with self.assertLogs(
                "pestifer.charmmff.charmmffcontent", level="WARNING"
            ) as caplog:
                c.add_custom_directory(d2)
            self.assertTrue(any("overrides earlier definition" in m
                                for m in caplog.output))
            self.assertEqual(c.resi_to_topfile_map["SAME"], "from_b.str")
        finally:
            import shutil
            shutil.rmtree(d1); shutil.rmtree(d2)
