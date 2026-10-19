# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The flat ``charmmff.files`` list, and the two views a build derives from it.

Pure list logic plus (for one test) parsing two shipped parameter files, so all of this
runs on CI -- no VMD, no NAMD.
"""
import os
import tempfile
import unittest

from pestifer.charmmff import fileset


#: What `scripters/psfgen.py` and `scripters/namd.py` built from
#: ``charmmff.{standard,custom}`` before the flat list replaced them, on 2026-10-19.
#: These are golden: the refactor was meant to change the *keywords*, not the force
#: field, and load order decides both terminal patches (topology) and which duplicate
#: parameter wins (last-wins merge).
LEGACY_TOPOLOGIES = [
    'top_all36_prot.rtf', 'top_all36_cgenff.rtf', 'top_all36_lipid.rtf',
    'top_all36_carb.rtf', 'top_all36_na.rtf',
    'toppar_water_ions.str', 'toppar_all36_carb_glycopeptide.str',
    'toppar_all36_carb_imlab.str', 'toppar_all36_prot_modify_res.str',
    'toppar_all36_moreions.str',
    'top_all35_ethers.rtf',                      # psfgen wants the ethers rtf last
]
LEGACY_PARAMETERS = [
    'par_all36m_prot.prm', 'par_all36_lipid.prm', 'par_all36_carb.prm',
    'par_all36_na.prm', 'par_all36_cgenff.prm',
    'toppar_water_ions.str', 'toppar_all36_carb_glycopeptide.str',
    'toppar_all36_carb_imlab.str', 'toppar_all36_prot_modify_res.str',
    'toppar_pestifer_dihedral_fills.prm', 'toppar_pestifer_conflict_resolutions.prm',
    'toppar_all36_moreions.str',
]


class TestDefaultSetIsUnchanged(unittest.TestCase):
    """The flat list must reproduce the pre-refactor load order exactly."""

    def test_topologies(self):
        self.assertEqual(fileset.topologies(fileset.resolve({})), LEGACY_TOPOLOGIES)

    def test_parameters(self):
        self.assertEqual(fileset.parameters(fileset.resolve({})), LEGACY_PARAMETERS)

    def test_the_golden_lists_are_not_trivially_equal(self):
        """A bug that returned the same list for both views would satisfy one assertion
        above by accident, so pin that the two views genuinely differ."""
        self.assertNotEqual(LEGACY_TOPOLOGIES, LEGACY_PARAMETERS)
        self.assertIn('par_all36m_prot.prm', LEGACY_PARAMETERS)
        self.assertNotIn('par_all36m_prot.prm', LEGACY_TOPOLOGIES)


class TestAddingAndReplacing(unittest.TestCase):

    def test_a_user_file_is_actually_loaded(self):
        """The defect this refactor fixes.  ``charmmff.user_custom.{rtf,str,prm}`` were
        declared in the schema, documented, and read by *no* scripter -- naming a file
        there validated, built, and loaded nothing.  A name in ``files`` must reach both
        derived views."""
        f = fileset.resolve({'files': ['myligand.str']})
        self.assertIn('myligand.str', fileset.topologies(f))
        self.assertIn('myligand.str', fileset.parameters(f))

    def test_an_added_file_loads_last(self):
        """Last-wins merging means an explicitly requested file must win, including over
        pestifer's own corrections."""
        f = fileset.resolve({'files': ['myligand.str']})
        self.assertEqual(fileset.parameters(f)[-1], 'myligand.str')
        self.assertGreater(fileset.parameters(f).index('myligand.str'),
                           fileset.parameters(f).index('toppar_pestifer_conflict_resolutions.prm'))

    def test_an_added_topology_loads_after_every_other_stream(self):
        """...with one exception carried over from the pre-refactor code: psfgen wants
        ``top_all35_ethers.rtf`` read last of all, so it stays last."""
        t = fileset.topologies(fileset.resolve({'files': ['myligand.str']}))
        self.assertEqual(t[-1], 'top_all35_ethers.rtf')
        self.assertEqual(t[-2], 'myligand.str')

    def test_a_prm_addition_reaches_parameters_but_not_topologies(self):
        f = fileset.resolve({'files': ['mine.prm']})
        self.assertIn('mine.prm', fileset.parameters(f))
        self.assertNotIn('mine.prm', fileset.topologies(f))

    def test_replace_defaults(self):
        f = fileset.resolve({'replace_defaults': True,
                             'files': ['top_all36_prot.rtf', 'par_all36m_prot.prm']})
        self.assertEqual(f, ['top_all36_prot.rtf', 'par_all36m_prot.prm'])
        self.assertEqual(fileset.topologies(f), ['top_all36_prot.rtf'])
        self.assertEqual(fileset.parameters(f), ['par_all36m_prot.prm'])

    def test_renaming_a_default_does_not_move_it(self):
        """Naming a file the defaults already carry is a no-op, not a reorder -- moving
        it to the end would silently change which duplicate wins."""
        self.assertEqual(fileset.resolve({'files': ['toppar_water_ions.str']}),
                         fileset.resolve({}))

    def test_a_non_charmm_filename_is_rejected(self):
        with self.assertRaises(ValueError):
            fileset.resolve({'files': ['notes.txt']})

    def test_add_is_idempotent_and_skips_defaults(self):
        cfg = {}
        fileset.add(cfg, 'toppar_water_ions.str', 'str')   # already a default
        self.assertEqual(cfg.get('files', []), [])
        fileset.add(cfg, 'extra.str', 'str')
        fileset.add(cfg, 'extra.str', 'str')
        self.assertEqual(cfg['files'], ['extra.str'])


class TestLegacyGroupsFold(unittest.TestCase):
    """Pre-2026-10-19 configs must keep working, and must fold to the same load set."""

    def test_standard_str_addition_folds(self):
        """What example 29 used to say."""
        cfg = {'standard': {'rtf': [], 'str': ['toppar_all36_lipid_sphingo.str'], 'prm': []}}
        fileset.fold_legacy_groups(cfg)
        self.assertEqual(cfg['files'], ['toppar_all36_lipid_sphingo.str'])
        self.assertIn('toppar_all36_lipid_sphingo.str',
                      fileset.topologies(fileset.resolve(cfg)))

    def test_the_untouched_default_groups_fold_to_nothing(self):
        """A config that merely carries the schema defaults for all three groups must
        produce no additions and no deprecation warning."""
        cfg = {
            'standard': {'rtf': ['top_all36_prot.rtf', 'top_all35_ethers.rtf',
                                 'top_all36_cgenff.rtf', 'top_all36_lipid.rtf',
                                 'top_all36_carb.rtf', 'top_all36_na.rtf'],
                         'str': ['toppar_water_ions.str',
                                 'toppar_all36_carb_glycopeptide.str',
                                 'toppar_all36_carb_imlab.str',
                                 'toppar_all36_prot_modify_res.str'],
                         'prm': ['par_all36m_prot.prm', 'par_all36_lipid.prm',
                                 'par_all36_carb.prm', 'par_all36_na.prm',
                                 'par_all36_cgenff.prm']},
            'custom': {'rtf': [], 'str': ['toppar_all36_moreions.str'],
                       'prm': ['toppar_pestifer_dihedral_fills.prm',
                               'toppar_pestifer_conflict_resolutions.prm']},
        }
        fileset.fold_legacy_groups(cfg)
        self.assertEqual(cfg.get('files', []), [])
        self.assertEqual(fileset.topologies(fileset.resolve(cfg)), LEGACY_TOPOLOGIES)
        self.assertEqual(fileset.parameters(fileset.resolve(cfg)), LEGACY_PARAMETERS)

    def test_user_custom_file_lists_now_load(self):
        """These keys had no consumer at all before the fold existed."""
        cfg = {'user_custom': {'searchpath': ['/some/dir'], 'rtf': [], 'prm': [],
                               'str': ['myligand.str'], 'segtypes': {'ligand': ['LIG']}}}
        fileset.fold_legacy_groups(cfg)
        self.assertIn('myligand.str', fileset.topologies(fileset.resolve(cfg)))
        self.assertEqual(cfg['searchpath'], ['/some/dir'])
        self.assertEqual(cfg['segtypes'], {'ligand': ['LIG']})

    def test_folding_twice_changes_nothing(self):
        """Config folds once and ResourceManager folds again; the second must be a
        no-op, not a duplicate append."""
        cfg = {'standard': {'str': ['toppar_all36_lipid_sphingo.str']}}
        fileset.fold_legacy_groups(cfg)
        first = list(cfg['files'])
        fileset.fold_legacy_groups(cfg)
        self.assertEqual(cfg['files'], first)

    def test_new_keys_win_over_legacy(self):
        cfg = {'searchpath': ['/new'], 'segtypes': {'cofactor': ['HEM']},
               'user_custom': {'searchpath': ['/old'], 'segtypes': {'ligand': ['LIG']}}}
        fileset.fold_legacy_groups(cfg)
        self.assertEqual(cfg['searchpath'], ['/new'])
        self.assertEqual(cfg['segtypes'], {'cofactor': ['HEM']})


class TestSearchpathExpansion(unittest.TestCase):

    def test_tilde_and_envvar(self):
        os.environ['PESTIFER_TEST_FFDIR'] = '/tmp/ffdir'
        try:
            got = fileset.searchpath({'searchpath': ['~/mine', '$PESTIFER_TEST_FFDIR/x']})
        finally:
            del os.environ['PESTIFER_TEST_FFDIR']
        self.assertEqual(got, [os.path.expanduser('~/mine'), '/tmp/ffdir/x'])


class TestParameterOrderIsLoadBearing(unittest.TestCase):
    """Why :data:`fileset.DEFAULT_FILES` is written in true load order rather than
    grouped by extension.

    ``CharmmParamFile.merge`` is last-wins.  ``toppar_pestifer_conflict_resolutions.prm``
    deliberately overrides a dihedral that ``toppar_all36_carb_imlab.str`` also defines,
    annotated in the file as *"CGenFF v5.0 value, chosen over imlab's analogy value"*.
    Put the pestifer corrections before the release streams -- which is what "all .prm,
    then all .str" does -- and that arbitration silently reverts.
    """

    QUARTET = ('NG2O1', 'CG2R61', 'CG2R61', 'NG2S3')
    IMLAB = 'toppar_all36_carb_imlab.str'
    FIXES = 'toppar_pestifer_conflict_resolutions.prm'

    @classmethod
    def setUpClass(cls):
        from pestifer.core.resourcemanager import ResourceManager
        cls.cc = ResourceManager(charmmff_config={}).charmmff_content
        cls.tmp = tempfile.TemporaryDirectory()
        cwd = os.getcwd()
        os.chdir(cls.tmp.name)
        try:
            # copy_charmmfile_local writes into the cwd; take absolute paths away with
            # us, because conftest re-chdirs between setUpClass and each test
            for n in (cls.IMLAB, cls.FIXES):
                cls.cc.copy_charmmfile_local(n)
            cls.path = {n: os.path.join(cls.tmp.name, n) for n in (cls.IMLAB, cls.FIXES)}
            for n, f in cls.path.items():
                assert os.path.exists(f), f'{n} was not copied to {f}'
        finally:
            os.chdir(cwd)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def _kchi(self, order):
        from pestifer.charmmff.charmmffprm import CharmmParamFile
        merged = CharmmParamFile()
        for n in order:
            merged.merge(CharmmParamFile.from_file(self.path[n]))
        hits = [d for d in merged.dihedrals
                if merged._dihedral_key(d)[0] in (self.QUARTET, tuple(reversed(self.QUARTET)))
                and d.n == 2]
        self.assertEqual(len(hits), 1, f'expected exactly one n=2 term, got {len(hits)}')
        return hits[0].Kchi

    def test_the_two_orders_disagree(self):
        pestifer_last = self._kchi([self.IMLAB, self.FIXES])
        imlab_last = self._kchi([self.FIXES, self.IMLAB])
        self.assertNotEqual(pestifer_last, imlab_last,
                            'the two files do not actually conflict on this quartet, so '
                            'this test proves nothing about ordering')
        self.assertAlmostEqual(pestifer_last, 1.25, places=4)
        self.assertAlmostEqual(imlab_last, 3.1, places=4)

    def test_the_default_set_puts_the_corrections_after_the_stream(self):
        params = fileset.parameters(fileset.resolve({}))
        self.assertLess(params.index(self.IMLAB), params.index(self.FIXES),
                        f'{self.FIXES} must load after {self.IMLAB} or its overrides revert')


class TestTheScripterActuallyEmitsAUserFile(unittest.TestCase):
    """The regression test for the defect, **at the level the defect lived**.

    ``charmmff.user_custom.{rtf,str,prm}`` were declared, documented, and read by no
    scripter: ``scripters/psfgen.py`` concatenated ``standard.rtf + custom.rtf +
    standard.str + custom.str`` and ``scripters/namd.py`` concatenated ``standard.prm +
    standard.str`` then ``custom.prm + custom.str``.  ``user_custom`` appeared in
    neither, so a user's file validated, built, and loaded nothing.

    Asserting on :mod:`pestifer.charmmff.fileset` alone would not catch a reintroduction
    -- that is one level above the call sites that were wrong.  This drives the real
    scripters, end to end from a config dict, and is the test to revert the
    ``fileset.resolve`` call in either scripter against.
    """

    STREAM = ('* a test ligand stream\n*\n\nRESI ZZTL   0.000\nGROUP\n'
              'ATOM C1 CG331 -0.270\n\nEND\n')

    @classmethod
    def setUpClass(cls):
        from pestifer.core.config import Config
        cls.tmp = tempfile.TemporaryDirectory()
        with open(os.path.join(cls.tmp.name, 'zztestlig.str'), 'w') as f:
            f.write(cls.STREAM)
        cls.cfg = Config(userdict={'charmmff': {'searchpath': [cls.tmp.name],
                                                'files': ['zztestlig.str']},
                                   'tasks': []}, quiet=True).configure_new()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def setUp(self):
        # fetch_standard_* copy every file in the set into the CURRENT directory, and
        # conftest puts that at this test module's own directory -- which would leave the
        # whole CHARMM set untracked under tests/.  Run in scratch instead.
        self._scratch = tempfile.TemporaryDirectory()
        self._cwd = os.getcwd()
        os.chdir(self._scratch.name)

    def tearDown(self):
        os.chdir(self._cwd)
        self._scratch.cleanup()

    def test_psfgen_scripter_loads_it(self):
        got = self.cfg.scripters['psfgen'].fetch_standard_charmm_topologies()
        self.assertIn('zztestlig.str', got,
                      'a file named in charmmff.files never reached psfgen')
        self.assertIn('top_all36_prot.rtf', got,
                      'the defaults vanished -- this would pass as two empty piles agreeing')

    def test_namd_scripter_loads_it(self):
        got = self.cfg.scripters['namd'].fetch_standard_charmm_parameters()
        self.assertIn('zztestlig.str', got,
                      'a file named in charmmff.files never reached NAMD')
        self.assertIn('par_all36m_prot.prm', got)

    def test_legacy_user_custom_str_loads_it_too(self):
        """The same file named the old way must now load, which it never did before."""
        from pestifer.core.config import Config
        cfg = Config(userdict={'charmmff': {'user_custom': {'searchpath': [self.tmp.name],
                                                            'str': ['zztestlig.str']}},
                               'tasks': []}, quiet=True).configure_new()
        self.assertIn('zztestlig.str',
                      cfg.scripters['psfgen'].fetch_standard_charmm_topologies())
