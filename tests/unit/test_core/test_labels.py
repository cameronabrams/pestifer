import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path

from pestifer.core import labels


class TestRegisterResnameSegtype(unittest.TestCase):
    """Source-editing of ``_segtypes`` in labels.py used by ``modify-package``."""

    def setUp(self):
        self.d = Path(tempfile.mkdtemp())
        self.tmp = self.d / 'labels.py'
        shutil.copy(labels.__file__, self.tmp)

    def _load(self):
        spec = importlib.util.spec_from_file_location('editedlabels', self.tmp)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_insert_registers_and_classifies(self):
        path, added = labels.register_resnames_segtype(['MYLIG', 'myl2'], 'ligand', labels_path=self.tmp)
        self.assertEqual(added, ['MYLIG', 'MYL2'])   # upper-cased
        m = self._load()
        self.assertIn('MYLIG', m._segtypes['ligand']['resnames'])
        self.assertIn('MYL2', m._segtypes['ligand']['resnames'])
        # the segtype map derived at import time picks it up
        self.assertEqual(m.Labels.segtype_of_resname['MYLIG'], 'ligand')
        # existing curated names are preserved
        self.assertIn('HEM', m._segtypes['ligand']['resnames'])

    def test_idempotent_for_existing_name(self):
        # HEM is a curated ligand; re-registering it is a no-op
        _, added = labels.register_resnames_segtype('HEM', 'ligand', labels_path=self.tmp)
        self.assertEqual(added, [])

    def test_unknown_segtype_raises(self):
        with self.assertRaises(ValueError):
            labels.register_resnames_segtype('MYLIG', 'nonesuch', labels_path=self.tmp)

    def test_result_still_parses(self):
        labels.register_resnames_segtype('NEWION', 'ion', labels_path=self.tmp)
        # loading the edited module is itself a parse+exec check
        m = self._load()
        self.assertIn('NEWION', m._segtypes['ion']['resnames'])

    def test_segtype_names_lists_keys(self):
        names = labels.segtype_names()
        self.assertIn('ligand', names)
        self.assertIn('ion', names)


if __name__ == '__main__':
    unittest.main()


class TestConfigSegtypeOverridesWin(unittest.TestCase):
    """`psfgen.segtypes` is an explicit classification, so it must win over pestifer's own -- derived
    or curated.  It used to apply only to names pestifer did not already know, so example 5's
    `other: [ACET, ACT]` changed ACT and silently left ACET as it was."""

    def setUp(self):
        import copy
        from pestifer.core import labels as L
        self._saved = (copy.deepcopy(L._segtypes), dict(L.Labels.segtype_of_resname))

    def tearDown(self):
        from pestifer.core import labels as L
        segtypes, mapping = self._saved
        L._segtypes.clear(); L._segtypes.update(segtypes)
        L.Labels.segtypes = L._segtypes
        L.Labels.segtype_of_resname.clear(); L.Labels.segtype_of_resname.update(mapping)

    def test_an_already_classified_residue_takes_the_configured_segtype(self):
        from pestifer.core.labels import Labels
        self.assertEqual(Labels.segtype_of_resname.get('ACET'), 'ligand')      # derived
        Labels.update_segtypes({'other': ['ACET', 'ACT']})
        self.assertEqual(Labels.segtype_of_resname['ACET'], 'other')
        self.assertEqual(Labels.segtype_of_resname['ACT'], 'other')

    def test_a_curated_name_moves_rather_than_being_listed_twice(self):
        from pestifer.core.labels import Labels
        self.assertIn('HEME', Labels.segtypes['ligand']['resnames'])
        Labels.update_segtypes({'cofactor': ['HEME']})
        self.assertEqual(Labels.segtype_of_resname['HEME'], 'cofactor')
        self.assertNotIn('HEME', Labels.segtypes['ligand']['resnames'])
        self.assertEqual(Labels.segtypes['cofactor']['resnames'].count('HEME'), 1)

    def test_auto_classification_of_unknown_residues_is_unchanged(self):
        from pestifer.core.labels import Labels
        Labels.update_segtypes({'ligand': ['ZZZ9']})
        self.assertEqual(Labels.segtype_of_resname['ZZZ9'], 'ligand')


class TestSugarAliasesAreComplete(unittest.TestCase):
    """A sugar reachable through a residue alias must also be classified as a glycan by its PDB code,
    and must have atom aliases for any atom whose PDB name CHARMM does not use.  A2G, RM4, BDP, IDR,
    SLB and ANE5 were aliased or aliasable with no glycan segtype; AGLCNA had a residue alias (from
    NDG) but none of the acetyl atom aliases BGLCNA had."""

    def test_every_aliased_sugar_code_is_a_glycan(self):
        from pestifer.core.labels import Labels, _residue_aliases
        for alias in _residue_aliases:
            pdb, charmm = alias.split()[:2]
            if Labels.segtype_of_resname.get(charmm) == 'glycan':
                with self.subTest(alias=alias):
                    self.assertEqual(Labels.segtype_of_resname.get(pdb), 'glycan')

    def test_n_acetyl_and_carboxylate_atoms_are_aliased_for_every_reachable_residue(self):
        from pestifer.core.labels import _residue_aliases, _atom_aliases
        reached = {a.split()[1] for a in _residue_aliases}
        atom_pairs = {tuple(a.split()[:3]) for a in _atom_aliases}
        needs = {'C7': 'C', 'O7': 'O', 'C8': 'CT', 'N2': 'N'}
        for res in reached & {'AGLCNA', 'BGLCNA', 'AGALNA', 'BGALNA'}:
            for pdb_atom, charmm_atom in needs.items():
                with self.subTest(residue=res, atom=pdb_atom):
                    self.assertIn((res, pdb_atom, charmm_atom), atom_pairs)
        for res in reached & {'AGLCA', 'BGLCA', 'AIDOA'}:
            for pdb_atom, charmm_atom in (('O6A', 'O61'), ('O6B', 'O62')):
                with self.subTest(residue=res, atom=pdb_atom):
                    self.assertIn((res, pdb_atom, charmm_atom), atom_pairs)


class TestResidueAliasesDoNotReclassifyCharmmResidues(unittest.TestCase):
    """A residue alias renames EVERY residue of that name that psfgen reads, including one pestifer
    wrote itself.  An alias from a PDB code that is also a CHARMM residue of a different segtype
    silently turns one molecule into another: "GLA AGAL" (alpha-galactose in the PDB) would have
    turned CHARMM's gamma-linolenic acid into a sugar.  Found 2026-09-14 before it shipped."""

    # Deliberate, same-segtype or already-curated collisions, each with its reason.
    KNOWN = {
        'PO4': 'phosphate ion -> H2PO4; both ions',
        'BGLC': 'CHARMM beta-glucose -> BGLCNA: restores a 6-char name truncated to 4 columns (see CLAUDE.md)',
        'GCU': 'also a CHARMM modified-RNA RESI; curated as glycan (glucuronic acid) before this test existed',
        'SIA': 'also a CHARMM modified-RNA RESI; curated as glycan (sialic acid) before this test existed',
    }

    def test_no_alias_moves_a_derived_charmm_residue_to_another_segtype(self):
        from pestifer.core.labels import _residue_aliases, _load_derived_segtypes, Labels
        derived = {r: st for st, names in _load_derived_segtypes().items() for r in names}
        for alias in _residue_aliases:
            pdb, charmm = alias.split()[:2]
            if pdb in self.KNOWN or pdb not in derived:
                continue
            with self.subTest(alias=alias):
                self.assertEqual(derived[pdb], Labels.segtype_of_resname.get(charmm),
                                 f'{alias} renames CHARMM RESI {pdb} ({derived[pdb]}) into a {Labels.segtype_of_resname.get(charmm)}')


class TestAtomselectMacrosAreGeneratedAndParseable(unittest.TestCase):
    """resources/tcl/macros.tcl drifted from the classification for months: it was generated from the
    curated lists, classification moved to a derived table, and regenerating would have shrunk the
    lipid and glycan keywords pestifer's own Tcl selects on.  When it was regenerated from the full
    classification, a derived name (SB3-10) made VMD reject the whole lipid macro and quietly fall
    back to its built-in keyword."""

    def test_macros_file_matches_the_generator(self):
        import os
        from pestifer.core.resourcemanager import ResourceManager
        RM = ResourceManager()
        path = os.path.join(RM.resource_path['tcl'], 'macros.tcl')
        self.assertEqual(open(path).read(), RM.atomselect_macros_text(),
                         'macros.tcl is stale: run `pestifer modify-package charmmff update-atomselect-macros`')

    def test_a_runtime_addition_never_reaches_the_generated_macros(self):
        # a config's psfgen.segtypes, or a residue auto-classified from ~/.pestifer/toppar, extends
        # the live tables; the package file must not depend on either (AP5 from a user toppar dir
        # once made this generator's output differ between a full suite run and a lone test)
        from pestifer.core.labels import Labels
        from pestifer.core.resourcemanager import ResourceManager
        Labels.update_segtypes({'ligand': ['ZZRT9']})
        self.assertEqual(Labels.segtype_of_resname['ZZRT9'], 'ligand')
        self.assertNotIn('ZZRT9', ResourceManager().atomselect_macros_text())

    def test_every_unquoted_name_is_one_vmd_can_parse(self):
        import os, re
        from pestifer.core.resourcemanager import ResourceManager
        text = open(os.path.join(ResourceManager().resource_path['tcl'], 'macros.tcl')).read()
        bodies = re.findall(r'update_atomselect_macro \w+ "resname ([^"]*)"', text)
        self.assertTrue(bodies)
        for body in bodies:
            bare = [t for t in body.split() if not (t.startswith("'") and t.endswith("'"))]
            bad = [t for t in bare if not re.fullmatch(r'[A-Za-z0-9_]+', t)]
            self.assertEqual(bad, [])
