# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the config-derived citation report.

No toolchain needed: the predicates read a config dict, so the tests hand them one.
"""

import unittest

from pestifer.util import citations
from pestifer.util.citations import citations_for, log_citations


def _cfg(tasks=None, namd_version=3, user_custom=None):
    """Minimal config-shaped mapping: only what the predicates read."""
    return {'user': {
        'tasks': tasks or [],
        'namd': {'namd-version': namd_version},
        'charmmff': {'user_custom': user_custom or {}},
        'psfgen': {'segtypes': {}, 'aliases': {}},
    }}


def _subjects(cfg):
    return [c.subject for c in citations_for(cfg)]


class TestUnconditional(unittest.TestCase):

    def test_core_software_always_cited(self):
        got = _subjects(_cfg())
        for expected in ('NAMD', 'VMD', 'psfgen', 'CHARMM36m (protein)'):
            self.assertIn(expected, got)

    def test_namd3_cited_by_default(self):
        self.assertIn('NAMD 3', _subjects(_cfg()))

    def test_namd2_build_does_not_cite_namd3(self):
        self.assertNotIn('NAMD 3', _subjects(_cfg(namd_version=2)))

    def test_a_bare_build_owes_nothing_conditional(self):
        cfg = _cfg(tasks=[{'fetch': {}}, {'psfgen': {}}, {'solvate': {}}, {'md': {}}])
        self.assertEqual([c for c in citations_for(cfg) if c.subject not in
                          ('NAMD', 'VMD', 'psfgen', 'CHARMM36m (protein)',
                           'CHARMM force field', 'NAMD 3')], [])


class TestConditional(unittest.TestCase):

    def test_membrane_build_cites_lipid_parameters(self):
        cfg = _cfg(tasks=[{'make_membrane_system': {}}])
        self.assertIn('CHARMM36 (lipids)', _subjects(cfg))

    def test_membrane_equilibrate_alone_also_counts(self):
        self.assertIn('CHARMM36 (lipids)', _subjects(_cfg(tasks=[{'membrane_equilibrate': {}}])))

    def test_grafts_cite_carbohydrate_parameters(self):
        cfg = _cfg(tasks=[{'psfgen': {'mods': {'grafts': ['A_1304:4b7i,C_1-8']}}}])
        self.assertIn('CHARMM36 (carbohydrates)', _subjects(cfg))

    def test_schema_default_glycan_block_does_not_trigger_carbohydrates(self):
        # the schema materializes source.sequence.glycans on every psfgen task, so its mere
        # presence says nothing about this build and must not pull in a citation
        cfg = _cfg(tasks=[{'psfgen': {'source': {'sequence': {'glycans': {'declash':
                                                                         {'maxcycles': 500}}}},
                                      'mods': {'grafts': []}}}])
        self.assertNotIn('CHARMM36 (carbohydrates)', _subjects(cfg))

    def test_nonwater_solvent_cites_cgenff(self):
        cfg = _cfg(tasks=[{'solvate': {'solvent': 'DMSO'}}])
        self.assertIn('CGenFF', _subjects(cfg))

    def test_water_solvent_does_not_cite_cgenff(self):
        for solvent in ('TIP3', 'tip3', 'TIP3P'):
            self.assertNotIn('CGenFF', _subjects(_cfg(tasks=[{'solvate': {'solvent': solvent}}])))

    def test_default_solvate_is_water(self):
        self.assertNotIn('CGenFF', _subjects(_cfg(tasks=[{'solvate': {}}])))

    def test_user_custom_parameters_cite_cgenff_and_the_program(self):
        got = _subjects(_cfg(user_custom={'searchpath': ['/some/dir']}))
        self.assertIn('CGenFF', got)
        self.assertEqual(got.count('CGenFF program'), 2)   # atom typing + charges

    def test_pdb2pqr_task_cites_pdb2pqr_apbs_and_propka(self):
        got = _subjects(_cfg(tasks=[{'pdb2pqr': {'pH': 7.0}}]))
        for expected in ('PDB2PQR', 'APBS/PDB2PQR', 'PROPKA'):
            self.assertIn(expected, got)

    def test_conditional_citations_explain_themselves(self):
        cfg = _cfg(tasks=[{'make_membrane_system': {}}])
        lipid = next(c for c in citations_for(cfg) if c.subject == 'CHARMM36 (lipids)')
        self.assertIn('make_membrane_system', lipid.reason)


class TestRendering(unittest.TestCase):

    def test_render_includes_doi_when_known(self):
        self.assertIn('doi:10.1002/jcc.20289', citations.NAMD.render())

    def test_render_omits_doi_when_absent(self):
        self.assertNotIn('doi:', citations.PSFGEN.render())

    def test_log_block_is_uniformly_prefixed(self):
        lines = []
        log_citations(_cfg(), log=lines.append)
        self.assertTrue(lines and all(l.startswith('citations: ') for l in lines))

    def test_log_never_raises_on_a_broken_config(self):
        # a citation report must not be the thing that fails a finished build; an unreadable
        # config degrades to the unconditional citations rather than to nothing, since NAMD,
        # VMD, psfgen and CHARMM are owed regardless of what the config said
        got = log_citations(None, log=lambda _: None)
        self.assertEqual([c.subject for c in got],
                         [c.subject for c in citations.ALWAYS] + ['NAMD 3'])

    def test_unreadable_config_yields_only_unconditional(self):
        self.assertEqual(_subjects({}), [c.subject for c in citations.ALWAYS] + ['NAMD 3'])
