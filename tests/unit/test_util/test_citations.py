# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the config-derived citation report.

No toolchain needed: the predicates read a config dict, so the tests hand them one.
"""

import os
import unittest
from unittest import mock

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


class _FakeStructureFile:
    """Stand-in for ``citations._structure_file``: callable, and a reusable context manager
    yielding a fixed path (or None)."""

    def __init__(self, value):
        self.value = value

    def __call__(self, _source_id):
        return self

    def __enter__(self):
        return self.value

    def __exit__(self, *exc):
        return False


class _FakeCid:
    def __init__(self, doi=None, pmid=None):
        self.doi, self.pmid = doi, pmid


class TestStructureSources(unittest.TestCase):

    def test_fetch_task_is_the_starting_structure(self):
        src = citations.structure_sources(_cfg(tasks=[{'fetch': {'sourceID': '6pti',
                                                                 'source': 'rcsb'}}]))
        self.assertEqual([(s.id, s.origin, s.citable) for s in src], [('6pti', 'rcsb', True)])

    def test_alphafold_model_is_not_citable(self):
        src = citations.structure_sources(_cfg(tasks=[{'fetch': {'sourceID': 'P22033',
                                                                 'source': 'alphafold'}}]))
        self.assertFalse(src[0].citable)

    def test_graft_donors_are_collected(self):
        cfg = _cfg(tasks=[{'psfgen': {'mods': {'grafts': ['A_1304:4b7i,C_1-8',
                                                          'B_1310:2wah,C_1-6']}}}])
        self.assertEqual(sorted(s.id for s in citations.structure_sources(cfg)), ['2wah', '4b7i'])

    def test_fusion_donors_are_collected(self):
        cfg = _cfg(tasks=[{'psfgen': {'mods': {'Cfusions': ['1ema:A:2-229,A']}}}])
        self.assertEqual([s.id for s in citations.structure_sources(cfg)], ['1ema'])

    def test_repeated_donor_reported_once(self):
        cfg = _cfg(tasks=[{'fetch': {'sourceID': '7xix'}},
                          {'psfgen': {'mods': {'grafts': ['A_1:4b7i,C_1-8',
                                                          'B_1:4b7i,C_1-8',
                                                          'C_1:4B7I,C_1-8']}}}])
        self.assertEqual([s.id for s in citations.structure_sources(cfg)], ['7xix', '4b7i'])

    def test_local_donor_file_is_not_citable(self):
        # a graft from a file the user made has no deposited paper behind it
        cfg = _cfg(tasks=[{'psfgen': {'mods': {'grafts': ['A_1:mydonor,C_1-8']}}}])
        self.assertFalse(citations.structure_sources(cfg)[0].citable)

    def test_unparsable_shortcode_is_skipped_not_fatal(self):
        cfg = _cfg(tasks=[{'psfgen': {'mods': {'grafts': ['total nonsense']}}}])
        self.assertEqual(citations.structure_sources(cfg), [])


class TestStructureCitations(unittest.TestCase):

    _cfg_one = staticmethod(lambda: _cfg(tasks=[{'fetch': {'sourceID': '6pti'}}]))

    def test_named_but_unresolved_without_pidibble_support(self):
        # an old pidibble must not make the report silently drop structures
        with mock.patch.object(citations, '_pidibble_import', return_value=object):
            cites, notes = citations.structure_citations(self._cfg_one())
        self.assertEqual([c.subject for c in cites], ['PDB 6PTI'])
        self.assertIn('unresolved', cites[0].text)
        self.assertTrue(any('1.10.0' in n for n in notes))

    def test_resolved_identifiers_are_reported(self):
        with mock.patch.object(citations, '_pidibble_import',
                                        return_value=type('P', (), {'citation_ids': 1})), \
             mock.patch.object(citations, '_structure_file',
                                        new=_FakeStructureFile('6pti.pdb')), \
             mock.patch.object(citations, '_pidibble_citation_ids',
                                        return_value=[_FakeCid('10.1/x', 12345)]):
            cites, notes = citations.structure_citations(self._cfg_one())
        self.assertEqual(cites[0].text, 'doi:10.1/x pmid:12345')
        self.assertEqual(notes, [])

    def test_unreadable_file_and_absent_citation_are_different_notes(self):
        api = type('P', (), {'citation_ids': 1})
        with mock.patch.object(citations, '_pidibble_import', return_value=api), \
             mock.patch.object(citations, '_structure_file',
                                        new=_FakeStructureFile('6pti.pdb')), \
             mock.patch.object(citations, '_pidibble_citation_ids', return_value=None):
            _, unreadable = citations.structure_citations(self._cfg_one())
        with mock.patch.object(citations, '_pidibble_import', return_value=api), \
             mock.patch.object(citations, '_structure_file',
                                        new=_FakeStructureFile('6pti.pdb')), \
             mock.patch.object(citations, '_pidibble_citation_ids', return_value=[]):
            _, none_listed = citations.structure_citations(self._cfg_one())
        self.assertIn('could not be read', unreadable[0])
        self.assertIn('lists no citation', none_listed[0])
        self.assertNotEqual(unreadable, none_listed)

    def test_missing_structure_file_is_noted(self):
        api = type('P', (), {'citation_ids': 1})
        with mock.patch.object(citations, '_pidibble_import', return_value=api), \
             mock.patch.object(citations, '_structure_file', new=_FakeStructureFile(None)):
            cites, notes = citations.structure_citations(self._cfg_one())
        self.assertEqual(cites, [])
        self.assertIn('neither in the working directory nor in a', notes[0])

    def test_noncitable_source_explained_not_dropped(self):
        cfg = _cfg(tasks=[{'fetch': {'sourceID': 'P22033', 'source': 'alphafold'}}])
        cites, notes = citations.structure_citations(cfg)
        self.assertEqual(cites, [])
        self.assertTrue(any('no deposited citation' in n for n in notes))

    def test_no_sources_means_no_notes(self):
        self.assertEqual(citations.structure_citations(_cfg()), ([], []))


class TestStructureFileLookup(unittest.TestCase):
    """Where the report finds the structure file it reads identifiers from.

    The ``terminate`` task sweeps every intermediate -- fetched structures included -- into
    ``<basename>-artifacts.tar.gz``, and the citation report runs after the run, so by the time
    it looks the originals are usually gone.  Reading the archive is not a nicety; without it the
    coordinate half of the report is empty for every build that ran to completion.
    """

    def setUp(self):
        import tempfile, os
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)

    def tearDown(self):
        import os
        os.chdir(self._cwd)
        self._tmp.cleanup()

    @staticmethod
    def _write_tarball(name, member):
        import tarfile, os
        os.makedirs('stage', exist_ok=True)
        with open(os.path.join('stage', os.path.basename(member)), 'w') as f:
            f.write('HEADER    TEST\n')
        with tarfile.open(name, 'w:gz') as tf:
            tf.add(os.path.join('stage', os.path.basename(member)), arcname=member)

    def test_prefers_a_file_still_in_the_working_directory(self):
        with open('6pti.pdb', 'w') as f:
            f.write('HEADER\n')
        with citations._structure_file('6pti') as path:
            self.assertEqual(path, '6pti.pdb')

    def test_falls_back_to_the_artifacts_archive(self):
        self._write_tarball('demo-artifacts.tar.gz', 'demo-artifacts/6pti.pdb')
        with citations._structure_file('6pti') as path:
            self.assertIsNotNone(path)
            self.assertTrue(path.endswith('6pti.pdb'))
            self.assertTrue(os.path.exists(path))

    def test_extracted_file_is_cleaned_up(self):
        self._write_tarball('demo-artifacts.tar.gz', 'demo-artifacts/6pti.pdb')
        with citations._structure_file('6pti') as path:
            extracted = path
        self.assertFalse(os.path.exists(extracted))

    def test_case_insensitive_match_in_the_archive(self):
        self._write_tarball('demo-artifacts.tar.gz', 'demo-artifacts/6PTI.pdb')
        with citations._structure_file('6pti') as path:
            self.assertIsNotNone(path)

    def test_cif_preferred_over_pdb(self):
        for name in ('6pti.pdb', '6pti.cif'):
            with open(name, 'w') as f:
                f.write('x\n')
        with citations._structure_file('6pti') as path:
            self.assertEqual(path, '6pti.cif')

    def test_absent_everywhere_yields_none(self):
        with citations._structure_file('6pti') as path:
            self.assertIsNone(path)

    def test_note_mentions_both_places_searched(self):
        api = type('P', (), {'citation_ids': 1})
        cfg = _cfg(tasks=[{'fetch': {'sourceID': '6pti'}}])
        with mock.patch.object(citations, '_pidibble_import', return_value=api):
            _, notes = citations.structure_citations(cfg)
        self.assertIn('artifacts.tar.gz', notes[0])
