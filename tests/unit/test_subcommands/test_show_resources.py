# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Tests for show-resources resname formatting (pestifer.subcommands.show_resources)."""
import unittest

from pestifer.subcommands.show_resources import _report_resname, _report_resname_compact


def _base(**over):
    info = {
        'resname': 'BORO', 'in_topology': True, 'kind': 'residue (RESI)', 'is_patch': False,
        'topfile': 'top_all36_cgenff.rtf', 'source': 'standard', 'segtype': 'ligand',
        'charmm_alias': None, 'charmm_synonym': None, 'in_pdbrepository': False,
        'longname': None, 'nconformers': 0, 'charge': None, 'head_tail_length': None,
    }
    info.update(over)
    return info


class TestReportResnameSynonym(unittest.TestCase):

    def _capture(self, fn, info):
        lines = []
        fn(info, out_stream=lines.append)
        return '\n'.join(lines)

    def test_block_shows_synonym_when_present(self):
        out = self._capture(_report_resname,
                            _base(charmm_synonym='B1O2C1H5, methyl boronic acid, neutral'))
        self.assertIn('methyl boronic acid', out)
        # rendered on its own quoted line under the topology line
        self.assertIn('"B1O2C1H5, methyl boronic acid, neutral"', out)

    def test_block_omits_synonym_when_absent(self):
        out = self._capture(_report_resname, _base(charmm_synonym=None))
        self.assertNotIn('"', out)

    def test_compact_shows_synonym(self):
        out = self._capture(_report_resname_compact, _base(charmm_synonym='phenol, adm jr.'))
        self.assertIn('BORO', out)
        self.assertTrue(out.rstrip().endswith('phenol, adm jr.'))

    def test_compact_without_synonym_has_no_trailing_text(self):
        out = self._capture(_report_resname_compact, _base(charmm_synonym=None))
        self.assertTrue(out.rstrip().endswith('pdb: -'))


if __name__ == '__main__':
    unittest.main()
