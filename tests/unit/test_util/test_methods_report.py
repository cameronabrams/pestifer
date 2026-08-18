# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the drafted Methods section.

The value of this feature is that its statements are true, so these test the statements, not the
prose style: that a claim appears only when the record supports it, that a claim which would be
false does not appear, and that the LaTeX compiles.
"""

import json
import os
import tempfile
import unittest

from pestifer.util import methods_report as MR
from pestifer.core.run_record import RUN_RECORD_NAME, RUN_RECORD_VERSION, load_run_records


def _record(seed=27021972, atoms=14078, box=(55.9, 49.4, 50.0), protocol=None, citations=None):
    return {
        'run_record_version': RUN_RECORD_VERSION,
        'pestifer_version': '3.18.0',
        'title': 'test system',
        'config_file': 'demo.yaml',
        'seed': seed,
        'environment': {'charmmff': 'feb26',
                        'python': {'version': '3.13.13'},
                        'executables': {'namd3': {'version': 'NAMD 3.0.2 for Linux-x86_64-multicore'},
                                        'vmd': {'version': 'unknown'}}},
        'citations': {'entries': citations if citations is not None else [
            {'subject': 'NAMD', 'text': 'Phillips ...', 'doi': '10.1002/jcc.20289',
             'reason': '', 'key': 'namd2005'},
            {'subject': 'PDB 6PTI', 'text': 'doi:10.1016/x pmid:1', 'doi': '10.1016/x',
             'reason': 'the structure this build starts from', 'key': ''},
        ]},
        'system': {'atoms': atoms, 'total_charge': 0.0,
                   'box': {'a': box[0], 'b': box[1], 'c': box[2]}},
        'protocol': protocol if protocol is not None else [
            {'index': 2, 'task': 'md', 'ensemble': 'minimize', 'steps': 100},
        ],
    }


_ADAPTIVE_OK = [{'index': 5, 'task': 'density_equilibrate', 'ensemble': 'npt', 'steps': 82210,
                 'adaptive': True, 'converged': True, 'stopped_because': 'CONVERGED: ...'}]
_ADAPTIVE_CEILING = [{'index': 5, 'task': 'density_equilibrate', 'ensemble': 'npt', 'steps': 8000,
                      'adaptive': True, 'converged': False,
                      'stopped_because': 'CEILING: reached max_steps ...'}]


class TestTruthfulness(unittest.TestCase):
    """A Methods draft may omit; it may not misstate."""

    def test_converged_run_says_so(self):
        tex = MR.render_tex([_record(protocol=_ADAPTIVE_OK)])
        self.assertIn('until their convergence criterion was met', tex)
        self.assertNotIn('stopped at its step ceiling', tex)

    def test_unconverged_run_never_claims_convergence(self):
        # the defect this test exists for: a build that hit its ceiling must not be described
        # as having converged
        tex = MR.render_tex([_record(protocol=_ADAPTIVE_CEILING)])
        self.assertNotIn('until their convergence criterion was met', tex)
        self.assertIn('stopped at its step ceiling', tex)
        self.assertIn('TODO', tex)

    def test_steps_come_from_the_run_not_the_config(self):
        tex = MR.render_tex([_record(protocol=_ADAPTIVE_CEILING)])
        self.assertIn('8,000 steps', tex)

    def test_unknown_probe_never_reaches_prose(self):
        tex = MR.render_tex([_record()])
        self.assertNotIn('unknown', tex)

    def test_no_coordinate_citation_yields_a_todo_not_a_claim(self):
        tex = MR.render_tex([_record(citations=[])])
        self.assertIn('\\TODO{state the source', tex)

    def test_production_is_always_a_todo(self):
        tex = MR.render_tex([_record()])
        self.assertIn('did not run the production simulation', tex)

    def test_draft_banner_present(self):
        for text in (MR.render_tex([_record()]), MR.render_bib([_record()]),
                     MR.render_markdown([_record()])):
            self.assertIn(MR.DRAFT_BANNER, text)


class TestLatexValidity(unittest.TestCase):

    def test_task_underscores_are_escaped(self):
        tex = MR.render_tex([_record(protocol=_ADAPTIVE_CEILING)])
        # a bare density_equilibrate is a math subscript and will not compile
        for line in tex.splitlines():
            if line.startswith('%'):
                continue
            self.assertNotIn('density_equilibrate', line)
        self.assertIn('density\\_equilibrate', tex)

    def test_braces_balance(self):
        tex = MR.render_tex([_record(protocol=_ADAPTIVE_CEILING)])
        self.assertEqual(tex.count('{'), tex.count('}'))


class TestReplicaMerge(unittest.TestCase):

    def test_replicas_are_described_as_a_set(self):
        tex = MR.render_tex([_record(seed=1), _record(seed=2), _record(seed=3)])
        self.assertIn('3 independent replicas', tex)
        self.assertIn('1, 2, 3', tex)

    def test_single_run_states_its_seed(self):
        self.assertIn('The random number seed was 27021972', MR.render_tex([_record()]))

    def test_agreed_facts_are_stated_once(self):
        tex = MR.render_tex([_record(seed=1), _record(seed=2)])
        self.assertEqual(tex.count('14,078 atoms'), 1)

    def test_differing_boxes_are_reported_as_a_range_not_dropped(self):
        tex = MR.render_tex([_record(seed=1, box=(55.9, 49.4, 50.0)),
                             _record(seed=2, box=(56.2, 49.6, 50.3))])
        self.assertIn('ranged from', tex)
        self.assertIn('55.9', tex)
        self.assertIn('56.2', tex)

    def test_replica_step_counts_carry_a_caveat(self):
        tex = MR.render_tex([_record(seed=1, protocol=_ADAPTIVE_OK),
                             _record(seed=2, protocol=_ADAPTIVE_OK)])
        self.assertIn('Step counts are from replica 1', tex)


class TestBibliography(unittest.TestCase):

    def test_curated_entry_used_for_known_software(self):
        bib = MR.render_bib([_record()])
        self.assertIn('@article{namd2005', bib)
        self.assertIn('Scalable molecular dynamics', bib)

    def test_structure_entry_is_doi_only(self):
        # a PDB-sourced title cannot be capitalized correctly, so no metadata is invented
        bib = MR.render_bib([_record()])
        self.assertIn('@misc{pdb-6pti', bib)
        self.assertIn('10.1016/x', bib)

    def test_every_cited_key_exists_in_the_bib(self):
        records = [_record(protocol=_ADAPTIVE_OK)]
        tex, bib = MR.render_tex(records), MR.render_bib(records)
        import re
        cited = set()
        for m in re.finditer(r'\\cite\{([^}]*)\}', tex):
            cited.update(k.strip() for k in m.group(1).split(','))
        for key in cited:
            self.assertIn('{' + key + ',', bib, f'cited key {key} has no bib entry')

    def test_deduplicated_across_replicas(self):
        bib = MR.render_bib([_record(seed=1), _record(seed=2)])
        self.assertEqual(bib.count('@article{namd2005'), 1)


class TestWriteAndLoad(unittest.TestCase):

    def test_round_trip_through_disk(self):
        with tempfile.TemporaryDirectory() as td:
            rundir = os.path.join(td, 'rep-01')
            os.makedirs(rundir)
            with open(os.path.join(rundir, RUN_RECORD_NAME), 'w') as f:
                json.dump(_record(), f)
            records = load_run_records([rundir])
            written = MR.write_report(records, outdir=td, formats=('tex', 'bib', 'md'))
            self.assertEqual(len(written), 3)
            for p in written:
                self.assertTrue(os.path.getsize(p) > 0)


def _have(prog):
    import shutil
    return shutil.which(prog) is not None


@unittest.skipUnless(_have('pdflatex'), 'no pdflatex on PATH')
class TestGeneratedLatexCompiles(unittest.TestCase):
    """Compile the generated fragment for real.

    Brace-balance and escaping assertions are proxies; this is the thing they stand in for.  It
    earns its keep: the escaping tests passed while ``NAMD 3.0.2 for Linux-x86_64-multicore`` --
    an unescaped underscore in a *version string* rather than a task name -- made the fragment
    fail to compile.
    """

    DRIVER = (
        '\\documentclass{article}\n'
        '\\usepackage{url}\n'
        '\\newcommand{\\TODO}[1]{\\textbf{[TODO: #1]}}\n'
        '\\begin{document}\n'
        '\\input{methods}\n'
        '\\bibliographystyle{plain}\n'
        '\\bibliography{methods}\n'
        '\\end{document}\n'
    )

    def _build(self, records):
        import subprocess
        with tempfile.TemporaryDirectory() as td:
            MR.write_report(records, outdir=td, formats=('tex', 'bib'))
            with open(os.path.join(td, 'driver.tex'), 'w') as f:
                f.write(self.DRIVER)
            def run(cmd):
                return subprocess.run(cmd, cwd=td, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, text=True, timeout=180)
            p = run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'driver.tex'])
            self.assertEqual(p.returncode, 0, f'pdflatex failed:\n{p.stdout[-2500:]}')
            if _have('bibtex'):
                b = run(['bibtex', 'driver'])
                self.assertNotIn('Warning--to sort', b.stdout)
                run(['pdflatex', '-interaction=nonstopmode', 'driver.tex'])
                final = run(['pdflatex', '-interaction=nonstopmode', 'driver.tex'])
                self.assertNotIn('undefined', final.stdout.lower(),
                                 'a \\cite has no matching bib entry')
            self.assertTrue(os.path.exists(os.path.join(td, 'driver.pdf')))

    def test_single_run_compiles(self):
        self._build([_record(protocol=_ADAPTIVE_OK)])

    def test_ceilinged_run_compiles(self):
        # the branch carrying escaped task names and the extra TODO
        self._build([_record(protocol=_ADAPTIVE_CEILING)])

    def test_replica_set_compiles(self):
        self._build([_record(seed=1, box=(55.9, 49.4, 50.0), protocol=_ADAPTIVE_OK),
                     _record(seed=2, box=(56.2, 49.6, 50.3), protocol=_ADAPTIVE_OK)])

    def test_awkward_names_do_not_break_the_build(self):
        # underscores and percent signs reach the fragment from config filenames and versions
        rec = _record(protocol=_ADAPTIVE_CEILING)
        rec['config_file'] = 'my_system_v2.yaml'
        rec['environment']['executables']['namd3']['version'] = 'NAMD 3.0.2 for Linux-x86_64-multicore'
        rec['pestifer_version'] = '3.18.0_rc1'
        self._build([rec])
