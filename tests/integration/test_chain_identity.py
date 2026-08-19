# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Regression test for chain identity surviving the PSF -> PDB round trip after MD.

A build writes a PDB from psfgen, runs NAMD (which returns binary coordinates carrying no
chain or segment information at all), and then regenerates a PDB by pairing those coordinates
back with the PSF.  Chain IDs live only in the PDB, so they have to be re-derived on every
regeneration -- and a bug there is silent.  The system still builds, still runs, and every atom
is in the right place; only its *identity* is wrong, which surfaces much later as a selection
that quietly matches nothing, or matches the wrong subunit.

This has been a real defect here ("chain IDs survive PSF->PDB regeneration after MD").

2WAH is used because its chain assignment is non-trivial: it has two protein chains, two glycan
segments, and two further chains, and the glycan segments deliberately *share* the chain ID of
the protein they are attached to (``AG01`` -> ``A``).  A regeneration that naively derived chain
from segname would produce a plausible-looking file with the glycans on their own chains, and no
other check in this suite would notice.

The comparison is against the psfgen-stage PDB from the same build, read out of the artifacts
archive ``terminate`` writes -- not against a hard-coded expectation, so the test states the
invariant (identity is *preserved*) rather than a snapshot that would need updating whenever the
chain-assignment policy legitimately changes.
"""

import collections
import glob
import os
import tarfile
import tempfile
import textwrap

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

from tests.integration.helpers import assert_psf_sane, parse_psf, parse_pdb_coords

pytestmark = pytest.mark.needs_tools


def _chain_of_line(line):
    return line[21]


def _segment_to_chains(psf_path, pdb_path):
    """{segname -> set of chain IDs} by pairing the PSF and PDB in atom order."""
    atoms, _bonds = parse_psf(psf_path)
    lines = [l for l in open(pdb_path) if l.startswith(('ATOM', 'HETATM'))]
    assert len(lines) == len(atoms), (
        f'{os.path.basename(pdb_path)}: {len(lines)} atoms, but its PSF declares {len(atoms)}')
    out = collections.defaultdict(set)
    for i, l in enumerate(lines, 1):
        out[atoms[i]['segname']].add(_chain_of_line(l))
    return dict(out)


def _extract_psfgen_stage(workdir, into):
    """The psfgen-stage PSF/PDB from the build's artifacts archive.

    ``terminate`` sweeps every intermediate into ``<basename>-artifacts.tar.gz``, so the file the
    build wrote before MD is there rather than loose in the run directory.
    """
    tarballs = glob.glob(os.path.join(workdir, '*-artifacts.tar.gz'))
    assert tarballs, f'no artifacts archive in {workdir}'
    with tarfile.open(tarballs[0]) as tf:
        names = tf.getnames()
        psf = next((n for n in names if n.endswith('psfgen-build.psf')), None)
        pdb = next((n for n in names if n.endswith('psfgen-build.pdb')), None)
        assert psf and pdb, f'psfgen-stage files not found in {tarballs[0]}'
        tf.extract(psf, path=into, filter='data')
        tf.extract(pdb, path=into, filter='data')
    return os.path.join(into, psf), os.path.join(into, pdb)


@pytest.mark.slow
def test_chain_ids_survive_regeneration_after_md(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = tmp_path / 'chains.yaml'
    cfg.write_text(textwrap.dedent("""
        title: chain identity regression (2wah, protein + glycans)
        tasks:
          - fetch:
              sourceID: 2wah
              source_format: pdb
          - psfgen:
              source:
                biological_assembly: 1
          - md:
              ensemble: minimize
              nsteps: 100
          - terminate:
              basename: chains
    """).lstrip())
    config = Config(userfile=str(cfg)).configure_new()
    report = Controller().configure(config).do_tasks()
    failed = {r['taskname']: r['result'] for r in report.values() if r['result'] != 0}
    assert not failed, f'build task(s) failed: {failed}'

    final_psf, final_pdb = tmp_path / 'chains.psf', tmp_path / 'chains.pdb'
    after = _segment_to_chains(final_psf, final_pdb)

    with tempfile.TemporaryDirectory() as td:
        stage_psf, stage_pdb = _extract_psfgen_stage(str(tmp_path), td)
        before = _segment_to_chains(stage_psf, stage_pdb)

    # (a) nothing lost its chain ID in regeneration
    blanks = [s for s, cs in after.items() if ' ' in cs or '' in cs]
    assert not blanks, (
        f'segments {sorted(blanks)} have blank chain IDs after MD -- chain identity was dropped '
        f'when the PDB was regenerated from the PSF')

    # (b) each segment still carries a single, unambiguous chain ID
    split = {s: sorted(cs) for s, cs in after.items() if len(cs) != 1}
    assert not split, f'segments span more than one chain ID after MD: {split}'

    # (c) the assignment is the same one psfgen made.  This is the actual regression: glycan
    #     segments share their parent protein's chain ID, so a regeneration that derived chain
    #     from segname would pass (a) and (b) and still be wrong.
    assert after == before, (
        f'chain assignment changed across MD.\n  psfgen stage: '
        f'{ {k: sorted(v) for k, v in sorted(before.items())} }\n  after MD:     '
        f'{ {k: sorted(v) for k, v in sorted(after.items())} }')

    # guard: this structure must actually exercise the shared-chain-ID case
    shared = [c for c, n in collections.Counter(
        next(iter(cs)) for cs in after.values()).items() if n > 1]
    assert shared, (
        'no chain ID is shared by two segments in this build, so it would not detect a '
        'regeneration that derived chain IDs from segment names')

    assert_psf_sane(final_psf, final_pdb, unminimized=True, context='chain identity')
