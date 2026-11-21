# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A conformer's ID is its file number: ``conf: 7`` selects ``<RESI>-07.pdb``.

The ID is the conformer's position in its entry's ``info.yaml`` ``conformers`` list, and
``make_pdb_collection`` wrote that list in ``glob`` order -- the filesystem's.  So ``conf: 0``
named an arbitrary file, and a different one after each regeneration.  Found 2026-11-21: the
2026-10-08 collection's POPC "conformer 0" was ``POPC-04.pdb``, and the regenerated collection's
was another file entirely, which turned ``test_makemembranesystem_popc_pope`` (pinned to conf 0)
from a 33.5 A patch that passed into a 47.4 A patch that collapsed under NPT.

The reader now sorts the list (``pdbrepository._order_conformers``), which fixes collections
already shipped and cached; the writer sorts too, so new ``info.yaml`` files agree on disk.
"""
import os
import tarfile

import yaml

from pestifer.charmmff.pdbrepository import PDBCollection

_SHUFFLED = ['X-02.pdb', 'X-00.pdb', 'X-10.pdb', 'X-01.pdb']   # 10 after 02: numeric, not lexical


def _write_entry(root):
    """A one-entry collection ``<root>/stream/X`` whose info.yaml lists conformers out of order.
    Each pdb's content names its own file, so a lookup can be traced to the file it read."""
    d = os.path.join(root, 'stream', 'X')
    os.makedirs(d)
    for fn in _SHUFFLED:
        with open(os.path.join(d, fn), 'w') as f:
            f.write(f'REMARK {fn}\nEND\n')
    with open(os.path.join(d, 'info.yaml'), 'w') as f:
        yaml.safe_dump({'charge': 0.0,
                        'conformers': [{'pdb': fn, 'head-tail-length': float(i)}
                                       for i, fn in enumerate(_SHUFFLED)]}, f)
    return os.path.join(root, 'stream')


def _check(coll):
    entry = coll.contents['X']
    want = ['X-00.pdb', 'X-01.pdb', 'X-02.pdb', 'X-10.pdb']
    assert [c['pdb'] for c in entry.info['conformers']] == want
    for i, fn in enumerate(want):
        assert entry.pdbcontents[i].startswith(f'REMARK {fn}'), (i, entry.pdbcontents[i])
        assert entry.get_conformer_data(i)['pdb'] == fn
    # the metadata travels with its file: X-00's head-tail-length was written at list position 1
    assert entry.get_conformer_data(0)['head-tail-length'] == 1.0


def test_directory_collection_ids_follow_file_numbers(tmp_path):
    _check(PDBCollection.build_from_resources(_write_entry(str(tmp_path))))


def test_tarball_collection_ids_follow_file_numbers(tmp_path):
    src = _write_entry(str(tmp_path / 'src'))
    tgz = tmp_path / 'stream.tgz'
    with tarfile.open(tgz, 'w:gz') as tf:
        tf.add(src, arcname='stream')
    _check(PDBCollection.build_from_resources(str(tgz)))
