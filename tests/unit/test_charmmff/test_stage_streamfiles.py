# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""``CHARMMFFContent.stage_streamfiles_locally``: keep the stream files that resolve, drop the rest.

A PSF written outside pestifer records its topology sources in ``REMARKS topology`` lines, and
those can name a file pestifer does not ship.  The recurring one is VMD's NAMD-specific
``toppar_water_ions_namd.str`` -- pestifer ships ``toppar_water_ions.str``, its equivalent.
``copy_charmmfile_local`` only *warns* on a miss and writes nothing, so passing the name onward
puts a ``topology <missing>`` line in the psfgen script and psfgen dies with "Unable to open
topology file / MOLECULE DESTROYED BY FATAL ERROR".

This function exists because that filter was written three times: ``continuation`` and ``merge``
each had their own copy, and ``make_membrane_system``'s prebuilt branch -- added later -- had
none, which is the bug reported 2026-09-30 from an Env membrane restart.  The tests below pin the
behaviour once, for all three.
"""
import os
import shutil
import tempfile
import unittest
from unittest import mock

from pestifer.charmmff.charmmffcontent import CHARMMFFContent

SHIPPED = 'toppar_water_ions.str'
NOT_SHIPPED = 'toppar_water_ions_namd.str'     # VMD's; the exact file that aborted a real build


class _CC:
    """A CHARMMFFContent stand-in: `copy_charmmfile_local` writes a file only for known names."""

    def __init__(self, resolvable):
        self.resolvable = set(resolvable)
        self.asked = []

    def copy_charmmfile_local(self, basename):
        self.asked.append(basename)
        if basename in self.resolvable:
            with open(basename, 'w') as fh:
                fh.write('* stub\n')
        return basename

    stage_streamfiles_locally = CHARMMFFContent.stage_streamfiles_locally


class TestStageStreamfilesLocally(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.d)

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.d, ignore_errors=True)

    def test_an_unresolvable_file_is_dropped_and_a_resolvable_one_is_kept(self):
        cc = _CC({SHIPPED})
        available, dropped = cc.stage_streamfiles_locally([SHIPPED, NOT_SHIPPED])
        self.assertEqual(available, [SHIPPED])
        self.assertEqual(dropped, [NOT_SHIPPED])

    def test_it_tries_every_name_before_deciding(self):
        """The decision is made on what landed on disk, not on a guess from the name."""
        cc = _CC({SHIPPED})
        cc.stage_streamfiles_locally([SHIPPED, NOT_SHIPPED])
        self.assertEqual(cc.asked, [SHIPPED, NOT_SHIPPED])

    def test_order_is_preserved(self):
        cc = _CC({'a.str', 'c.str'})
        available, _ = cc.stage_streamfiles_locally(['a.str', 'b.str', 'c.str'])
        self.assertEqual(available, ['a.str', 'c.str'])

    def test_everything_resolving_drops_nothing(self):
        cc = _CC({'a.str', 'b.str'})
        available, dropped = cc.stage_streamfiles_locally(['a.str', 'b.str'])
        self.assertEqual((available, dropped), (['a.str', 'b.str'], []))

    def test_nothing_resolving_returns_empty_rather_than_raising(self):
        """A build with no usable templates must continue: the caller opts into this helper
        precisely where the structure comes from the PSF, so an empty list is correct."""
        cc = _CC(set())
        available, dropped = cc.stage_streamfiles_locally([NOT_SHIPPED])
        self.assertEqual(available, [])
        self.assertEqual(dropped, [NOT_SHIPPED])

    def test_the_warning_names_the_file_and_the_caller(self):
        """The message has to say which file and which task, or a user cannot act on it."""
        cc = _CC(set())
        with mock.patch('pestifer.charmmff.charmmffcontent.logger') as log:
            cc.stage_streamfiles_locally([NOT_SHIPPED], context='make_membrane_system (prebuilt bilayer)')
        msg = str(log.warning.call_args[0][0])
        self.assertIn(NOT_SHIPPED, msg)
        self.assertIn('prebuilt bilayer', msg)
        self.assertIn('place the file in the run directory', msg)

    def test_a_resolvable_file_produces_no_warning(self):
        cc = _CC({SHIPPED})
        with mock.patch('pestifer.charmmff.charmmffcontent.logger') as log:
            cc.stage_streamfiles_locally([SHIPPED])
        log.warning.assert_not_called()
