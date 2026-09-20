# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""What happens to atom serials past the 5-column PDB ceiling, pinned against the installed pidibble.

A PDB serial field is five columns.  Writers switch to hex at 100000 (``186a0``) and, past
``0xFFFFF`` = 1,048,575, give up and write ``*****``.  pidibble before 1.12.1 tested for ``*``
only *after* its hex branch -- and hex trips at atom 100000, long before the first marker -- so
re-reading any pestifer-written PDB of more than 1,048,575 atoms raised ``ValueError``.  That cost
a 1,579,027-atom membrane build its chain map, package and run record after 59 hours of MD.

These tests fail against pidibble 1.12.0, which is why ``pyproject.toml`` floors at 1.12.1.  They
also pin what the fix does *not* do: the serials are unrepresentable, not merely unparsed, so the
file re-reads with every atom past the ceiling at serial 0.  Losing them is survivable only
because nothing on the build's critical path needs them -- see ``TerminateTask.do``'s chain-map
guard, added for the same incident.
"""
import unittest

from pidibble.hex import AtomSerialParser

CEILING = 0xFFFFF   # 1,048,575: the largest serial five hex columns hold


class TestOverflowSerialsAreReadable(unittest.TestCase):

    def test_the_marker_parses_after_the_hex_switch(self):
        """The case that reaches a real file: hex has tripped long before the first marker."""
        p = AtomSerialParser()
        self.assertEqual(p('186a0'), 100000)          # the switch, at atom 100000
        self.assertEqual(p('*****'), 0)               # raised ValueError before pidibble 1.12.1

    def test_the_marker_parses_on_its_own(self):
        self.assertEqual(AtomSerialParser()('*****'), 0)

    def test_serials_below_the_ceiling_still_round_trip(self):
        p = AtomSerialParser()
        self.assertEqual(p('99999'), 99999)
        self.assertEqual(p(format(CEILING, 'x')), CEILING)

    def test_overflow_serials_are_lost_not_recovered(self):
        """Not a fix to want harder: five columns cannot hold these numbers.  Anything reading a
        PDB past the ceiling gets 0 for every atom beyond it, so it must not key on serials."""
        p = AtomSerialParser()
        p('186a0')
        self.assertEqual([p('*****'), p('*****')], [0, 0])
