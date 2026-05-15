# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for :mod:`pestifer.charmmff.ligand_paramgen.orchestrator`.
"""
import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock

from pestifer.charmmff.ligand_paramgen.orchestrator import (
    _drop_altlocs,
    find_unknown_ligand_residues,
)
from pestifer.molecule.atom import Hetatm
from pestifer.objs.resid import ResID


_serial = [0]


def _fake_lite(name, altloc=" ", occ=1.0):
    """Lightweight stand-in for the atom-level ``_drop_altlocs`` tests:
    we only need ``.name``, ``.altloc``, and ``.occ`` there."""
    return SimpleNamespace(name=name, altloc=altloc, occ=occ)


def _real_hetatm(name, altloc=" ", occ=1.0, chainID="A", resid_seq=1,
                 resname="LIG", elem="C", x=0.0, y=0.0, z=0.0):
    """Real pestifer Hetatm so it survives :class:`Residue` construction."""
    _serial[0] += 1
    return Hetatm(
        serial=_serial[0], name=name, altloc=altloc, resname=resname,
        chainID=chainID,
        resid=ResID(resseqnum=resid_seq, insertion=" "),
        x=x, y=y, z=z, occ=occ, beta=0.0, elem=elem, charge="0",
    )


class TestDropAltlocs(unittest.TestCase):
    """Verify altloc dedup rules used to clean up disordered ligand atoms
    like AP5 in 1AKE."""

    def test_no_altlocs_passes_through(self):
        atoms = [_fake_lite("PA"), _fake_lite("O1A"), _fake_lite("O2A")]
        self.assertEqual(_drop_altlocs(atoms), atoms)

    def test_higher_occupancy_wins(self):
        a_hi = _fake_lite("O1G", altloc="A", occ=0.7)
        a_lo = _fake_lite("O1G", altloc="B", occ=0.3)
        kept = _drop_altlocs([a_hi, a_lo])
        self.assertEqual(len(kept), 1)
        self.assertIs(kept[0], a_hi)

    def test_blank_altloc_beats_letter_at_equal_occ(self):
        blank = _fake_lite("OXT", altloc=" ", occ=0.5)
        letter = _fake_lite("OXT", altloc="A", occ=0.5)
        kept = _drop_altlocs([letter, blank])
        self.assertEqual(len(kept), 1)
        self.assertIs(kept[0], blank)

    def test_A_beats_B_at_equal_occ(self):
        a = _fake_lite("PD", altloc="A", occ=0.5)
        b = _fake_lite("PD", altloc="B", occ=0.5)
        kept = _drop_altlocs([a, b])
        self.assertEqual(len(kept), 1)
        self.assertIs(kept[0], a)

    def test_only_altloc_B_kept_when_alone(self):
        # If no blank/A exists, we still keep something rather than dropping
        # the atom entirely.
        b_only = _fake_lite("OXT", altloc="B", occ=1.0)
        kept = _drop_altlocs([b_only])
        self.assertEqual(kept, [b_only])

    def test_preserves_original_order(self):
        a = _fake_lite("PA", altloc=" ", occ=1.0)
        b_hi = _fake_lite("OB", altloc="A", occ=0.7)
        b_lo = _fake_lite("OB", altloc="B", occ=0.3)
        c = _fake_lite("PC", altloc=" ", occ=1.0)
        kept = _drop_altlocs([a, b_hi, b_lo, c])
        self.assertEqual([k.name for k in kept], ["PA", "OB", "PC"])


class TestFindUnknownLigandResidues(unittest.TestCase):
    """Exercise the resname filter + altloc dedup + representative-residue
    selection. PDB parsing itself is mocked away."""

    def _make_cc(self, known: set[str]):
        cc = MagicMock()
        cc.get_topfile_of_resname.side_effect = (
            lambda r: "known.str" if r in known else None
        )
        return cc

    def _parsed(self, hetatms):
        """Build a minimal pidibble-like PDBRecordDict stand-in. We bypass
        ``Hetatm(record)`` reconstruction by patching ``Hetatm`` in the
        orchestrator to a passthrough that returns each pre-built record
        unchanged."""
        return {"HETATM": hetatms}

    def _passthrough_patch(self):
        """Patch the Hetatm name in orchestrator's namespace so that
        ``Hetatm(rec)`` returns ``rec`` itself, but preserve ``_PDB_keyword``
        so the ``parsed.get(Hetatm._PDB_keyword, [])`` lookup still works."""
        from pestifer.charmmff.ligand_paramgen import orchestrator as orch
        passthrough = MagicMock(side_effect=lambda r: r)
        passthrough._PDB_keyword = "HETATM"
        return unittest.mock.patch.object(orch, "Hetatm", passthrough)

    def test_skips_resnames_known_to_charmmff(self):
        rec_known = _real_hetatm("O", resname="HOH", elem="O")
        rec_unknown = _real_hetatm("C1", resname="LIG", elem="C")
        with self._passthrough_patch():
            unknowns, skipped = find_unknown_ligand_residues(
                self._parsed([rec_known, rec_unknown]),
                self._make_cc(known={"HOH"}),
            )
        self.assertIn("HOH", skipped)
        self.assertIn("LIG", unknowns)
        self.assertNotIn("HOH", unknowns)

    def test_groups_first_residue_only(self):
        # Same resname at two distinct (chain, resid) keys -- we keep only
        # the first occurrence's atoms.
        records = [
            _real_hetatm("C1", resname="LIG", chainID="A", resid_seq=1),
            _real_hetatm("C2", resname="LIG", chainID="A", resid_seq=1),
            _real_hetatm("C1", resname="LIG", chainID="B", resid_seq=2),
            _real_hetatm("C2", resname="LIG", chainID="B", resid_seq=2),
        ]
        with self._passthrough_patch():
            unknowns, _ = find_unknown_ligand_residues(
                self._parsed(records), self._make_cc(known=set()),
            )
        self.assertIn("LIG", unknowns)
        self.assertEqual(unknowns["LIG"].chainID, "A")
        self.assertEqual(len(unknowns["LIG"].atoms.data), 2)

    def test_altlocs_are_dropped_in_representative(self):
        records = [
            _real_hetatm("PA", altloc=" ", resname="LIG"),
            _real_hetatm("O1A", altloc="A", occ=0.6, resname="LIG"),
            _real_hetatm("O1A", altloc="B", occ=0.4, resname="LIG"),
        ]
        with self._passthrough_patch():
            unknowns, _ = find_unknown_ligand_residues(
                self._parsed(records), self._make_cc(known=set()),
            )
        # PA + (one of O1A -- the altloc=A, higher-occ one)
        self.assertEqual(len(unknowns["LIG"].atoms.data), 2)
        names = [a.name for a in unknowns["LIG"].atoms.data]
        self.assertIn("PA", names)
        self.assertIn("O1A", names)

    def test_empty_pdb_returns_empty(self):
        unknowns, skipped = find_unknown_ligand_residues(
            {}, self._make_cc(known=set())
        )
        self.assertEqual(unknowns, {})
        self.assertEqual(skipped, [])


if __name__ == "__main__":
    unittest.main()
