# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for :mod:`pestifer.charmmff.ligand_paramgen.rcsb`.

HTTP calls to ``data.rcsb.org`` are mocked so the suite stays offline-clean.
"""
import unittest
from unittest.mock import MagicMock, patch

import requests

from pestifer.charmmff.ligand_paramgen.rcsb import (
    RCSBLookupError,
    fetch_ligand_smiles,
)


def _make_response(status_code: int = 200, json_body=None, text: str = ""):
    resp = MagicMock(spec=requests.Response)
    resp.status_code = status_code
    resp.ok = 200 <= status_code < 300
    resp.text = text
    resp.json.return_value = json_body or {}
    return resp


class TestFetchLigandSMILES(unittest.TestCase):

    def setUp(self):
        # Clear the in-process LRU cache so each test starts clean.
        fetch_ligand_smiles.cache_clear()

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_prefers_openeye_canonical(self, mock_get):
        mock_get.return_value = _make_response(json_body={
            "pdbx_chem_comp_descriptor": [
                {"type": "SMILES_CANONICAL", "program": "OpenEye OEToolkits",
                 "descriptor": "OPENEYE_STEREO"},
                {"type": "SMILES_CANONICAL", "program": "CACTVS",
                 "descriptor": "CACTVS_STEREO"},
                {"type": "SMILES", "program": "ACDLabs",
                 "descriptor": "ACD_FLAT"},
            ],
        })
        self.assertEqual(fetch_ligand_smiles("XYZ"), "OPENEYE_STEREO")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_falls_back_to_cactvs_canonical(self, mock_get):
        mock_get.return_value = _make_response(json_body={
            "pdbx_chem_comp_descriptor": [
                {"type": "SMILES_CANONICAL", "program": "CACTVS",
                 "descriptor": "CACTVS_STEREO"},
                {"type": "SMILES", "program": "ACDLabs",
                 "descriptor": "ACD_FLAT"},
            ],
        })
        self.assertEqual(fetch_ligand_smiles("XYZ"), "CACTVS_STEREO")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_falls_back_to_any_smiles(self, mock_get):
        mock_get.return_value = _make_response(json_body={
            "pdbx_chem_comp_descriptor": [
                {"type": "SMILES", "program": "ACDLabs",
                 "descriptor": "ACD_FLAT"},
            ],
        })
        self.assertEqual(fetch_ligand_smiles("XYZ"), "ACD_FLAT")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_404_raises_RCSBLookupError(self, mock_get):
        mock_get.return_value = _make_response(
            status_code=404, text='{"status":404,"message":"not found"}'
        )
        with self.assertRaises(RCSBLookupError):
            fetch_ligand_smiles("ZZZZ")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_http_error_raises_RCSBLookupError(self, mock_get):
        mock_get.return_value = _make_response(
            status_code=500, text="internal server error"
        )
        with self.assertRaises(RCSBLookupError):
            fetch_ligand_smiles("XYZ")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_network_error_raises_RCSBLookupError(self, mock_get):
        mock_get.side_effect = requests.ConnectionError("boom")
        with self.assertRaises(RCSBLookupError):
            fetch_ligand_smiles("XYZ")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_missing_smiles_raises_RCSBLookupError(self, mock_get):
        mock_get.return_value = _make_response(json_body={
            "pdbx_chem_comp_descriptor": [
                {"type": "InChI", "descriptor": "InChI=..."},
            ],
        })
        with self.assertRaises(RCSBLookupError):
            fetch_ligand_smiles("XYZ")

    def test_empty_id_raises(self):
        with self.assertRaises(RCSBLookupError):
            fetch_ligand_smiles("   ")

    @patch("pestifer.charmmff.ligand_paramgen.rcsb.requests.get")
    def test_case_insensitive_and_strips_whitespace(self, mock_get):
        mock_get.return_value = _make_response(json_body={
            "pdbx_chem_comp_descriptor": [
                {"type": "SMILES_CANONICAL", "program": "OpenEye OEToolkits",
                 "descriptor": "X"},
            ],
        })
        fetch_ligand_smiles("  atp  ")
        # URL should have ATP (uppercase, stripped)
        called_url = mock_get.call_args.args[0]
        self.assertIn("/ATP", called_url)


if __name__ == "__main__":
    unittest.main()
