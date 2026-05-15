# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Look up SMILES for a PDB chemical-component ID via the RCSB Data API.

The RCSB Data API exposes per-component records at
``https://data.rcsb.org/rest/v1/core/chemcomp/{COMP_ID}``. Each record's
``pdbx_chem_comp_descriptor`` array carries several SMILES variants
(OpenEye/CACTVS/ACDLabs, with and without stereo); we prefer the
OpenEye canonical-stereo form when present.

Used by the CGenFF ligand-parameterization pipeline so a HETATM resname
can be turned into a usable SMILES without the user having to supply
one for any ligand the PDB already knows about.
"""
from __future__ import annotations

import logging
from functools import lru_cache

import requests

logger = logging.getLogger(__name__)

RCSB_CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"

# Ordered preference: (program-substring, type). Earlier entries win.
# OpenEye canonical SMILES carry stereo; CACTVS canonical do too. The
# bare "SMILES" entries are non-canonical fallbacks. ACDLabs entries are
# last-ditch fallbacks for very old components that lack the others.
_SMILES_PREFERENCE: tuple[tuple[str, str], ...] = (
    ("OpenEye", "SMILES_CANONICAL"),
    ("CACTVS",  "SMILES_CANONICAL"),
    ("OpenEye", "SMILES"),
    ("CACTVS",  "SMILES"),
    ("",        "SMILES_CANONICAL"),
    ("",        "SMILES"),
)


class RCSBLookupError(RuntimeError):
    """Raised when an RCSB SMILES lookup fails."""


@lru_cache(maxsize=512)
def fetch_ligand_smiles(comp_id: str, *, timeout: float = 10.0) -> str:
    """
    Fetch the best-available SMILES for a PDB chemical component.

    Parameters
    ----------
    comp_id
        PDB chemical-component ID (e.g. ``"ATP"``, ``"NAG"``). Case-insensitive;
        whitespace stripped.
    timeout
        HTTP request timeout in seconds.

    Returns
    -------
    str
        SMILES string. Preference order: OpenEye canonical (stereo) > CACTVS
        canonical (stereo) > OpenEye non-canonical > CACTVS non-canonical >
        any remaining ``SMILES_CANONICAL`` > any remaining ``SMILES``.

    Raises
    ------
    RCSBLookupError
        If the component ID is not found, the network call fails, or the
        response carries no SMILES descriptor.
    """
    key = comp_id.strip().upper()
    if not key:
        raise RCSBLookupError("Empty PDB component ID.")
    url = RCSB_CHEMCOMP_URL.format(comp_id=key)
    try:
        resp = requests.get(url, timeout=timeout)
    except requests.RequestException as exc:
        raise RCSBLookupError(f"RCSB request failed for {key!r}: {exc}") from exc
    if resp.status_code == 404:
        raise RCSBLookupError(f"No PDB chemical component found for {key!r}.")
    if not resp.ok:
        raise RCSBLookupError(
            f"RCSB returned HTTP {resp.status_code} for {key!r}: "
            f"{resp.text[:200]!r}"
        )
    descriptors = resp.json().get("pdbx_chem_comp_descriptor") or []
    smiles = _pick_best_smiles(descriptors)
    if smiles is None:
        raise RCSBLookupError(f"No SMILES descriptor in RCSB response for {key!r}.")
    logger.debug("RCSB SMILES for %s: %s", key, smiles)
    return smiles


def _pick_best_smiles(descriptors: list[dict]) -> str | None:
    """Pick the highest-priority SMILES from a pdbx_chem_comp_descriptor list."""
    for prog_substr, type_match in _SMILES_PREFERENCE:
        for d in descriptors:
            if d.get("type") != type_match:
                continue
            if prog_substr and prog_substr not in (d.get("program") or ""):
                continue
            descriptor = d.get("descriptor")
            if descriptor:
                return descriptor
    return None
