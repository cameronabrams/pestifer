# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
An alignment coormod: computes the rigid-body transformation that minimizes
RMSD between a selection in the pipeline system and the same (or a different)
selection in a reference coordinate file, then applies that transformation
to the pipeline system.

VMD's ``measure fit`` is used to compute the 4×4 homogeneous transformation
matrix; the reference molecule is loaded and deleted within the generated
script so it leaves no persistent state.
"""
import logging
logger = logging.getLogger(__name__)
from typing import ClassVar
from pydantic import Field, model_validator
from ..core.baseobj import BaseObj, BaseObjList


class Align(BaseObj):
    """
    Align the pipeline system to a reference coordinate file.

    Parameters (from YAML ``specs``)
    ---------------------------------
    ref_pdb : str, optional
        Path or filename of the reference PDB file.  Mutually exclusive with
        ``ref_sourceID``.
    ref_sourceID : str, optional
        RCSB PDB ID for the reference structure.  The PDB file is downloaded
        to the working directory if not already present.  Mutually exclusive
        with ``ref_pdb``.  Exactly one of ``ref_pdb`` / ``ref_sourceID`` must
        be supplied.
    ref_psf : str, optional
        Path to the reference PSF file.  When omitted the reference is loaded
        as a plain PDB (sufficient for coordinate fitting).
    mobile_sel : str, optional
        VMD atomselect string for the pipeline atoms used in the fit.
        Defaults to ``"all"``.
    ref_sel : str, optional
        VMD atomselect string for the reference atoms used in the fit.
        Defaults to ``mobile_sel`` when omitted.
    apply_to : str, optional
        VMD atomselect string for the pipeline atoms that are actually moved
        by the fitted transformation.  Defaults to ``"all"``.

    Notes
    -----
    Both selections passed to ``measure fit`` must contain the same number
    of atoms in the same order.  Mismatched counts will cause a VMD error.
    """

    _required_fields = set()
    _optional_fields = {'ref_pdb', 'ref_sourceID', 'ref_psf', 'mobile_sel', 'ref_sel', 'apply_to'}
    _mutually_exclusive = {frozenset({'ref_pdb', 'ref_sourceID'})}

    ref_pdb: str | None = Field(None, description="Filename of the reference PDB file.")
    ref_sourceID: str | None = Field(None, description="RCSB PDB ID; the file is downloaded to CWD if not already present.")
    ref_psf: str | None = Field(None, description="Path to the reference PSF file (optional).")
    mobile_sel: str = Field('all', description="VMD atomselect string for pipeline atoms used in the fit.")
    ref_sel: str | None = Field(None, description="VMD atomselect string for reference atoms used in the fit (defaults to mobile_sel).")
    apply_to: str = Field('all', description="VMD atomselect string for pipeline atoms to move after fitting.")

    @model_validator(mode='after')
    def _require_one_ref(self):
        if self.ref_pdb is None and self.ref_sourceID is None:
            raise ValueError("Align requires exactly one of 'ref_pdb' or 'ref_sourceID'.")
        return self

    @property
    def effective_ref_pdb(self) -> str:
        """Resolved filename of the reference PDB (local path or derived from ref_sourceID)."""
        return self.ref_pdb if self.ref_pdb is not None else f'{self.ref_sourceID}.pdb'

    _yaml_header: ClassVar[str] = 'align'
    _objcat: ClassVar[str] = 'coord'


class AlignList(BaseObjList[Align]):
    def describe(self):
        return f'<AlignList with {len(self)} alignments>'
