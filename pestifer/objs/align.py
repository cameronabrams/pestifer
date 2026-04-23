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
from pydantic import Field
from ..core.baseobj import BaseObj, BaseObjList


class Align(BaseObj):
    """
    Align the pipeline system to a reference coordinate file.

    Parameters (from YAML ``specs``)
    ---------------------------------
    ref_pdb : str
        Path to the reference PDB file.
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

    _required_fields = {'ref_pdb'}
    _optional_fields = {'ref_psf', 'mobile_sel', 'ref_sel', 'apply_to'}

    ref_pdb: str = Field(..., description="Path to the reference PDB file.")
    ref_psf: str | None = Field(None, description="Path to the reference PSF file (optional).")
    mobile_sel: str = Field('all', description="VMD atomselect string for pipeline atoms used in the fit.")
    ref_sel: str | None = Field(None, description="VMD atomselect string for reference atoms used in the fit (defaults to mobile_sel).")
    apply_to: str = Field('all', description="VMD atomselect string for pipeline atoms to move after fitting.")

    _yaml_header: ClassVar[str] = 'align'
    _objcat: ClassVar[str] = 'coord'


class AlignList(BaseObjList[Align]):
    def describe(self):
        return f'<AlignList with {len(self)} alignments>'
