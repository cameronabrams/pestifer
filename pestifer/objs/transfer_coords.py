# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A coordinate-transfer coormod: copies atomic coordinates from a selected
subset of a donor system directly onto the corresponding atoms of the
pipeline system.

Unlike :class:`~pestifer.objs.align.Align`, which computes a rigid-body
transformation, ``transfer_coords`` overwrites positions atom-by-atom and
therefore handles non-rigid changes (re-modeled loops, isolated
energy-minimization, ligand repositioning, etc.).

An optional pre-alignment step can be requested by supplying both
``align_donor_sel`` and ``align_mobile_sel``.  When present, the entire
donor molecule is rigidly fitted to the pipeline system using those
selections before the coordinate transfer is performed.  Supplying only
one of the two raises a ``ValueError``.

All selection strings must be supplied explicitly; no defaults are inferred
from each other in order to avoid hidden decision-making.  Congruent
selections (same atom count, same order) are required for both the
alignment pair and the transfer pair.
"""
import logging
logger = logging.getLogger(__name__)
from typing import ClassVar
from pydantic import Field, model_validator
from ..core.baseobj import BaseObj, BaseObjList


class TransferCoords(BaseObj):
    """
    Copy coordinates from a donor PDB onto a selection in the pipeline system,
    with an optional prior rigid-body alignment of the donor to the pipeline.

    Parameters (from YAML ``specs``)
    ---------------------------------
    donor_pdb : str
        Path to the donor PDB file.
    donor_psf : str, optional
        Path to a donor PSF file.  When provided both files are loaded so
        that VMD selection syntax can reference topology attributes.
    donor_sel : str
        VMD atomselect string for the atoms to copy **from** in the donor.
    mobile_sel : str
        VMD atomselect string for the atoms to overwrite **in** the pipeline
        system.  Must select the same number of atoms as ``donor_sel``,
        in the same order.
    align_donor_sel : str, optional
        VMD atomselect string for donor atoms used in the pre-alignment fit.
        Must be supplied together with ``align_mobile_sel``; omitting both
        skips the alignment step entirely.
    align_mobile_sel : str, optional
        VMD atomselect string for pipeline atoms used in the pre-alignment
        fit.  Must be supplied together with ``align_donor_sel``.
    """

    _required_fields = {'donor_pdb', 'donor_sel', 'mobile_sel'}
    _optional_fields = {'donor_psf', 'align_donor_sel', 'align_mobile_sel'}

    donor_pdb: str = Field(..., description="Path to the donor PDB file.")
    donor_psf: str | None = Field(None, description="Path to the donor PSF file (optional).")
    donor_sel: str = Field(..., description="VMD atomselect string for atoms to copy from in the donor.")
    mobile_sel: str = Field(..., description="VMD atomselect string for atoms to overwrite in the pipeline system.")
    align_donor_sel: str | None = Field(None, description="Donor selection for pre-alignment fit (requires align_mobile_sel).")
    align_mobile_sel: str | None = Field(None, description="Pipeline selection for pre-alignment fit (requires align_donor_sel).")

    _yaml_header: ClassVar[str] = 'transfer_coords'
    _objcat: ClassVar[str] = 'coord'

    @model_validator(mode='after')
    def _check_align_pair(self):
        has_donor = self.align_donor_sel is not None
        has_mobile = self.align_mobile_sel is not None
        if has_donor != has_mobile:
            missing = 'align_mobile_sel' if has_donor else 'align_donor_sel'
            present = 'align_donor_sel' if has_donor else 'align_mobile_sel'
            raise ValueError(
                f"transfer_coords: '{present}' was supplied but '{missing}' is missing. "
                f"Both alignment selections must be provided together, or neither."
            )
        return self

    @property
    def pre_align(self) -> bool:
        return self.align_donor_sel is not None


class TransferCoordsList(BaseObjList[TransferCoords]):
    def describe(self):
        return f'<TransferCoordsList with {len(self)} transfers>'
