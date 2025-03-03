.. _config_ref tasks bilayer embed:

``embed``
=========

parameters controlling protein embedding

Single-valued parameters:

  * ``xydist``: distance from perimeter of protein to box edge in x and y (A) (default: 15)

  * ``zdist``: distance from perimeter of protein to box edge in z (A) (default: 20)

  * ``protein_radius_scaling``: factor by which to multiply computed protein TM radius to compute membrane excluded area (default: 1.01)

  * ``n_ter``: in or out; N-termini are usually 'in' for membrane proteins (default: in)

  * ``z_head_group``: VMD atomselect string defining head-group of z-axis

  * ``z_tail_group``: VMD atomselect string defining tail-group of z-axis

  * ``z_ref_group``: VMD atomselect string defining a center of mass at a specific z



