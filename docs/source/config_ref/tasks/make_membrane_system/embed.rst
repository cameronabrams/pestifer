.. _config_ref tasks make_membrane_system embed:

``embed``
=========

parameters controlling protein embedding

Single-valued parameters:

  * ``q_tolerance``: tolerance for net charge neutralization (default: 0.0001)

  * ``xydist``: distance from perimeter of protein to box edge in x and y (Å) (default: 15)

  * ``zdist``: distance from perimeter of protein to box edge in z (Å) (default: 20)

  * ``margin``: distance from any protein atom in which no lipid atoms are permitted when embedding (Å) (default: 2.4)

  * ``z_head_group``: VMD atomselect string defining head-group of z-axis

  * ``z_tail_group``: VMD atomselect string defining tail-group of z-axis

  * ``z_ref_group``: VMD atomselect string defining a center of mass at a specific z



