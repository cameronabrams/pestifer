####### TRANSFORM 0 BEGINS #####
# The following mappings of A.U. chains to chains is used:
### BEGIN SEGMENT G ###
set G [atomselect 0 "serial 1 to 1174"]
$G writepdb PROTEIN_G_34_to_185.pdb
set G [atomselect 0 "serial 1175 to 2795"]
$G writepdb PROTEIN_G_187_to_398.pdb
set G [atomselect 0 "serial 2796 to 3543"]
$G writepdb PROTEIN_G_411_to_505.pdb
segment G {
    pdb PROTEIN_G_34_to_185.pdb
    residue 185A GLU G
    residue 185B ASN G
    residue 185C GLN G
    residue 185D GLY G
    residue 185E ASN G
    residue 185F ARG G
    residue 185G SER G
    residue 185H ASN G
    residue 185I ASN G
    residue 185J GLY G
    pdb PROTEIN_G_187_to_398.pdb
    residue 400 THR G
    residue 401 SER G
    residue 402 VAL G
    residue 403 GLN G
    residue 404 GLY G
    residue 405 SER G
    residue 406 ASN G
    residue 407 SER G
    residue 408 THR G
    residue 409 GLY G
    residue 410 SER G
    residue 410A GLY G
    pdb PROTEIN_G_411_to_505.pdb
}
coordpdb PROTEIN_G_34_to_185.pdb G
coordpdb PROTEIN_G_187_to_398.pdb G
coordpdb PROTEIN_G_411_to_505.pdb G
coord G 185A N [cacoIn_nOut 185 G 0]
coord G 400 N [cacoIn_nOut 398 G 0]
### END SEGMENT G ###
### BEGIN SEGMENT B ###
set B [atomselect 0 "serial 3545 to 3722"]
$B writepdb PROTEIN_B_521_to_547.pdb
set B [atomselect 0 "serial 3723 to 4519"]
$B writepdb PROTEIN_B_569_to_664.pdb
segment B {
    pdb PROTEIN_B_521_to_547.pdb
    residue 548 ILE B
    residue 549 VAL B
    residue 550 GLN B
    residue 551 GLN B
    residue 552 GLN B
    residue 553 SER B
    residue 554 ASN B
    residue 555 LEU B
    residue 556 LEU B
    residue 557 ARG B
    residue 558 ALA B
    residue 559 PRO B
    residue 560 GLU B
    residue 561 ALA B
    residue 562 GLN B
    residue 563 GLN B
    residue 564 HSE B
    residue 565 LEU B
    residue 566 LEU B
    residue 567 LYS B
    residue 568 LEU B
    residue 568A GLY B
    pdb PROTEIN_B_569_to_664.pdb
}
coordpdb PROTEIN_B_521_to_547.pdb B
coordpdb PROTEIN_B_569_to_664.pdb B
coord B 548 N [cacoIn_nOut 547 B 0]
### END SEGMENT B ###
####### TRANSFORM 0 ENDS ######
####### TRANSFORM 1 BEGINS #####
# The following mappings of A.U. chains to chains is used:
#   A: E
#   B: F
#   C: H
#   D: I
#   G: J
### BEGIN SEGMENT G ###
set G [atomselect 0 "serial 1 to 1174"]
$G set chain J
set G_orig_x [$G get x]
set G_orig_y [$G get y]
set G_orig_z [$G get z]
set G_orig_resid [$G get resid]
set G_orig_resname [$G get resname]
set G_orig_name [$G get name]
$G move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G writepdb PROTEIN_J_34_to_185.pdb
$G set x $G_orig_x
$G set y $G_orig_y
$G set z $G_orig_z
$G set resid $G_orig_resid
$G set resname $G_orig_resname
$G set name $G_orig_name
set G [atomselect 0 "serial 1175 to 2795"]
$G set chain J
set G_orig_x [$G get x]
set G_orig_y [$G get y]
set G_orig_z [$G get z]
set G_orig_resid [$G get resid]
set G_orig_resname [$G get resname]
set G_orig_name [$G get name]
$G move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G writepdb PROTEIN_J_187_to_398.pdb
$G set x $G_orig_x
$G set y $G_orig_y
$G set z $G_orig_z
$G set resid $G_orig_resid
$G set resname $G_orig_resname
$G set name $G_orig_name
set G [atomselect 0 "serial 2796 to 3543"]
$G set chain J
set G_orig_x [$G get x]
set G_orig_y [$G get y]
set G_orig_z [$G get z]
set G_orig_resid [$G get resid]
set G_orig_resname [$G get resname]
set G_orig_name [$G get name]
$G move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G writepdb PROTEIN_J_411_to_505.pdb
$G set x $G_orig_x
$G set y $G_orig_y
$G set z $G_orig_z
$G set resid $G_orig_resid
$G set resname $G_orig_resname
$G set name $G_orig_name
segment J {
    pdb PROTEIN_J_34_to_185.pdb
    residue 185A GLU G
    residue 185B ASN G
    residue 185C GLN G
    residue 185D GLY G
    residue 185E ASN G
    residue 185F ARG G
    residue 185G SER G
    residue 185H ASN G
    residue 185I ASN G
    residue 185J GLY G
    pdb PROTEIN_J_187_to_398.pdb
    residue 400 THR G
    residue 401 SER G
    residue 402 VAL G
    residue 403 GLN G
    residue 404 GLY G
    residue 405 SER G
    residue 406 ASN G
    residue 407 SER G
    residue 408 THR G
    residue 409 GLY G
    residue 410 SER G
    residue 410A GLY G
    pdb PROTEIN_J_411_to_505.pdb
}
coordpdb PROTEIN_J_34_to_185.pdb J
coordpdb PROTEIN_J_187_to_398.pdb J
coordpdb PROTEIN_J_411_to_505.pdb J
coord J 185A N [cacoIn_nOut 185 J 0]
coord J 400 N [cacoIn_nOut 398 J 0]
### END SEGMENT J ###
### BEGIN SEGMENT B ###
set B [atomselect 0 "serial 3545 to 3722"]
$B set chain F
set B_orig_x [$B get x]
set B_orig_y [$B get y]
set B_orig_z [$B get z]
set B_orig_resid [$B get resid]
set B_orig_resname [$B get resname]
set B_orig_name [$B get name]
$B move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B writepdb PROTEIN_F_521_to_547.pdb
$B set x $B_orig_x
$B set y $B_orig_y
$B set z $B_orig_z
$B set resid $B_orig_resid
$B set resname $B_orig_resname
$B set name $B_orig_name
set B [atomselect 0 "serial 3723 to 4519"]
$B set chain F
set B_orig_x [$B get x]
set B_orig_y [$B get y]
set B_orig_z [$B get z]
set B_orig_resid [$B get resid]
set B_orig_resname [$B get resname]
set B_orig_name [$B get name]
$B move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B writepdb PROTEIN_F_569_to_664.pdb
$B set x $B_orig_x
$B set y $B_orig_y
$B set z $B_orig_z
$B set resid $B_orig_resid
$B set resname $B_orig_resname
$B set name $B_orig_name
segment F {
    pdb PROTEIN_F_521_to_547.pdb
    residue 548 ILE B
    residue 549 VAL B
    residue 550 GLN B
    residue 551 GLN B
    residue 552 GLN B
    residue 553 SER B
    residue 554 ASN B
    residue 555 LEU B
    residue 556 LEU B
    residue 557 ARG B
    residue 558 ALA B
    residue 559 PRO B
    residue 560 GLU B
    residue 561 ALA B
    residue 562 GLN B
    residue 563 GLN B
    residue 564 HSE B
    residue 565 LEU B
    residue 566 LEU B
    residue 567 LYS B
    residue 568 LEU B
    residue 568A GLY B
    pdb PROTEIN_F_569_to_664.pdb
}
coordpdb PROTEIN_F_521_to_547.pdb F
coordpdb PROTEIN_F_569_to_664.pdb F
coord F 548 N [cacoIn_nOut 547 F 0]
### END SEGMENT F ###
####### TRANSFORM 1 ENDS ######
####### TRANSFORM 2 BEGINS #####
# The following mappings of A.U. chains to chains is used:
#   A: K
#   B: L
#   C: M
#   D: N
#   G: O
### BEGIN SEGMENT G ###
set G [atomselect 0 "serial 1 to 1174"]
$G set chain O
set G_orig_x [$G get x]
set G_orig_y [$G get y]
set G_orig_z [$G get z]
set G_orig_resid [$G get resid]
set G_orig_resname [$G get resname]
set G_orig_name [$G get name]
$G move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G writepdb PROTEIN_O_34_to_185.pdb
$G set x $G_orig_x
$G set y $G_orig_y
$G set z $G_orig_z
$G set resid $G_orig_resid
$G set resname $G_orig_resname
$G set name $G_orig_name
set G [atomselect 0 "serial 1175 to 2795"]
$G set chain O
set G_orig_x [$G get x]
set G_orig_y [$G get y]
set G_orig_z [$G get z]
set G_orig_resid [$G get resid]
set G_orig_resname [$G get resname]
set G_orig_name [$G get name]
$G move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G writepdb PROTEIN_O_187_to_398.pdb
$G set x $G_orig_x
$G set y $G_orig_y
$G set z $G_orig_z
$G set resid $G_orig_resid
$G set resname $G_orig_resname
$G set name $G_orig_name
set G [atomselect 0 "serial 2796 to 3543"]
$G set chain O
set G_orig_x [$G get x]
set G_orig_y [$G get y]
set G_orig_z [$G get z]
set G_orig_resid [$G get resid]
set G_orig_resname [$G get resname]
set G_orig_name [$G get name]
$G move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G writepdb PROTEIN_O_411_to_505.pdb
$G set x $G_orig_x
$G set y $G_orig_y
$G set z $G_orig_z
$G set resid $G_orig_resid
$G set resname $G_orig_resname
$G set name $G_orig_name
segment O {
    pdb PROTEIN_O_34_to_185.pdb
    residue 185A GLU G
    residue 185B ASN G
    residue 185C GLN G
    residue 185D GLY G
    residue 185E ASN G
    residue 185F ARG G
    residue 185G SER G
    residue 185H ASN G
    residue 185I ASN G
    residue 185J GLY G
    pdb PROTEIN_O_187_to_398.pdb
    residue 400 THR G
    residue 401 SER G
    residue 402 VAL G
    residue 403 GLN G
    residue 404 GLY G
    residue 405 SER G
    residue 406 ASN G
    residue 407 SER G
    residue 408 THR G
    residue 409 GLY G
    residue 410 SER G
    residue 410A GLY G
    pdb PROTEIN_O_411_to_505.pdb
}
coordpdb PROTEIN_O_34_to_185.pdb O
coordpdb PROTEIN_O_187_to_398.pdb O
coordpdb PROTEIN_O_411_to_505.pdb O
coord O 185A N [cacoIn_nOut 185 O 0]
coord O 400 N [cacoIn_nOut 398 O 0]
### END SEGMENT O ###
### BEGIN SEGMENT B ###
set B [atomselect 0 "serial 3545 to 3722"]
$B set chain L
set B_orig_x [$B get x]
set B_orig_y [$B get y]
set B_orig_z [$B get z]
set B_orig_resid [$B get resid]
set B_orig_resname [$B get resname]
set B_orig_name [$B get name]
$B move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B writepdb PROTEIN_L_521_to_547.pdb
$B set x $B_orig_x
$B set y $B_orig_y
$B set z $B_orig_z
$B set resid $B_orig_resid
$B set resname $B_orig_resname
$B set name $B_orig_name
set B [atomselect 0 "serial 3723 to 4519"]
$B set chain L
set B_orig_x [$B get x]
set B_orig_y [$B get y]
set B_orig_z [$B get z]
set B_orig_resid [$B get resid]
set B_orig_resname [$B get resname]
set B_orig_name [$B get name]
$B move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B writepdb PROTEIN_L_569_to_664.pdb
$B set x $B_orig_x
$B set y $B_orig_y
$B set z $B_orig_z
$B set resid $B_orig_resid
$B set resname $B_orig_resname
$B set name $B_orig_name
segment L {
    pdb PROTEIN_L_521_to_547.pdb
    residue 548 ILE B
    residue 549 VAL B
    residue 550 GLN B
    residue 551 GLN B
    residue 552 GLN B
    residue 553 SER B
    residue 554 ASN B
    residue 555 LEU B
    residue 556 LEU B
    residue 557 ARG B
    residue 558 ALA B
    residue 559 PRO B
    residue 560 GLU B
    residue 561 ALA B
    residue 562 GLN B
    residue 563 GLN B
    residue 564 HSE B
    residue 565 LEU B
    residue 566 LEU B
    residue 567 LYS B
    residue 568 LEU B
    residue 568A GLY B
    pdb PROTEIN_L_569_to_664.pdb
}
coordpdb PROTEIN_L_521_to_547.pdb L
coordpdb PROTEIN_L_569_to_664.pdb L
coord L 548 N [cacoIn_nOut 547 L 0]
### END SEGMENT L ###
####### TRANSFORM 2 ENDS ######