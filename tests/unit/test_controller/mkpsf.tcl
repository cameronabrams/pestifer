############################ BEGIN PESTIFER PSFGEN #############################
################################# BEGIN HEADER #################################
source /home/cfa/Git/pestifer/pestifer/Resources/tcl/modules/src/loopmc.tcl
source /home/cfa/Git/pestifer/pestifer/Resources/tcl/vmdrc.tcl
package require psfgen
psfcontext mixedcase
topology /home/cfa/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all35_ethers.rtf
topology /home/cfa/charmm/toppar/top_all36_cgenff.rtf
topology /home/cfa/charmm/toppar/top_all36_lipid.rtf
topology /home/cfa/charmm/toppar/top_all36_na.rtf
topology /home/cfa/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology /home/cfa/Git/pestifer/pestifer/Resources/charmm/top_all36_carb.rtf
topology /home/cfa/Git/pestifer/pestifer/Resources/charmm/toppar_water_ions.str
pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias residue NAG BGNA
pdbalias atom BGNA C7 C
pdbalias atom BGNA O7 O
pdbalias atom BGNA C8 CT
pdbalias atom BGNA N2 N
pdbalias residue SIA ANE5
pdbalias atom ANE5 C10 C
pdbalias atom ANE5 C11 CT
pdbalias atom ANE5 N5 N
pdbalias atom ANE5 O1A O11
pdbalias atom ANE5 O1B O12
pdbalias atom ANE5 O10 O
pdbalias atom VCG C01 C1
pdbalias atom VCG C01 C1
pdbalias atom VCG C02 C2
pdbalias atom VCG C03 C3
pdbalias atom VCG C04 C4
pdbalias atom VCG C05 C5
pdbalias atom VCG C06 C6
pdbalias atom VCG C07 C7
pdbalias atom VCG C08 C8
pdbalias atom VCG C09 C9
pdbalias residue EIC LIN
set RESDICT(HIS) HSE
set RESDICT(ZN) ZN2
set RESDICT(HOH) TIP3
set RESDICT(CL) CLA
set RESDICT(NAG) BGNA
set RESDICT(MAN) AMAN
set RESDICT(BMA) BMAN
set RESDICT(FUC) AFUC
set RESDICT(GAL) BGAL
set RESDICT(SIA) ANE5AC
set RESDICT(ANE5) ANE5AC
set RESDICT(EIC) LIN
set RESDICT(VCG) VCG
set ANAMEDICT(CL) CLA
set ANAMEDICT(O) OH2
################################## END HEADER ##################################
mol new 4zmj waitfor all
set m0 [molinfo top get id]
mol top $m0
############################## TRANSFORM 0 BEGINS ##############################
############################### BEGIN SEGMENT G ################################
set G [atomselect 0 "serial 1 to 1174"]
############ Atom with serial 1174 in PDB needs serial 1174 for VMD ############
$G writepdb PROTEIN_G_34_to_185.pdb
set G [atomselect 0 "serial 1175 to 2795"]
############ Atom with serial 2795 in PDB needs serial 2795 for VMD ############
$G writepdb PROTEIN_G_187_to_398.pdb
set G [atomselect 0 "serial 2796 to 3543"]
############ Atom with serial 3543 in PDB needs serial 3543 for VMD ############
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
################################ END SEGMENT G #################################
############################### BEGIN SEGMENT B ################################
set B [atomselect 0 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
$B writepdb PROTEIN_B_521_to_547.pdb
set B [atomselect 0 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
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
################################ END SEGMENT B #################################
patch DISU G:54 G:74
patch DISU G:119 G:205
patch DISU G:126 G:196
patch DISU G:131 G:157
patch DISU G:218 G:247
patch DISU G:228 G:239
patch DISU G:296 G:331
patch DISU G:378 G:445
patch DISU G:385 G:418
patch DISU G:501 B:605
patch DISU B:598 B:604
############################### TRANSFORM 0 ENDS ###############################
############################## TRANSFORM 1 BEGINS ##############################
################ The following mappings of A.U. chains is used: ################
#   A: E
#   B: F
#   C: H
#   D: I
#   G: J
############################### BEGIN SEGMENT J ################################
set G [atomselect 0 "serial 1 to 1174"]
############ Atom with serial 1174 in PDB needs serial 1174 for VMD ############
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
############ Atom with serial 2795 in PDB needs serial 2795 for VMD ############
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
############ Atom with serial 3543 in PDB needs serial 3543 for VMD ############
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
    residue 185A GLU J
    residue 185B ASN J
    residue 185C GLN J
    residue 185D GLY J
    residue 185E ASN J
    residue 185F ARG J
    residue 185G SER J
    residue 185H ASN J
    residue 185I ASN J
    residue 185J GLY G
    pdb PROTEIN_J_187_to_398.pdb
    residue 400 THR J
    residue 401 SER J
    residue 402 VAL J
    residue 403 GLN J
    residue 404 GLY J
    residue 405 SER J
    residue 406 ASN J
    residue 407 SER J
    residue 408 THR J
    residue 409 GLY J
    residue 410 SER J
    residue 410A GLY G
    pdb PROTEIN_J_411_to_505.pdb
}
coordpdb PROTEIN_J_34_to_185.pdb J
coordpdb PROTEIN_J_187_to_398.pdb J
coordpdb PROTEIN_J_411_to_505.pdb J
coord J 185A N [cacoIn_nOut 185 J 0]
coord J 400 N [cacoIn_nOut 398 J 0]
################################ END SEGMENT J #################################
############################### BEGIN SEGMENT F ################################
set B [atomselect 0 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
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
set B [atomselect 0 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
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
    residue 548 ILE F
    residue 549 VAL F
    residue 550 GLN F
    residue 551 GLN F
    residue 552 GLN F
    residue 553 SER F
    residue 554 ASN F
    residue 555 LEU F
    residue 556 LEU F
    residue 557 ARG F
    residue 558 ALA F
    residue 559 PRO F
    residue 560 GLU F
    residue 561 ALA F
    residue 562 GLN F
    residue 563 GLN F
    residue 564 HSE F
    residue 565 LEU F
    residue 566 LEU F
    residue 567 LYS F
    residue 568 LEU F
    residue 568A GLY B
    pdb PROTEIN_F_569_to_664.pdb
}
coordpdb PROTEIN_F_521_to_547.pdb F
coordpdb PROTEIN_F_569_to_664.pdb F
coord F 548 N [cacoIn_nOut 547 F 0]
################################ END SEGMENT F #################################
patch DISU J:54 J:74
patch DISU J:119 J:205
patch DISU J:126 J:196
patch DISU J:131 J:157
patch DISU J:218 J:247
patch DISU J:228 J:239
patch DISU J:296 J:331
patch DISU J:378 J:445
patch DISU J:385 J:418
patch DISU J:501 F:605
patch DISU F:598 F:604
############################### TRANSFORM 1 ENDS ###############################
############################## TRANSFORM 2 BEGINS ##############################
################ The following mappings of A.U. chains is used: ################
#   A: K
#   B: L
#   C: M
#   D: N
#   G: O
############################### BEGIN SEGMENT O ################################
set G [atomselect 0 "serial 1 to 1174"]
############ Atom with serial 1174 in PDB needs serial 1174 for VMD ############
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
############ Atom with serial 2795 in PDB needs serial 2795 for VMD ############
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
############ Atom with serial 3543 in PDB needs serial 3543 for VMD ############
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
    residue 185A GLU O
    residue 185B ASN O
    residue 185C GLN O
    residue 185D GLY O
    residue 185E ASN O
    residue 185F ARG O
    residue 185G SER O
    residue 185H ASN O
    residue 185I ASN O
    residue 185J GLY G
    pdb PROTEIN_O_187_to_398.pdb
    residue 400 THR O
    residue 401 SER O
    residue 402 VAL O
    residue 403 GLN O
    residue 404 GLY O
    residue 405 SER O
    residue 406 ASN O
    residue 407 SER O
    residue 408 THR O
    residue 409 GLY O
    residue 410 SER O
    residue 410A GLY G
    pdb PROTEIN_O_411_to_505.pdb
}
coordpdb PROTEIN_O_34_to_185.pdb O
coordpdb PROTEIN_O_187_to_398.pdb O
coordpdb PROTEIN_O_411_to_505.pdb O
coord O 185A N [cacoIn_nOut 185 O 0]
coord O 400 N [cacoIn_nOut 398 O 0]
################################ END SEGMENT O #################################
############################### BEGIN SEGMENT L ################################
set B [atomselect 0 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
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
set B [atomselect 0 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
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
    residue 548 ILE L
    residue 549 VAL L
    residue 550 GLN L
    residue 551 GLN L
    residue 552 GLN L
    residue 553 SER L
    residue 554 ASN L
    residue 555 LEU L
    residue 556 LEU L
    residue 557 ARG L
    residue 558 ALA L
    residue 559 PRO L
    residue 560 GLU L
    residue 561 ALA L
    residue 562 GLN L
    residue 563 GLN L
    residue 564 HSE L
    residue 565 LEU L
    residue 566 LEU L
    residue 567 LYS L
    residue 568 LEU L
    residue 568A GLY B
    pdb PROTEIN_L_569_to_664.pdb
}
coordpdb PROTEIN_L_521_to_547.pdb L
coordpdb PROTEIN_L_569_to_664.pdb L
coord L 548 N [cacoIn_nOut 547 L 0]
################################ END SEGMENT L #################################
patch DISU O:54 O:74
patch DISU O:119 O:205
patch DISU O:126 O:196
patch DISU O:131 O:157
patch DISU O:218 O:247
patch DISU O:228 O:239
patch DISU O:296 O:331
patch DISU O:378 O:445
patch DISU O:385 O:418
patch DISU O:501 L:605
patch DISU L:598 L:604
############################### TRANSFORM 2 ENDS ###############################

guesscoord
regenerate angles dihedrals
writepsf cmap step1.psf
writepdb step1.pdb
exit
############################# END PESTIFER PSFGEN ##############################
######################## Thank you for using pestifer! #########################
