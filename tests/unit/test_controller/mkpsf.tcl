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
################ The following mappings of A.U. chains is used: ################
# A: A
# B: B
# C: C
# D: D
# G: G
############################### SEGMENTS FOLLOW ################################
############################### BEGIN SEGMENT G ################################
set G01 [atomselect 0 "serial 1 to 1174"]
############ Atom with serial 1174 in PDB needs serial 1174 for VMD ############
$G01 writepdb PROTEIN_G_34_to_185.pdb
set G03 [atomselect 0 "serial 1175 to 2795"]
############ Atom with serial 2795 in PDB needs serial 2795 for VMD ############
$G03 writepdb PROTEIN_G_187_to_398.pdb
set G05 [atomselect 0 "serial 2796 to 3543"]
############ Atom with serial 3543 in PDB needs serial 3543 for VMD ############
$G05 writepdb PROTEIN_G_411_to_505.pdb
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
    mutate 332 THR
    mutate 501 ALA
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_G_34_to_185.pdb G
coordpdb PROTEIN_G_187_to_398.pdb G
coordpdb PROTEIN_G_411_to_505.pdb G
######## Seeding orientation of model-built loop starting at G-GLU185A #########
coord G 185A N [cacoIn_nOut 185 G 0]
######### Seeding orientation of model-built loop starting at G-THR400 #########
coord G 400 N [cacoIn_nOut 398 G 0]
####################### Intra-segmental terminal patches #######################
patch CTER G:185I
patch NTER G:187
delatom G 185J
patch CTER G:410
patch NTER G:411
delatom G 410A
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT G #################################
############################### BEGIN SEGMENT B ################################
set B01 [atomselect 0 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
$B01 writepdb PROTEIN_B_521_to_547.pdb
set B03 [atomselect 0 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
$B03 writepdb PROTEIN_B_569_to_664.pdb
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
    mutate 559 ILE
    mutate 605 THR
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_B_521_to_547.pdb B
coordpdb PROTEIN_B_569_to_664.pdb B
######### Seeding orientation of model-built loop starting at B-ILE548 #########
coord B 548 N [cacoIn_nOut 547 B 0]
####################### Intra-segmental terminal patches #######################
patch CTER B:568
patch NTER B:569
delatom B 568A
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT B #################################
############################# DISU PATCHES FOLLOW ##############################
patch DISU G:54 G:74
patch DISU G:119 G:205
patch DISU G:126 G:196
patch DISU G:131 G:157
patch DISU G:218 G:247
patch DISU G:228 G:239
patch DISU G:296 G:331
patch DISU G:378 G:445
patch DISU G:385 G:418
patch DISU B:598 B:604
############################# LINK PATCHES FOLLOW ##############################
############################### TRANSFORM 0 ENDS ###############################
############################## TRANSFORM 1 BEGINS ##############################
################ The following mappings of A.U. chains is used: ################
# A: E
# B: F
# C: H
# D: I
# G: J
############################### SEGMENTS FOLLOW ################################
############################### BEGIN SEGMENT J ################################
set G00 [atomselect 0 "serial 1 to 1174"]
############ Atom with serial 1174 in PDB needs serial 1174 for VMD ############
set G00_orig_chain [$G00 get chain]
set G00_orig_x [$G00 get x]
set G00_orig_y [$G00 get y]
set G00_orig_z [$G00 get z]
set G00_orig_resid [$G00 get resid]
set G00_orig_resname [$G00 get resname]
set G00_orig_name [$G00 get name]
$G00 set chain J
$G00 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G00 writepdb PROTEIN_J_34_to_185.pdb
set G02 [atomselect 0 "serial 1175 to 2795"]
############ Atom with serial 2795 in PDB needs serial 2795 for VMD ############
set G02_orig_chain [$G02 get chain]
set G02_orig_x [$G02 get x]
set G02_orig_y [$G02 get y]
set G02_orig_z [$G02 get z]
set G02_orig_resid [$G02 get resid]
set G02_orig_resname [$G02 get resname]
set G02_orig_name [$G02 get name]
$G02 set chain J
$G02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G02 writepdb PROTEIN_J_187_to_398.pdb
set G04 [atomselect 0 "serial 2796 to 3543"]
############ Atom with serial 3543 in PDB needs serial 3543 for VMD ############
set G04_orig_chain [$G04 get chain]
set G04_orig_x [$G04 get x]
set G04_orig_y [$G04 get y]
set G04_orig_z [$G04 get z]
set G04_orig_resid [$G04 get resid]
set G04_orig_resname [$G04 get resname]
set G04_orig_name [$G04 get name]
$G04 set chain J
$G04 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G04 writepdb PROTEIN_J_411_to_505.pdb
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
    residue 185J GLY J
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
    residue 410A GLY J
    pdb PROTEIN_J_411_to_505.pdb
    mutate 332 THR
    mutate 501 ALA
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_J_34_to_185.pdb J
coordpdb PROTEIN_J_187_to_398.pdb J
coordpdb PROTEIN_J_411_to_505.pdb J
######## Seeding orientation of model-built loop starting at G-GLU185A #########
coord J 185A N [cacoIn_nOut 185 J 0]
######### Seeding orientation of model-built loop starting at G-THR400 #########
coord J 400 N [cacoIn_nOut 398 J 0]
####################### Intra-segmental terminal patches #######################
patch CTER J:185I
patch NTER J:187
delatom J 185J
patch CTER J:410
patch NTER J:411
delatom J 410A
############## Restoring A.U. state for all resolved subsegments ###############
$G00 set chain $G00_orig_chain
$G00 set x $G00_orig_x
$G00 set y $G00_orig_y
$G00 set z $G00_orig_z
$G00 set resid $G00_orig_resid
$G00 set resname $G00_orig_resname
$G00 set name $G00_orig_name
$G02 set chain $G02_orig_chain
$G02 set x $G02_orig_x
$G02 set y $G02_orig_y
$G02 set z $G02_orig_z
$G02 set resid $G02_orig_resid
$G02 set resname $G02_orig_resname
$G02 set name $G02_orig_name
$G04 set chain $G04_orig_chain
$G04 set x $G04_orig_x
$G04 set y $G04_orig_y
$G04 set z $G04_orig_z
$G04 set resid $G04_orig_resid
$G04 set resname $G04_orig_resname
$G04 set name $G04_orig_name
################################ END SEGMENT J #################################
############################### BEGIN SEGMENT F ################################
set B00 [atomselect 0 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
set B00_orig_chain [$B00 get chain]
set B00_orig_x [$B00 get x]
set B00_orig_y [$B00 get y]
set B00_orig_z [$B00 get z]
set B00_orig_resid [$B00 get resid]
set B00_orig_resname [$B00 get resname]
set B00_orig_name [$B00 get name]
$B00 set chain F
$B00 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B00 writepdb PROTEIN_F_521_to_547.pdb
set B02 [atomselect 0 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
set B02_orig_chain [$B02 get chain]
set B02_orig_x [$B02 get x]
set B02_orig_y [$B02 get y]
set B02_orig_z [$B02 get z]
set B02_orig_resid [$B02 get resid]
set B02_orig_resname [$B02 get resname]
set B02_orig_name [$B02 get name]
$B02 set chain F
$B02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B02 writepdb PROTEIN_F_569_to_664.pdb
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
    residue 568A GLY F
    pdb PROTEIN_F_569_to_664.pdb
    mutate 559 ILE
    mutate 605 THR
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_F_521_to_547.pdb F
coordpdb PROTEIN_F_569_to_664.pdb F
######### Seeding orientation of model-built loop starting at B-ILE548 #########
coord F 548 N [cacoIn_nOut 547 F 0]
####################### Intra-segmental terminal patches #######################
patch CTER F:568
patch NTER F:569
delatom F 568A
############## Restoring A.U. state for all resolved subsegments ###############
$B00 set chain $B00_orig_chain
$B00 set x $B00_orig_x
$B00 set y $B00_orig_y
$B00 set z $B00_orig_z
$B00 set resid $B00_orig_resid
$B00 set resname $B00_orig_resname
$B00 set name $B00_orig_name
$B02 set chain $B02_orig_chain
$B02 set x $B02_orig_x
$B02 set y $B02_orig_y
$B02 set z $B02_orig_z
$B02 set resid $B02_orig_resid
$B02 set resname $B02_orig_resname
$B02 set name $B02_orig_name
################################ END SEGMENT F #################################
############################# DISU PATCHES FOLLOW ##############################
patch DISU J:54 J:74
patch DISU J:119 J:205
patch DISU J:126 J:196
patch DISU J:131 J:157
patch DISU J:218 J:247
patch DISU J:228 J:239
patch DISU J:296 J:331
patch DISU J:378 J:445
patch DISU J:385 J:418
patch DISU F:598 F:604
############################# LINK PATCHES FOLLOW ##############################
############################### TRANSFORM 1 ENDS ###############################
############################## TRANSFORM 2 BEGINS ##############################
################ The following mappings of A.U. chains is used: ################
# A: K
# B: L
# C: M
# D: N
# G: O
############################### SEGMENTS FOLLOW ################################
############################### BEGIN SEGMENT O ################################
set G00 [atomselect 0 "serial 1 to 1174"]
############ Atom with serial 1174 in PDB needs serial 1174 for VMD ############
set G00_orig_chain [$G00 get chain]
set G00_orig_x [$G00 get x]
set G00_orig_y [$G00 get y]
set G00_orig_z [$G00 get z]
set G00_orig_resid [$G00 get resid]
set G00_orig_resname [$G00 get resname]
set G00_orig_name [$G00 get name]
$G00 set chain O
$G00 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G00 writepdb PROTEIN_O_34_to_185.pdb
set G02 [atomselect 0 "serial 1175 to 2795"]
############ Atom with serial 2795 in PDB needs serial 2795 for VMD ############
set G02_orig_chain [$G02 get chain]
set G02_orig_x [$G02 get x]
set G02_orig_y [$G02 get y]
set G02_orig_z [$G02 get z]
set G02_orig_resid [$G02 get resid]
set G02_orig_resname [$G02 get resname]
set G02_orig_name [$G02 get name]
$G02 set chain O
$G02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G02 writepdb PROTEIN_O_187_to_398.pdb
set G04 [atomselect 0 "serial 2796 to 3543"]
############ Atom with serial 3543 in PDB needs serial 3543 for VMD ############
set G04_orig_chain [$G04 get chain]
set G04_orig_x [$G04 get x]
set G04_orig_y [$G04 get y]
set G04_orig_z [$G04 get z]
set G04_orig_resid [$G04 get resid]
set G04_orig_resname [$G04 get resname]
set G04_orig_name [$G04 get name]
$G04 set chain O
$G04 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$G04 writepdb PROTEIN_O_411_to_505.pdb
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
    residue 185J GLY O
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
    residue 410A GLY O
    pdb PROTEIN_O_411_to_505.pdb
    mutate 332 THR
    mutate 501 ALA
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_O_34_to_185.pdb O
coordpdb PROTEIN_O_187_to_398.pdb O
coordpdb PROTEIN_O_411_to_505.pdb O
######## Seeding orientation of model-built loop starting at G-GLU185A #########
coord O 185A N [cacoIn_nOut 185 O 0]
######### Seeding orientation of model-built loop starting at G-THR400 #########
coord O 400 N [cacoIn_nOut 398 O 0]
####################### Intra-segmental terminal patches #######################
patch CTER O:185I
patch NTER O:187
delatom O 185J
patch CTER O:410
patch NTER O:411
delatom O 410A
############## Restoring A.U. state for all resolved subsegments ###############
$G00 set chain $G00_orig_chain
$G00 set x $G00_orig_x
$G00 set y $G00_orig_y
$G00 set z $G00_orig_z
$G00 set resid $G00_orig_resid
$G00 set resname $G00_orig_resname
$G00 set name $G00_orig_name
$G02 set chain $G02_orig_chain
$G02 set x $G02_orig_x
$G02 set y $G02_orig_y
$G02 set z $G02_orig_z
$G02 set resid $G02_orig_resid
$G02 set resname $G02_orig_resname
$G02 set name $G02_orig_name
$G04 set chain $G04_orig_chain
$G04 set x $G04_orig_x
$G04 set y $G04_orig_y
$G04 set z $G04_orig_z
$G04 set resid $G04_orig_resid
$G04 set resname $G04_orig_resname
$G04 set name $G04_orig_name
################################ END SEGMENT O #################################
############################### BEGIN SEGMENT L ################################
set B00 [atomselect 0 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
set B00_orig_chain [$B00 get chain]
set B00_orig_x [$B00 get x]
set B00_orig_y [$B00 get y]
set B00_orig_z [$B00 get z]
set B00_orig_resid [$B00 get resid]
set B00_orig_resname [$B00 get resname]
set B00_orig_name [$B00 get name]
$B00 set chain L
$B00 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B00 writepdb PROTEIN_L_521_to_547.pdb
set B02 [atomselect 0 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
set B02_orig_chain [$B02 get chain]
set B02_orig_x [$B02 get x]
set B02_orig_y [$B02 get y]
set B02_orig_z [$B02 get z]
set B02_orig_resid [$B02 get resid]
set B02_orig_resname [$B02 get resname]
set B02_orig_name [$B02 get name]
$B02 set chain L
$B02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$B02 writepdb PROTEIN_L_569_to_664.pdb
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
    residue 568A GLY L
    pdb PROTEIN_L_569_to_664.pdb
    mutate 559 ILE
    mutate 605 THR
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_L_521_to_547.pdb L
coordpdb PROTEIN_L_569_to_664.pdb L
######### Seeding orientation of model-built loop starting at B-ILE548 #########
coord L 548 N [cacoIn_nOut 547 L 0]
####################### Intra-segmental terminal patches #######################
patch CTER L:568
patch NTER L:569
delatom L 568A
############## Restoring A.U. state for all resolved subsegments ###############
$B00 set chain $B00_orig_chain
$B00 set x $B00_orig_x
$B00 set y $B00_orig_y
$B00 set z $B00_orig_z
$B00 set resid $B00_orig_resid
$B00 set resname $B00_orig_resname
$B00 set name $B00_orig_name
$B02 set chain $B02_orig_chain
$B02 set x $B02_orig_x
$B02 set y $B02_orig_y
$B02 set z $B02_orig_z
$B02 set resid $B02_orig_resid
$B02 set resname $B02_orig_resname
$B02 set name $B02_orig_name
################################ END SEGMENT L #################################
############################# DISU PATCHES FOLLOW ##############################
patch DISU O:54 O:74
patch DISU O:119 O:205
patch DISU O:126 O:196
patch DISU O:131 O:157
patch DISU O:218 O:247
patch DISU O:228 O:239
patch DISU O:296 O:331
patch DISU O:378 O:445
patch DISU O:385 O:418
patch DISU L:598 L:604
############################# LINK PATCHES FOLLOW ##############################
############################### TRANSFORM 2 ENDS ###############################

guesscoord
regenerate angles dihedrals
writepsf cmap step1.psf
writepdb step1.pdb
exit
############################# END PESTIFER PSFGEN ##############################
######################## Thank you for using pestifer! #########################
