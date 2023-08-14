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
mol new 4zmj.pdb waitfor all
set m1 [molinfo top get id]
mol top $m1
############################## TRANSFORM 0 BEGINS ##############################
################ The following mappings of A.U. chains is used: ################
# A: A
# B: B
# C: C
# D: D
# G: G
############################### SEGMENTS FOLLOW ################################
# Writing seg G
############################### BEGIN SEGMENT G ################################
set G01 [atomselect $m1 "serial 1 to 1174"]
$G01 writepdb PROTEIN_G_34_to_185.pdb
set G03 [atomselect $m1 "serial 1175 to 2795"]
$G03 writepdb PROTEIN_G_187_to_398.pdb
set G05 [atomselect $m1 "serial 2796 to 3543"]
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
# Below reverts an engineered mutation:
    mutate 332 THR
# Below reverts an engineered mutation:
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
# Done with seg G
# Writing seg B
############################### BEGIN SEGMENT B ################################
set B01 [atomselect $m1 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
$B01 writepdb PROTEIN_B_521_to_547.pdb
set B03 [atomselect $m1 "serial 3722 to 4518"]
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
# Below reverts an engineered mutation:
    mutate 559 ILE
# Below reverts an engineered mutation:
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
# Done with seg B
# Writing seg AG
###################### Not generating segment AG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg AG
# Writing seg CG
###################### Not generating segment CG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg CG
# Writing seg DG
###################### Not generating segment DG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg DG
# Writing seg GG
###################### Not generating segment GG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg GG
# Writing seg BG
###################### Not generating segment BG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg BG
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
# Writing seg G
############################### BEGIN SEGMENT J ################################
set J00 [atomselect $m1 "serial 1 to 1174"]
set J00_orig_chain [$J00 get chain]
set J00_orig_x [$J00 get x]
set J00_orig_y [$J00 get y]
set J00_orig_z [$J00 get z]
set J00_orig_resid [$J00 get resid]
set J00_orig_resname [$J00 get resname]
set J00_orig_name [$J00 get name]
$J00 set chain J
$J00 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$J00 writepdb PROTEIN_J_34_to_185.pdb
set J02 [atomselect $m1 "serial 1175 to 2795"]
set J02_orig_chain [$J02 get chain]
set J02_orig_x [$J02 get x]
set J02_orig_y [$J02 get y]
set J02_orig_z [$J02 get z]
set J02_orig_resid [$J02 get resid]
set J02_orig_resname [$J02 get resname]
set J02_orig_name [$J02 get name]
$J02 set chain J
$J02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$J02 writepdb PROTEIN_J_187_to_398.pdb
set J04 [atomselect $m1 "serial 2796 to 3543"]
set J04_orig_chain [$J04 get chain]
set J04_orig_x [$J04 get x]
set J04_orig_y [$J04 get y]
set J04_orig_z [$J04 get z]
set J04_orig_resid [$J04 get resid]
set J04_orig_resname [$J04 get resname]
set J04_orig_name [$J04 get name]
$J04 set chain J
$J04 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$J04 writepdb PROTEIN_J_411_to_505.pdb
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
# Below reverts an engineered mutation:
    mutate 332 THR
# Below reverts an engineered mutation:
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
$J00 set chain $J00_orig_chain
$J00 set x $J00_orig_x
$J00 set y $J00_orig_y
$J00 set z $J00_orig_z
$J00 set resid $J00_orig_resid
$J00 set resname $J00_orig_resname
$J00 set name $J00_orig_name
$J02 set chain $J02_orig_chain
$J02 set x $J02_orig_x
$J02 set y $J02_orig_y
$J02 set z $J02_orig_z
$J02 set resid $J02_orig_resid
$J02 set resname $J02_orig_resname
$J02 set name $J02_orig_name
$J04 set chain $J04_orig_chain
$J04 set x $J04_orig_x
$J04 set y $J04_orig_y
$J04 set z $J04_orig_z
$J04 set resid $J04_orig_resid
$J04 set resname $J04_orig_resname
$J04 set name $J04_orig_name
################################ END SEGMENT J #################################
# Done with seg G
# Writing seg B
############################### BEGIN SEGMENT F ################################
set F00 [atomselect $m1 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
set F00_orig_chain [$F00 get chain]
set F00_orig_x [$F00 get x]
set F00_orig_y [$F00 get y]
set F00_orig_z [$F00 get z]
set F00_orig_resid [$F00 get resid]
set F00_orig_resname [$F00 get resname]
set F00_orig_name [$F00 get name]
$F00 set chain F
$F00 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$F00 writepdb PROTEIN_F_521_to_547.pdb
set F02 [atomselect $m1 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
set F02_orig_chain [$F02 get chain]
set F02_orig_x [$F02 get x]
set F02_orig_y [$F02 get y]
set F02_orig_z [$F02 get z]
set F02_orig_resid [$F02 get resid]
set F02_orig_resname [$F02 get resname]
set F02_orig_name [$F02 get name]
$F02 set chain F
$F02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$F02 writepdb PROTEIN_F_569_to_664.pdb
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
# Below reverts an engineered mutation:
    mutate 559 ILE
# Below reverts an engineered mutation:
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
$F00 set chain $F00_orig_chain
$F00 set x $F00_orig_x
$F00 set y $F00_orig_y
$F00 set z $F00_orig_z
$F00 set resid $F00_orig_resid
$F00 set resname $F00_orig_resname
$F00 set name $F00_orig_name
$F02 set chain $F02_orig_chain
$F02 set x $F02_orig_x
$F02 set y $F02_orig_y
$F02 set z $F02_orig_z
$F02 set resid $F02_orig_resid
$F02 set resname $F02_orig_resname
$F02 set name $F02_orig_name
################################ END SEGMENT F #################################
# Done with seg B
# Writing seg AG
###################### Not generating segment AG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg AG
# Writing seg CG
###################### Not generating segment CG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg CG
# Writing seg DG
###################### Not generating segment DG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg DG
# Writing seg GG
###################### Not generating segment GG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg GG
# Writing seg BG
###################### Not generating segment BG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg BG
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
# Writing seg G
############################### BEGIN SEGMENT O ################################
set O00 [atomselect $m1 "serial 1 to 1174"]
set O00_orig_chain [$O00 get chain]
set O00_orig_x [$O00 get x]
set O00_orig_y [$O00 get y]
set O00_orig_z [$O00 get z]
set O00_orig_resid [$O00 get resid]
set O00_orig_resname [$O00 get resname]
set O00_orig_name [$O00 get name]
$O00 set chain O
$O00 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$O00 writepdb PROTEIN_O_34_to_185.pdb
set O02 [atomselect $m1 "serial 1175 to 2795"]
set O02_orig_chain [$O02 get chain]
set O02_orig_x [$O02 get x]
set O02_orig_y [$O02 get y]
set O02_orig_z [$O02 get z]
set O02_orig_resid [$O02 get resid]
set O02_orig_resname [$O02 get resname]
set O02_orig_name [$O02 get name]
$O02 set chain O
$O02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$O02 writepdb PROTEIN_O_187_to_398.pdb
set O04 [atomselect $m1 "serial 2796 to 3543"]
set O04_orig_chain [$O04 get chain]
set O04_orig_x [$O04 get x]
set O04_orig_y [$O04 get y]
set O04_orig_z [$O04 get z]
set O04_orig_resid [$O04 get resid]
set O04_orig_resname [$O04 get resname]
set O04_orig_name [$O04 get name]
$O04 set chain O
$O04 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$O04 writepdb PROTEIN_O_411_to_505.pdb
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
# Below reverts an engineered mutation:
    mutate 332 THR
# Below reverts an engineered mutation:
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
$O00 set chain $O00_orig_chain
$O00 set x $O00_orig_x
$O00 set y $O00_orig_y
$O00 set z $O00_orig_z
$O00 set resid $O00_orig_resid
$O00 set resname $O00_orig_resname
$O00 set name $O00_orig_name
$O02 set chain $O02_orig_chain
$O02 set x $O02_orig_x
$O02 set y $O02_orig_y
$O02 set z $O02_orig_z
$O02 set resid $O02_orig_resid
$O02 set resname $O02_orig_resname
$O02 set name $O02_orig_name
$O04 set chain $O04_orig_chain
$O04 set x $O04_orig_x
$O04 set y $O04_orig_y
$O04 set z $O04_orig_z
$O04 set resid $O04_orig_resid
$O04 set resname $O04_orig_resname
$O04 set name $O04_orig_name
################################ END SEGMENT O #################################
# Done with seg G
# Writing seg B
############################### BEGIN SEGMENT L ################################
set L00 [atomselect $m1 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
set L00_orig_chain [$L00 get chain]
set L00_orig_x [$L00 get x]
set L00_orig_y [$L00 get y]
set L00_orig_z [$L00 get z]
set L00_orig_resid [$L00 get resid]
set L00_orig_resname [$L00 get resname]
set L00_orig_name [$L00 get name]
$L00 set chain L
$L00 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$L00 writepdb PROTEIN_L_521_to_547.pdb
set L02 [atomselect $m1 "serial 3722 to 4518"]
############ Atom with serial 4519 in PDB needs serial 4518 for VMD ############
set L02_orig_chain [$L02 get chain]
set L02_orig_x [$L02 get x]
set L02_orig_y [$L02 get y]
set L02_orig_z [$L02 get z]
set L02_orig_resid [$L02 get resid]
set L02_orig_resname [$L02 get resname]
set L02_orig_name [$L02 get name]
$L02 set chain L
$L02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$L02 writepdb PROTEIN_L_569_to_664.pdb
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
# Below reverts an engineered mutation:
    mutate 559 ILE
# Below reverts an engineered mutation:
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
$L00 set chain $L00_orig_chain
$L00 set x $L00_orig_x
$L00 set y $L00_orig_y
$L00 set z $L00_orig_z
$L00 set resid $L00_orig_resid
$L00 set resname $L00_orig_resname
$L00 set name $L00_orig_name
$L02 set chain $L02_orig_chain
$L02 set x $L02_orig_x
$L02 set y $L02_orig_y
$L02 set z $L02_orig_z
$L02 set resid $L02_orig_resid
$L02 set resname $L02_orig_resname
$L02 set name $L02_orig_name
################################ END SEGMENT L #################################
# Done with seg B
# Writing seg AG
###################### Not generating segment AG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg AG
# Writing seg CG
###################### Not generating segment CG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg CG
# Writing seg DG
###################### Not generating segment DG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg DG
# Writing seg GG
###################### Not generating segment GG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg GG
# Writing seg BG
###################### Not generating segment BG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################
# Done with seg BG
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
