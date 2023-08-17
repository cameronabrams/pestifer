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
pdbalias residue MAN AMAN
pdbalias residue BMA BMAN
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
# Subsegment [0] is a resolved run
coordpdb PROTEIN_G_34_to_185.pdb G
# Subsegment [2] is a resolved run
coordpdb PROTEIN_G_187_to_398.pdb G
# Subsegment [4] is a resolved run
coordpdb PROTEIN_G_411_to_505.pdb G
# Subsegment [1]/5 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at G-GLU185A from G-ASN185
coord G 185A N [cacoIn_nOut 185 G $m1])
# Subsegment [3]/5 is a missing loop
# ...attached to subsegment 2
# Seeding orientation of model-built loop starting at G-THR400 from G-ASN398
coord G 400 N [cacoIn_nOut 398 G $m1])
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
# Subsegment [0] is a resolved run
coordpdb PROTEIN_B_521_to_547.pdb B
# Subsegment [2] is a resolved run
coordpdb PROTEIN_B_569_to_664.pdb B
# Subsegment [1]/3 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B-ILE548 from B-GLY547
coord B 548 N [cacoIn_nOut 547 B $m1])
####################### Intra-segmental terminal patches #######################
patch CTER B:568
patch NTER B:569
delatom B 568A
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT B #################################
# Done with seg B
# Writing seg AG
############################## BEGIN SEGMENT AG02 ##############################
# This segment is referenced as chain A in the input structure
set AG02 [atomselect $m1 "serial 4519 to 4590"]
$AG02 writepdb GENERIC_AG02.pdb
segment AG02 {
    pdb GENERIC_AG02.pdb
}
coordpdb GENERIC_AG02.pdb AG02
############################### END SEGMENT AG02 ###############################
# Done with seg AG
# Writing seg CG
############################## BEGIN SEGMENT CG02 ##############################
# This segment is referenced as chain C in the input structure
set CG02 [atomselect $m1 "serial 4591 to 4618"]
$CG02 writepdb GENERIC_CG02.pdb
segment CG02 {
    pdb GENERIC_CG02.pdb
}
coordpdb GENERIC_CG02.pdb CG02
############################### END SEGMENT CG02 ###############################
# Done with seg CG
# Writing seg GG
############################## BEGIN SEGMENT GG02 ##############################
# This segment is referenced as chain G in the input structure
set GG02 [atomselect $m1 "serial 4647 to 4814"]
$GG02 writepdb GENERIC_GG02.pdb
segment GG02 {
    pdb GENERIC_GG02.pdb
}
coordpdb GENERIC_GG02.pdb GG02
############################### END SEGMENT GG02 ###############################
# Done with seg GG
# Writing seg BG
############################## BEGIN SEGMENT BG02 ##############################
# This segment is referenced as chain B in the input structure
set BG02 [atomselect $m1 "serial 4815 to 4856"]
$BG02 writepdb GENERIC_BG02.pdb
segment BG02 {
    pdb GENERIC_BG02.pdb
}
coordpdb GENERIC_BG02.pdb BG02
############################### END SEGMENT BG02 ###############################
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
patch NGLB G:156 GG02:615
patch NGLB G:160 GG02:616
patch NGLB G:197 CG02:1
patch NGLB G:234 GG02:601
patch NGLB G:262 AG02:1
patch NGLB G:276 GG02:608
patch NGLB G:295 GG02:619
patch NGLB G:301 GG02:620
patch NGLB G:339 GG02:609
patch NGLB G:355 GG02:610
patch NGLB G:363 GG02:611
patch NGLB G:386 GG02:612
patch NGLB G:392 GG02:613
patch NGLB G:448 GG02:614
patch NGLB B:611 BG02:701
patch NGLB B:618 BG02:702
patch NGLB B:637 BG02:703
set cn 4
set abi [axeq 2 0 A C1 1]
set abj [axeq 1 0 A O4 -1]
set pres "1$cn$abi$abj"
patch $pres AG02:1 AG02:1
set cn 4
set abi [axeq 3 0 A C1 2]
set abj [axeq 2 0 A O4 -1]
set pres "1$cn$abi$abj"
patch $pres AG02:2 AG02:2
set cn 3
set abi [axeq 4 0 A C1 3]
set abj [axeq 3 0 A O3 -1]
set pres "1$cn$abi$abj"
patch $pres AG02:3 AG02:3
set cn 6
set abi [axeq 6 0 A C1 3]
set abj [axeq 3 0 A O6 -1]
if { $abi == "a" } { set abi A }
if { $abi == "b" } { set abi B }
if { $abj == "b" } { set abj T }
set pres "1$cn$abi$abj"
patch $pres AG02:3 AG02:3
set cn 2
set abi [axeq 5 0 A C1 4]
set abj [axeq 4 0 A O2 -1]
set pres "1$cn$abi$abj"
patch $pres AG02:4 AG02:4
set cn 4
set abi [axeq 2 0 C C1 1]
set abj [axeq 1 0 C O4 -1]
set pres "1$cn$abi$abj"
patch $pres CG02:1 CG02:1
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
# Subsegment [0] is a resolved run
coordpdb PROTEIN_J_34_to_185.pdb J
# Subsegment [2] is a resolved run
coordpdb PROTEIN_J_187_to_398.pdb J
# Subsegment [4] is a resolved run
coordpdb PROTEIN_J_411_to_505.pdb J
# Subsegment [1]/5 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at G-GLU185A from G-ASN185
coord J 185A N [cacoIn_nOut 185 J $m1])
# Subsegment [3]/5 is a missing loop
# ...attached to subsegment 2
# Seeding orientation of model-built loop starting at G-THR400 from G-ASN398
coord J 400 N [cacoIn_nOut 398 J $m1])
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
# Subsegment [0] is a resolved run
coordpdb PROTEIN_F_521_to_547.pdb F
# Subsegment [2] is a resolved run
coordpdb PROTEIN_F_569_to_664.pdb F
# Subsegment [1]/3 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B-ILE548 from B-GLY547
coord F 548 N [cacoIn_nOut 547 F $m1])
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
############################## BEGIN SEGMENT EG02 ##############################
# This segment is referenced as chain A in the input structure
set EG02 [atomselect $m1 "serial 4519 to 4590"]
set EG02_orig_chain [$EG02 get chain]
set EG02_orig_x [$EG02 get x]
set EG02_orig_y [$EG02 get y]
set EG02_orig_z [$EG02 get z]
set EG02_orig_resid [$EG02 get resid]
set EG02_orig_resname [$EG02 get resname]
set EG02_orig_name [$EG02 get name]
$EG02 set chain E
$EG02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$EG02 writepdb GENERIC_EG02.pdb
segment EG02 {
    pdb GENERIC_EG02.pdb
}
coordpdb GENERIC_EG02.pdb EG02
######################## Restoring A.U. state for EG02 #########################
$EG02 set chain $EG02_orig_chain
$EG02 set x $EG02_orig_x
$EG02 set y $EG02_orig_y
$EG02 set z $EG02_orig_z
$EG02 set resid $EG02_orig_resid
$EG02 set resname $EG02_orig_resname
$EG02 set name $EG02_orig_name
############################### END SEGMENT EG02 ###############################
# Done with seg AG
# Writing seg CG
############################## BEGIN SEGMENT HG02 ##############################
# This segment is referenced as chain C in the input structure
set HG02 [atomselect $m1 "serial 4591 to 4618"]
set HG02_orig_chain [$HG02 get chain]
set HG02_orig_x [$HG02 get x]
set HG02_orig_y [$HG02 get y]
set HG02_orig_z [$HG02 get z]
set HG02_orig_resid [$HG02 get resid]
set HG02_orig_resname [$HG02 get resname]
set HG02_orig_name [$HG02 get name]
$HG02 set chain H
$HG02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$HG02 writepdb GENERIC_HG02.pdb
segment HG02 {
    pdb GENERIC_HG02.pdb
}
coordpdb GENERIC_HG02.pdb HG02
######################## Restoring A.U. state for HG02 #########################
$HG02 set chain $HG02_orig_chain
$HG02 set x $HG02_orig_x
$HG02 set y $HG02_orig_y
$HG02 set z $HG02_orig_z
$HG02 set resid $HG02_orig_resid
$HG02 set resname $HG02_orig_resname
$HG02 set name $HG02_orig_name
############################### END SEGMENT HG02 ###############################
# Done with seg CG
# Writing seg GG
############################## BEGIN SEGMENT JG02 ##############################
# This segment is referenced as chain G in the input structure
set JG02 [atomselect $m1 "serial 4647 to 4814"]
set JG02_orig_chain [$JG02 get chain]
set JG02_orig_x [$JG02 get x]
set JG02_orig_y [$JG02 get y]
set JG02_orig_z [$JG02 get z]
set JG02_orig_resid [$JG02 get resid]
set JG02_orig_resname [$JG02 get resname]
set JG02_orig_name [$JG02 get name]
$JG02 set chain J
$JG02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$JG02 writepdb GENERIC_JG02.pdb
segment JG02 {
    pdb GENERIC_JG02.pdb
}
coordpdb GENERIC_JG02.pdb JG02
######################## Restoring A.U. state for JG02 #########################
$JG02 set chain $JG02_orig_chain
$JG02 set x $JG02_orig_x
$JG02 set y $JG02_orig_y
$JG02 set z $JG02_orig_z
$JG02 set resid $JG02_orig_resid
$JG02 set resname $JG02_orig_resname
$JG02 set name $JG02_orig_name
############################### END SEGMENT JG02 ###############################
# Done with seg GG
# Writing seg BG
############################## BEGIN SEGMENT FG02 ##############################
# This segment is referenced as chain B in the input structure
set FG02 [atomselect $m1 "serial 4815 to 4856"]
set FG02_orig_chain [$FG02 get chain]
set FG02_orig_x [$FG02 get x]
set FG02_orig_y [$FG02 get y]
set FG02_orig_z [$FG02 get z]
set FG02_orig_resid [$FG02 get resid]
set FG02_orig_resname [$FG02 get resname]
set FG02_orig_name [$FG02 get name]
$FG02 set chain F
$FG02 move { { -0.5 -0.866025 0.0 107.18  } { 0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$FG02 writepdb GENERIC_FG02.pdb
segment FG02 {
    pdb GENERIC_FG02.pdb
}
coordpdb GENERIC_FG02.pdb FG02
######################## Restoring A.U. state for FG02 #########################
$FG02 set chain $FG02_orig_chain
$FG02 set x $FG02_orig_x
$FG02 set y $FG02_orig_y
$FG02 set z $FG02_orig_z
$FG02 set resid $FG02_orig_resid
$FG02 set resname $FG02_orig_resname
$FG02 set name $FG02_orig_name
############################### END SEGMENT FG02 ###############################
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
patch NGLB J:156 JG02:615
patch NGLB J:160 JG02:616
patch NGLB J:197 HG02:1
patch NGLB J:234 JG02:601
patch NGLB J:262 EG02:1
patch NGLB J:276 JG02:608
patch NGLB J:295 JG02:619
patch NGLB J:301 JG02:620
patch NGLB J:339 JG02:609
patch NGLB J:355 JG02:610
patch NGLB J:363 JG02:611
patch NGLB J:386 JG02:612
patch NGLB J:392 JG02:613
patch NGLB J:448 JG02:614
patch NGLB F:611 FG02:701
patch NGLB F:618 FG02:702
patch NGLB F:637 FG02:703
set cn 4
set abi [axeq 2 0 A C1 1]
set abj [axeq 1 0 A O4 -1]
set pres "1$cn$abi$abj"
patch $pres EG02:1 EG02:1
set cn 4
set abi [axeq 3 0 A C1 2]
set abj [axeq 2 0 A O4 -1]
set pres "1$cn$abi$abj"
patch $pres EG02:2 EG02:2
set cn 3
set abi [axeq 4 0 A C1 3]
set abj [axeq 3 0 A O3 -1]
set pres "1$cn$abi$abj"
patch $pres EG02:3 EG02:3
set cn 6
set abi [axeq 6 0 A C1 3]
set abj [axeq 3 0 A O6 -1]
if { $abi == "a" } { set abi A }
if { $abi == "b" } { set abi B }
if { $abj == "b" } { set abj T }
set pres "1$cn$abi$abj"
patch $pres EG02:3 EG02:3
set cn 2
set abi [axeq 5 0 A C1 4]
set abj [axeq 4 0 A O2 -1]
set pres "1$cn$abi$abj"
patch $pres EG02:4 EG02:4
set cn 4
set abi [axeq 2 0 C C1 1]
set abj [axeq 1 0 C O4 -1]
set pres "1$cn$abi$abj"
patch $pres HG02:1 HG02:1
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
# Subsegment [0] is a resolved run
coordpdb PROTEIN_O_34_to_185.pdb O
# Subsegment [2] is a resolved run
coordpdb PROTEIN_O_187_to_398.pdb O
# Subsegment [4] is a resolved run
coordpdb PROTEIN_O_411_to_505.pdb O
# Subsegment [1]/5 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at G-GLU185A from G-ASN185
coord O 185A N [cacoIn_nOut 185 O $m1])
# Subsegment [3]/5 is a missing loop
# ...attached to subsegment 2
# Seeding orientation of model-built loop starting at G-THR400 from G-ASN398
coord O 400 N [cacoIn_nOut 398 O $m1])
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
# Subsegment [0] is a resolved run
coordpdb PROTEIN_L_521_to_547.pdb L
# Subsegment [2] is a resolved run
coordpdb PROTEIN_L_569_to_664.pdb L
# Subsegment [1]/3 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B-ILE548 from B-GLY547
coord L 548 N [cacoIn_nOut 547 L $m1])
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
############################## BEGIN SEGMENT KG02 ##############################
# This segment is referenced as chain A in the input structure
set KG02 [atomselect $m1 "serial 4519 to 4590"]
set KG02_orig_chain [$KG02 get chain]
set KG02_orig_x [$KG02 get x]
set KG02_orig_y [$KG02 get y]
set KG02_orig_z [$KG02 get z]
set KG02_orig_resid [$KG02 get resid]
set KG02_orig_resname [$KG02 get resname]
set KG02_orig_name [$KG02 get name]
$KG02 set chain K
$KG02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$KG02 writepdb GENERIC_KG02.pdb
segment KG02 {
    pdb GENERIC_KG02.pdb
}
coordpdb GENERIC_KG02.pdb KG02
######################## Restoring A.U. state for KG02 #########################
$KG02 set chain $KG02_orig_chain
$KG02 set x $KG02_orig_x
$KG02 set y $KG02_orig_y
$KG02 set z $KG02_orig_z
$KG02 set resid $KG02_orig_resid
$KG02 set resname $KG02_orig_resname
$KG02 set name $KG02_orig_name
############################### END SEGMENT KG02 ###############################
# Done with seg AG
# Writing seg CG
############################## BEGIN SEGMENT MG02 ##############################
# This segment is referenced as chain C in the input structure
set MG02 [atomselect $m1 "serial 4591 to 4618"]
set MG02_orig_chain [$MG02 get chain]
set MG02_orig_x [$MG02 get x]
set MG02_orig_y [$MG02 get y]
set MG02_orig_z [$MG02 get z]
set MG02_orig_resid [$MG02 get resid]
set MG02_orig_resname [$MG02 get resname]
set MG02_orig_name [$MG02 get name]
$MG02 set chain M
$MG02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$MG02 writepdb GENERIC_MG02.pdb
segment MG02 {
    pdb GENERIC_MG02.pdb
}
coordpdb GENERIC_MG02.pdb MG02
######################## Restoring A.U. state for MG02 #########################
$MG02 set chain $MG02_orig_chain
$MG02 set x $MG02_orig_x
$MG02 set y $MG02_orig_y
$MG02 set z $MG02_orig_z
$MG02 set resid $MG02_orig_resid
$MG02 set resname $MG02_orig_resname
$MG02 set name $MG02_orig_name
############################### END SEGMENT MG02 ###############################
# Done with seg CG
# Writing seg GG
############################## BEGIN SEGMENT OG02 ##############################
# This segment is referenced as chain G in the input structure
set OG02 [atomselect $m1 "serial 4647 to 4814"]
set OG02_orig_chain [$OG02 get chain]
set OG02_orig_x [$OG02 get x]
set OG02_orig_y [$OG02 get y]
set OG02_orig_z [$OG02 get z]
set OG02_orig_resid [$OG02 get resid]
set OG02_orig_resname [$OG02 get resname]
set OG02_orig_name [$OG02 get name]
$OG02 set chain O
$OG02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$OG02 writepdb GENERIC_OG02.pdb
segment OG02 {
    pdb GENERIC_OG02.pdb
}
coordpdb GENERIC_OG02.pdb OG02
######################## Restoring A.U. state for OG02 #########################
$OG02 set chain $OG02_orig_chain
$OG02 set x $OG02_orig_x
$OG02 set y $OG02_orig_y
$OG02 set z $OG02_orig_z
$OG02 set resid $OG02_orig_resid
$OG02 set resname $OG02_orig_resname
$OG02 set name $OG02_orig_name
############################### END SEGMENT OG02 ###############################
# Done with seg GG
# Writing seg BG
############################## BEGIN SEGMENT LG02 ##############################
# This segment is referenced as chain B in the input structure
set LG02 [atomselect $m1 "serial 4815 to 4856"]
set LG02_orig_chain [$LG02 get chain]
set LG02_orig_x [$LG02 get x]
set LG02_orig_y [$LG02 get y]
set LG02_orig_z [$LG02 get z]
set LG02_orig_resid [$LG02 get resid]
set LG02_orig_resname [$LG02 get resname]
set LG02_orig_name [$LG02 get name]
$LG02 set chain L
$LG02 move { { -0.5 0.866025 0.0 -107.18  } { -0.866025 -0.5 0.0 185.64121  } { 0.0 0.0 1.0 0.0  } { 0 0 0 1 } }
$LG02 writepdb GENERIC_LG02.pdb
segment LG02 {
    pdb GENERIC_LG02.pdb
}
coordpdb GENERIC_LG02.pdb LG02
######################## Restoring A.U. state for LG02 #########################
$LG02 set chain $LG02_orig_chain
$LG02 set x $LG02_orig_x
$LG02 set y $LG02_orig_y
$LG02 set z $LG02_orig_z
$LG02 set resid $LG02_orig_resid
$LG02 set resname $LG02_orig_resname
$LG02 set name $LG02_orig_name
############################### END SEGMENT LG02 ###############################
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
patch NGLB O:156 OG02:615
patch NGLB O:160 OG02:616
patch NGLB O:197 MG02:1
patch NGLB O:234 OG02:601
patch NGLB O:262 KG02:1
patch NGLB O:276 OG02:608
patch NGLB O:295 OG02:619
patch NGLB O:301 OG02:620
patch NGLB O:339 OG02:609
patch NGLB O:355 OG02:610
patch NGLB O:363 OG02:611
patch NGLB O:386 OG02:612
patch NGLB O:392 OG02:613
patch NGLB O:448 OG02:614
patch NGLB L:611 LG02:701
patch NGLB L:618 LG02:702
patch NGLB L:637 LG02:703
set cn 4
set abi [axeq 2 0 A C1 1]
set abj [axeq 1 0 A O4 -1]
set pres "1$cn$abi$abj"
patch $pres KG02:1 KG02:1
set cn 4
set abi [axeq 3 0 A C1 2]
set abj [axeq 2 0 A O4 -1]
set pres "1$cn$abi$abj"
patch $pres KG02:2 KG02:2
set cn 3
set abi [axeq 4 0 A C1 3]
set abj [axeq 3 0 A O3 -1]
set pres "1$cn$abi$abj"
patch $pres KG02:3 KG02:3
set cn 6
set abi [axeq 6 0 A C1 3]
set abj [axeq 3 0 A O6 -1]
if { $abi == "a" } { set abi A }
if { $abi == "b" } { set abi B }
if { $abj == "b" } { set abj T }
set pres "1$cn$abi$abj"
patch $pres KG02:3 KG02:3
set cn 2
set abi [axeq 5 0 A C1 4]
set abj [axeq 4 0 A O2 -1]
set pres "1$cn$abi$abj"
patch $pres KG02:4 KG02:4
set cn 4
set abi [axeq 2 0 C C1 1]
set abj [axeq 1 0 C O4 -1]
set pres "1$cn$abi$abj"
patch $pres MG02:1 MG02:1
############################### TRANSFORM 2 ENDS ###############################
guesscoord
regenerate angles dihedrals
writepsf cmap step1-build.psf
writepdb step1-build.pdb
exit
############################# END PESTIFER PSFGEN ##############################
######################## Thank you for using pestifer! #########################
