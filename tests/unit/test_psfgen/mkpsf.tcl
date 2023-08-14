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
mol top $m1
############################## TRANSFORM 0 BEGINS ##############################
################ The following mappings of A.U. chains is used: ################
# A: A
# B: B
# C: C
# D: D
# E: E
# G: G
# H: H
# L: L
############################### SEGMENTS FOLLOW ################################
############################### BEGIN SEGMENT G ################################
set G01 [atomselect 1 "serial 1 to 1648"]
$G01 writepdb PROTEIN_G_90_to_396.pdb
set G03 [atomselect 1 "serial 1649 to 2300"]
$G03 writepdb PROTEIN_G_410_to_492.pdb
segment G {
    pdb PROTEIN_G_90_to_396.pdb
    residue 397 ASN G
    residue 398 SER G
    residue 399 THR G
    residue 400 TRP G
    residue 401 SER G
    residue 402 THR G
    residue 403 LYS G
    residue 404 GLY G
    residue 405 SER G
    residue 406 ASN G
    residue 407 ASN G
    residue 408 THR G
    residue 409 GLU G
    residue 409A GLY G
    pdb PROTEIN_G_410_to_492.pdb
# Below reverts a database conflict:
    mutate 403 GLU
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_G_90_to_396.pdb G
coordpdb PROTEIN_G_410_to_492.pdb G
######### Seeding orientation of model-built loop starting at G-ASN397 #########
coord G 397 N [cacoIn_nOut 396 G 0]
####################### Intra-segmental terminal patches #######################
patch CTER G:409
patch GLYP G:410
delatom G 409A
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT G #################################

############################### BEGIN SEGMENT C ################################
set C00 [atomselect 1 "serial 2301 to 3712"]
############ Atom with serial 3713 in PDB needs serial 3712 for VMD ############
$C00 writepdb PROTEIN_C_1_to_181.pdb
segment C {
    pdb PROTEIN_C_1_to_181.pdb
# Below reverts an engineered mutation:
    mutate 184 SER
# Below reverts an engineered mutation:
    mutate 185 ILE
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_C_1_to_181.pdb C
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT C #################################

############################### BEGIN SEGMENT L ################################
set L00 [atomselect 1 "serial 3713 to 5357"]
############ Atom with serial 5359 in PDB needs serial 5357 for VMD ############
$L00 writepdb PROTEIN_L_1_to_213.pdb
segment L {
    pdb PROTEIN_L_1_to_213.pdb
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_L_1_to_213.pdb L
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT L #################################

############################### BEGIN SEGMENT H ################################
set H00 [atomselect 1 "serial 5358 to 7080"]
############ Atom with serial 7083 in PDB needs serial 7080 for VMD ############
$H00 writepdb PROTEIN_H_1_to_229.pdb
segment H {
    pdb PROTEIN_H_1_to_229.pdb
}
###################### Coordinate-specification commands #######################
coordpdb PROTEIN_H_1_to_229.pdb H
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ END SEGMENT H #################################

################## Not generating generic segment GW (WATER) ###################

################## Not generating generic segment CW (WATER) ###################

################## Not generating generic segment LW (WATER) ###################

################## Not generating generic segment HW (WATER) ###################

###################### Not generating segment AG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################

###################### Not generating segment BG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################

###################### Not generating segment DG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################

###################### Not generating segment EG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################

###################### Not generating segment GG (GLYCAN) ######################
####################### GLYCANS ARE NOT YET IMPLEMENTED ########################

############################# DISU PATCHES FOLLOW ##############################
patch DISU G:119 G:205
patch DISU G:126 G:196
patch DISU G:218 G:247
patch DISU G:228 G:239
patch DISU G:296 G:331
patch DISU G:331 G:385
patch DISU G:378 G:445
patch DISU G:385 G:418
patch DISU C:16 C:84
patch DISU C:130 C:159
patch DISU L:23 L:88
patch DISU L:136 L:196
patch DISU H:22 H:96
patch DISU H:155 H:211
############################# LINK PATCHES FOLLOW ##############################
############################### TRANSFORM 0 ENDS ###############################

exit
############################# END PESTIFER PSFGEN ##############################
######################## Thank you for using pestifer! #########################
