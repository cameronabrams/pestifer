############################### BEGIN SEGMENT G ################################
set G01 [atomselect 1 "serial 1 to 1174"]
$G01 writepdb PROTEIN_G_34_to_185.pdb
set G03 [atomselect 1 "serial 1175 to 2795"]
$G03 writepdb PROTEIN_G_187_to_398.pdb
set G05 [atomselect 1 "serial 2796 to 3543"]
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
set B01 [atomselect 1 "serial 3544 to 3721"]
############ Atom with serial 3722 in PDB needs serial 3721 for VMD ############
$B01 writepdb PROTEIN_B_521_to_547.pdb
set B03 [atomselect 1 "serial 3722 to 4518"]
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
