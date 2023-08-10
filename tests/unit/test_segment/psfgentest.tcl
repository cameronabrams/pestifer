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
