# pestifer.scripters: 00-01-00_psfgen-build.tcl
####################### Created Thu Apr 16 16:05:21 2026 #######################
package require PestiferCRot
namespace import PestiferCRot::*
package require psfgen
psfcontext mixedcase
topology top_all36_prot.rtf
topology top_all36_cgenff.rtf
topology top_all36_lipid.rtf
topology top_all36_carb.rtf
topology top_all36_na.rtf
topology toppar_water_ions.str
topology toppar_all36_carb_glycopeptide.str
topology toppar_all36_prot_modify_res.str
topology toppar_all36_moreions.str
topology top_all35_ethers.rtf
pdbalias atom ILE CD1 CD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5 C10 C
pdbalias atom ANE5 C11 CT
pdbalias atom ANE5 N5 N
pdbalias atom ANE5 O1A O11
pdbalias atom ANE5 O1B O12
pdbalias atom ANE5 O10 O
pdbalias atom VCG C01 C1
pdbalias atom VCG C02 C2
pdbalias atom VCG C03 C3
pdbalias atom VCG C04 C4
pdbalias atom VCG C05 C5
pdbalias atom VCG C06 C6
pdbalias atom VCG C07 C7
pdbalias atom VCG C08 C8
pdbalias atom VCG C09 C9
pdbalias atom TIP3 O OH2
pdbalias atom ILE CD1 CD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5 C10 C
pdbalias atom ANE5 C11 CT
pdbalias atom ANE5 N5 N
pdbalias atom ANE5 O1A O11
pdbalias atom ANE5 O1B O12
pdbalias atom ANE5 O10 O
pdbalias atom VCG C01 C1
pdbalias atom VCG C02 C2
pdbalias atom VCG C03 C3
pdbalias atom VCG C04 C4
pdbalias atom VCG C05 C5
pdbalias atom VCG C06 C6
pdbalias atom VCG C07 C7
pdbalias atom VCG C08 C8
pdbalias atom VCG C09 C9
pdbalias atom TIP3 O OH2
pdbalias residue HIS HSD
pdbalias residue PO4 H2PO4
pdbalias residue H2PO H2PO4
pdbalias residue MAN AMAN
pdbalias residue BMA BMAN
pdbalias residue BGLC BGLCNA
pdbalias residue NAG BGLCNA
pdbalias residue NDG AGLCNA
pdbalias residue FUC AFUC
pdbalias residue FUL BFUC
pdbalias residue GAL BGAL
pdbalias residue SIA ANE5AC
pdbalias residue ANE5 ANE5AC
pdbalias residue EIC LIN
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias residue CL CLA
pdbalias residue C6DH C6DHPC
pdbalias residue C7DH C7DHPC
pdbalias residue DT THY
pdbalias residue DA ADE
pdbalias residue DC CYT
pdbalias residue DG GUA
pdbalias residue DU URA
pdbalias residue HEM HEME
pdbalias residue TOCL TOCL1
pdbalias residue HIS HSD
pdbalias residue PO4 H2PO4
pdbalias residue H2PO H2PO4
pdbalias residue MAN AMAN
pdbalias residue BMA BMAN
pdbalias residue BGLC BGLCNA
pdbalias residue NAG BGLCNA
pdbalias residue NDG AGLCNA
pdbalias residue FUC AFUC
pdbalias residue FUL BFUC
pdbalias residue GAL BGAL
pdbalias residue SIA ANE5AC
pdbalias residue ANE5 ANE5AC
pdbalias residue EIC LIN
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias residue CL CLA
pdbalias residue C6DH C6DHPC
pdbalias residue C7DH C7DHPC
pdbalias residue DT THY
pdbalias residue DA ADE
pdbalias residue DC CYT
pdbalias residue DG GUA
pdbalias residue DU URA
pdbalias residue HEM HEME
pdbalias residue TOCL TOCL1
mol new TMD-N2-open-sym-last-frame-dry.pdb waitfor all autobonds off
set m2 [molinfo top get id]
set nf [molinfo $m2 get numframes]
if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $m2 }
mol top $m2
############################## Transform 1 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
############################### Segments follow ################################
# Segment A begins
############################### Segment A begins ###############################
set A00 [atomselect $m2 "serial 1 to 2795"]
$A00 set segname A
$A00 writepdb segtype_polymer_A_515_to_685.pdb
segment A {
    pdb segtype_polymer_A_515_to_685.pdb
    residue 686 ILE A
    residue 687 ILE A
    residue 688 ILE A
    residue 689 VAL A
    residue 690 GLY A
    residue 691 SER A
    residue 692 LEU A
    residue 693 ILE A
    residue 694 GLY A
    residue 695 LEU A
    residue 696 ARG A
    residue 697 ILE A
    residue 698 VAL A
    residue 699 PHE A
    residue 700 ALA A
    residue 701 VAL A
    residue 702 LEU A
    residue 703 SER A
    residue 704 LEU A
    residue 705 VAL A
    residue 706 ASN A
    residue 707 ARG A
    residue 708 VAL A
    residue 709 ARG A
    residue 710 GLN A
}
################################ End segment A #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_A_515_to_685.pdb A
# Subsegment [1]/2 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at A_ILE686 from A_PHE685
coord A 686 N {140.08033 135.32376 55.34191}
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment A ends ################################
# Segment A ends
# Segment B begins
############################### Segment B begins ###############################
set B00 [atomselect $m2 "serial 11992 to 14786"]
$B00 set segname B
$B00 writepdb segtype_polymer_B_515_to_685.pdb
segment B {
    pdb segtype_polymer_B_515_to_685.pdb
    residue 686 ILE B
    residue 687 ILE B
    residue 688 ILE B
    residue 689 VAL B
    residue 690 GLY B
    residue 691 SER B
    residue 692 LEU B
    residue 693 ILE B
    residue 694 GLY B
    residue 695 LEU B
    residue 696 ARG B
    residue 697 ILE B
    residue 698 VAL B
    residue 699 PHE B
    residue 700 ALA B
    residue 701 VAL B
    residue 702 LEU B
    residue 703 SER B
    residue 704 LEU B
    residue 705 VAL B
    residue 706 ASN B
    residue 707 ARG B
    residue 708 VAL B
    residue 709 ARG B
    residue 710 GLN B
}
################################ End segment B #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_B_515_to_685.pdb B
# Subsegment [1]/2 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at B_ILE686 from B_PHE685
coord B 686 N {130.99358 136.17478 52.53613}
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment B ends ################################
# Segment B ends
# Segment C begins
############################### Segment C begins ###############################
set C00 [atomselect $m2 "serial 9223 to 11991"]
$C00 set segname C
$C00 writepdb segtype_polymer_C_1_to_176.pdb
segment C {
    pdb segtype_polymer_C_1_to_176.pdb
}
################################ End segment C #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_C_1_to_176.pdb C
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment C ends ################################
# Segment C ends
# Segment D begins
############################### Segment D begins ###############################
set D00 [atomselect $m2 "serial 23983 to 26777"]
$D00 set segname D
$D00 writepdb segtype_polymer_D_515_to_685.pdb
segment D {
    pdb segtype_polymer_D_515_to_685.pdb
    residue 686 ILE D
    residue 687 ILE D
    residue 688 ILE D
    residue 689 VAL D
    residue 690 GLY D
    residue 691 SER D
    residue 692 LEU D
    residue 693 ILE D
    residue 694 GLY D
    residue 695 LEU D
    residue 696 ARG D
    residue 697 ILE D
    residue 698 VAL D
    residue 699 PHE D
    residue 700 ALA D
    residue 701 VAL D
    residue 702 LEU D
    residue 703 SER D
    residue 704 LEU D
    residue 705 VAL D
    residue 706 ASN D
    residue 707 ARG D
    residue 708 VAL D
    residue 709 ARG D
    residue 710 GLN D
}
################################ End segment D #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_D_515_to_685.pdb D
# Subsegment [1]/2 is a missing loop
# ...attached to subsegment 0
# Seeding orientation of model-built loop starting at D_ILE686 from D_PHE685
coord D 686 N {127.57966 128.61762 54.12308}
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment D ends ################################
# Segment D ends
# Segment E begins
############################### Segment E begins ###############################
set E00 [atomselect $m2 "serial 21214 to 23982"]
$E00 set segname E
$E00 writepdb segtype_polymer_E_1_to_176.pdb
segment E {
    pdb segtype_polymer_E_1_to_176.pdb
}
################################ End segment E #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_E_1_to_176.pdb E
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment E ends ################################
# Segment E ends
# Segment F begins
############################### Segment F begins ###############################
set F00 [atomselect $m2 "serial 33205 to 35973"]
$F00 set segname F
$F00 writepdb segtype_polymer_F_1_to_176.pdb
segment F {
    pdb segtype_polymer_F_1_to_176.pdb
}
################################ End segment F #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_F_1_to_176.pdb F
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment F ends ################################
# Segment F ends
# Segment G begins
############################### Segment G begins ###############################
set G00 [atomselect $m2 "serial 2796 to 9222"]
$G00 set segname G
$G00 writepdb segtype_polymer_G_32_to_511.pdb
segment G {
    pdb segtype_polymer_G_32_to_511.pdb
}
################################ End segment G #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_G_32_to_511.pdb G
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment G ends ################################
# Segment G ends
# Segment H begins
############################### Segment H begins ###############################
set H00 [atomselect $m2 "serial 38720 to 42118"]
$H00 set segname H
$H00 writepdb segtype_polymer_H_1_to_213.pdb
segment H {
    pdb segtype_polymer_H_1_to_213.pdb
}
################################ End segment H #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_H_1_to_213.pdb H
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment H ends ################################
# Segment H ends
# Segment I begins
############################### Segment I begins ###############################
set I00 [atomselect $m2 "serial 14787 to 21213"]
$I00 set segname I
$I00 writepdb segtype_polymer_I_32_to_511.pdb
segment I {
    pdb segtype_polymer_I_32_to_511.pdb
}
################################ End segment I #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_I_32_to_511.pdb I
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment I ends ################################
# Segment I ends
# Segment J begins
############################### Segment J begins ###############################
set J00 [atomselect $m2 "serial 26778 to 33204"]
$J00 set segname J
$J00 writepdb segtype_polymer_J_32_to_511.pdb
segment J {
    pdb segtype_polymer_J_32_to_511.pdb
}
################################ End segment J #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_J_32_to_511.pdb J
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment J ends ################################
# Segment J ends
# Segment K begins
############################### Segment K begins ###############################
set K00 [atomselect $m2 "serial 42231 to 45629"]
$K00 set segname K
$K00 writepdb segtype_polymer_K_1_to_213.pdb
segment K {
    pdb segtype_polymer_K_1_to_213.pdb
}
################################ End segment K #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_K_1_to_213.pdb K
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment K ends ################################
# Segment K ends
# Segment L begins
############################### Segment L begins ###############################
set L00 [atomselect $m2 "serial 45742 to 48985"]
$L00 set segname L
$L00 writepdb segtype_polymer_L_1_to_212.pdb
segment L {
    pdb segtype_polymer_L_1_to_212.pdb
}
################################ End segment L #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_L_1_to_212.pdb L
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment L ends ################################
# Segment L ends
# Segment M begins
############################### Segment M begins ###############################
set M00 [atomselect $m2 "serial 49098 to 52496"]
$M00 set segname M
$M00 writepdb segtype_polymer_M_1_to_213.pdb
segment M {
    pdb segtype_polymer_M_1_to_213.pdb
}
################################ End segment M #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_M_1_to_213.pdb M
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment M ends ################################
# Segment M ends
# Segment N begins
############################### Segment N begins ###############################
set N00 [atomselect $m2 "serial 52609 to 55852"]
$N00 set segname N
$N00 writepdb segtype_polymer_N_1_to_212.pdb
segment N {
    pdb segtype_polymer_N_1_to_212.pdb
}
################################ End segment N #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_N_1_to_212.pdb N
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment N ends ################################
# Segment N ends
# Segment O begins
############################### Segment O begins ###############################
set O00 [atomselect $m2 "serial 55965 to 59208"]
$O00 set segname O
$O00 writepdb segtype_polymer_O_1_to_212.pdb
segment O {
    pdb segtype_polymer_O_1_to_212.pdb
}
################################ End segment O #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_O_1_to_212.pdb O
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment O ends ################################
# Segment O ends
# Segment XA1 begins
###################### Segment XA1 begins as image of XA1 ######################
set XA1 [atomselect $m2 "serial 38608 to 38635"]
$XA1 set segname XA1
$XA1 set resid [list 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701]
$XA1 writepdb segtype_generic_XA1.pdb
segment XA1 {
    first none
    last none
    pdb segtype_generic_XA1.pdb
}
coordpdb segtype_generic_XA1.pdb XA1
############################### Segment XA1 ends ###############################
# Segment XA1 ends
# Segment XA2 begins
###################### Segment XA2 begins as image of XA2 ######################
set XA2 [atomselect $m2 "serial 38636 to 38663"]
$XA2 set segname XA2
$XA2 set resid [list 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702]
$XA2 writepdb segtype_generic_XA2.pdb
segment XA2 {
    first none
    last none
    pdb segtype_generic_XA2.pdb
}
coordpdb segtype_generic_XA2.pdb XA2
############################### Segment XA2 ends ###############################
# Segment XA2 ends
# Segment XA3 begins
###################### Segment XA3 begins as image of XA3 ######################
set XA3 [atomselect $m2 "serial 38664 to 38691"]
$XA3 set segname XA3
$XA3 set resid [list 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703]
$XA3 writepdb segtype_generic_XA3.pdb
segment XA3 {
    first none
    last none
    pdb segtype_generic_XA3.pdb
}
coordpdb segtype_generic_XA3.pdb XA3
############################### Segment XA3 ends ###############################
# Segment XA3 ends
# Segment XA4 begins
###################### Segment XA4 begins as image of XA4 ######################
set XA4 [atomselect $m2 "serial 38692 to 38719"]
$XA4 set segname XA4
$XA4 set resid [list 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704]
$XA4 writepdb segtype_generic_XA4.pdb
segment XA4 {
    first none
    last none
    pdb segtype_generic_XA4.pdb
}
coordpdb segtype_generic_XA4.pdb XA4
############################### Segment XA4 ends ###############################
# Segment XA4 ends
# Segment XB1 begins
###################### Segment XB1 begins as image of XB1 ######################
set XB1 [atomselect $m2 "serial 45630 to 45657"]
$XB1 set segname XB1
$XB1 set resid [list 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701]
$XB1 writepdb segtype_generic_XB1.pdb
segment XB1 {
    first none
    last none
    pdb segtype_generic_XB1.pdb
}
coordpdb segtype_generic_XB1.pdb XB1
############################### Segment XB1 ends ###############################
# Segment XB1 ends
# Segment XB2 begins
###################### Segment XB2 begins as image of XB2 ######################
set XB2 [atomselect $m2 "serial 45658 to 45685"]
$XB2 set segname XB2
$XB2 set resid [list 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702]
$XB2 writepdb segtype_generic_XB2.pdb
segment XB2 {
    first none
    last none
    pdb segtype_generic_XB2.pdb
}
coordpdb segtype_generic_XB2.pdb XB2
############################### Segment XB2 ends ###############################
# Segment XB2 ends
# Segment XB3 begins
###################### Segment XB3 begins as image of XB3 ######################
set XB3 [atomselect $m2 "serial 45686 to 45713"]
$XB3 set segname XB3
$XB3 set resid [list 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703]
$XB3 writepdb segtype_generic_XB3.pdb
segment XB3 {
    first none
    last none
    pdb segtype_generic_XB3.pdb
}
coordpdb segtype_generic_XB3.pdb XB3
############################### Segment XB3 ends ###############################
# Segment XB3 ends
# Segment XB4 begins
###################### Segment XB4 begins as image of XB4 ######################
set XB4 [atomselect $m2 "serial 45714 to 45741"]
$XB4 set segname XB4
$XB4 set resid [list 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704]
$XB4 writepdb segtype_generic_XB4.pdb
segment XB4 {
    first none
    last none
    pdb segtype_generic_XB4.pdb
}
coordpdb segtype_generic_XB4.pdb XB4
############################### Segment XB4 ends ###############################
# Segment XB4 ends
# Segment XD1 begins
###################### Segment XD1 begins as image of XD1 ######################
set XD1 [atomselect $m2 "serial 52497 to 52524"]
$XD1 set segname XD1
$XD1 set resid [list 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701 701]
$XD1 writepdb segtype_generic_XD1.pdb
segment XD1 {
    first none
    last none
    pdb segtype_generic_XD1.pdb
}
coordpdb segtype_generic_XD1.pdb XD1
############################### Segment XD1 ends ###############################
# Segment XD1 ends
# Segment XD2 begins
###################### Segment XD2 begins as image of XD2 ######################
set XD2 [atomselect $m2 "serial 52525 to 52552"]
$XD2 set segname XD2
$XD2 set resid [list 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702 702]
$XD2 writepdb segtype_generic_XD2.pdb
segment XD2 {
    first none
    last none
    pdb segtype_generic_XD2.pdb
}
coordpdb segtype_generic_XD2.pdb XD2
############################### Segment XD2 ends ###############################
# Segment XD2 ends
# Segment XD3 begins
###################### Segment XD3 begins as image of XD3 ######################
set XD3 [atomselect $m2 "serial 52553 to 52580"]
$XD3 set segname XD3
$XD3 set resid [list 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703 703]
$XD3 writepdb segtype_generic_XD3.pdb
segment XD3 {
    first none
    last none
    pdb segtype_generic_XD3.pdb
}
coordpdb segtype_generic_XD3.pdb XD3
############################### Segment XD3 ends ###############################
# Segment XD3 ends
# Segment XD4 begins
###################### Segment XD4 begins as image of XD4 ######################
set XD4 [atomselect $m2 "serial 52581 to 52608"]
$XD4 set segname XD4
$XD4 set resid [list 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704 704]
$XD4 writepdb segtype_generic_XD4.pdb
segment XD4 {
    first none
    last none
    pdb segtype_generic_XD4.pdb
}
coordpdb segtype_generic_XD4.pdb XD4
############################### Segment XD4 ends ###############################
# Segment XD4 ends
# Segment XG1 begins
###################### Segment XG1 begins as image of XG1 ######################
set XG1 [atomselect $m2 "serial 35974 to 36000"]
$XG1 set segname XG1
$XG1 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG1 writepdb segtype_generic_XG1.pdb
segment XG1 {
    first none
    last none
    pdb segtype_generic_XG1.pdb
}
coordpdb segtype_generic_XG1.pdb XG1
############################### Segment XG1 ends ###############################
# Segment XG1 ends
# Segment XG10 begins
##################### Segment XG10 begins as image of XG10 #####################
set XG10 [atomselect $m2 "serial 36200 to 36220"]
$XG10 set segname XG10
$XG10 set resid [list 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5]
$XG10 writepdb segtype_generic_XG10.pdb
segment XG10 {
    first none
    last none
    pdb segtype_generic_XG10.pdb
}
coordpdb segtype_generic_XG10.pdb XG10
############################## Segment XG10 ends ###############################
# Segment XG10 ends
# Segment XG11 begins
##################### Segment XG11 begins as image of XG11 #####################
set XG11 [atomselect $m2 "serial 36221 to 36242"]
$XG11 set segname XG11
$XG11 set resid [list 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6]
$XG11 writepdb segtype_generic_XG11.pdb
segment XG11 {
    first none
    last none
    pdb segtype_generic_XG11.pdb
}
coordpdb segtype_generic_XG11.pdb XG11
############################## Segment XG11 ends ###############################
# Segment XG11 ends
# Segment XG12 begins
##################### Segment XG12 begins as image of XG12 #####################
set XG12 [atomselect $m2 "serial 36243 to 36264"]
$XG12 set segname XG12
$XG12 set resid [list 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7]
$XG12 writepdb segtype_generic_XG12.pdb
segment XG12 {
    first none
    last none
    pdb segtype_generic_XG12.pdb
}
coordpdb segtype_generic_XG12.pdb XG12
############################## Segment XG12 ends ###############################
# Segment XG12 ends
# Segment XG13 begins
##################### Segment XG13 begins as image of XG13 #####################
set XG13 [atomselect $m2 "serial 36265 to 36291"]
$XG13 set segname XG13
$XG13 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG13 writepdb segtype_generic_XG13.pdb
segment XG13 {
    first none
    last none
    pdb segtype_generic_XG13.pdb
}
coordpdb segtype_generic_XG13.pdb XG13
############################## Segment XG13 ends ###############################
# Segment XG13 ends
# Segment XG14 begins
##################### Segment XG14 begins as image of XG14 #####################
set XG14 [atomselect $m2 "serial 36292 to 36318"]
$XG14 set segname XG14
$XG14 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG14 writepdb segtype_generic_XG14.pdb
segment XG14 {
    first none
    last none
    pdb segtype_generic_XG14.pdb
}
coordpdb segtype_generic_XG14.pdb XG14
############################## Segment XG14 ends ###############################
# Segment XG14 ends
# Segment XG15 begins
##################### Segment XG15 begins as image of XG15 #####################
set XG15 [atomselect $m2 "serial 36319 to 36339"]
$XG15 set segname XG15
$XG15 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG15 writepdb segtype_generic_XG15.pdb
segment XG15 {
    first none
    last none
    pdb segtype_generic_XG15.pdb
}
coordpdb segtype_generic_XG15.pdb XG15
############################## Segment XG15 ends ###############################
# Segment XG15 ends
# Segment XG16 begins
##################### Segment XG16 begins as image of XG16 #####################
set XG16 [atomselect $m2 "serial 36340 to 36361"]
$XG16 set segname XG16
$XG16 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XG16 writepdb segtype_generic_XG16.pdb
segment XG16 {
    first none
    last none
    pdb segtype_generic_XG16.pdb
}
coordpdb segtype_generic_XG16.pdb XG16
############################## Segment XG16 ends ###############################
# Segment XG16 ends
# Segment XG17 begins
##################### Segment XG17 begins as image of XG17 #####################
set XG17 [atomselect $m2 "serial 36362 to 36388"]
$XG17 set segname XG17
$XG17 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG17 writepdb segtype_generic_XG17.pdb
segment XG17 {
    first none
    last none
    pdb segtype_generic_XG17.pdb
}
coordpdb segtype_generic_XG17.pdb XG17
############################## Segment XG17 ends ###############################
# Segment XG17 ends
# Segment XG18 begins
##################### Segment XG18 begins as image of XG18 #####################
set XG18 [atomselect $m2 "serial 36389 to 36416"]
$XG18 set segname XG18
$XG18 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG18 writepdb segtype_generic_XG18.pdb
segment XG18 {
    first none
    last none
    pdb segtype_generic_XG18.pdb
}
coordpdb segtype_generic_XG18.pdb XG18
############################## Segment XG18 ends ###############################
# Segment XG18 ends
# Segment XG19 begins
##################### Segment XG19 begins as image of XG19 #####################
set XG19 [atomselect $m2 "serial 36417 to 36443"]
$XG19 set segname XG19
$XG19 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG19 writepdb segtype_generic_XG19.pdb
segment XG19 {
    first none
    last none
    pdb segtype_generic_XG19.pdb
}
coordpdb segtype_generic_XG19.pdb XG19
############################## Segment XG19 ends ###############################
# Segment XG19 ends
# Segment XG2 begins
###################### Segment XG2 begins as image of XG2 ######################
set XG2 [atomselect $m2 "serial 36001 to 36028"]
$XG2 set segname XG2
$XG2 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG2 writepdb segtype_generic_XG2.pdb
segment XG2 {
    first none
    last none
    pdb segtype_generic_XG2.pdb
}
coordpdb segtype_generic_XG2.pdb XG2
############################### Segment XG2 ends ###############################
# Segment XG2 ends
# Segment XG20 begins
##################### Segment XG20 begins as image of XG20 #####################
set XG20 [atomselect $m2 "serial 36444 to 36471"]
$XG20 set segname XG20
$XG20 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG20 writepdb segtype_generic_XG20.pdb
segment XG20 {
    first none
    last none
    pdb segtype_generic_XG20.pdb
}
coordpdb segtype_generic_XG20.pdb XG20
############################## Segment XG20 ends ###############################
# Segment XG20 ends
# Segment XG21 begins
##################### Segment XG21 begins as image of XG21 #####################
set XG21 [atomselect $m2 "serial 36472 to 36498"]
$XG21 set segname XG21
$XG21 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG21 writepdb segtype_generic_XG21.pdb
segment XG21 {
    first none
    last none
    pdb segtype_generic_XG21.pdb
}
coordpdb segtype_generic_XG21.pdb XG21
############################## Segment XG21 ends ###############################
# Segment XG21 ends
# Segment XG22 begins
##################### Segment XG22 begins as image of XG22 #####################
set XG22 [atomselect $m2 "serial 36499 to 36525"]
$XG22 set segname XG22
$XG22 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG22 writepdb segtype_generic_XG22.pdb
segment XG22 {
    first none
    last none
    pdb segtype_generic_XG22.pdb
}
coordpdb segtype_generic_XG22.pdb XG22
############################## Segment XG22 ends ###############################
# Segment XG22 ends
# Segment XG23 begins
##################### Segment XG23 begins as image of XG23 #####################
set XG23 [atomselect $m2 "serial 36526 to 36546"]
$XG23 set segname XG23
$XG23 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG23 writepdb segtype_generic_XG23.pdb
segment XG23 {
    first none
    last none
    pdb segtype_generic_XG23.pdb
}
coordpdb segtype_generic_XG23.pdb XG23
############################## Segment XG23 ends ###############################
# Segment XG23 ends
# Segment XG24 begins
##################### Segment XG24 begins as image of XG24 #####################
set XG24 [atomselect $m2 "serial 36547 to 36568"]
$XG24 set segname XG24
$XG24 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XG24 writepdb segtype_generic_XG24.pdb
segment XG24 {
    first none
    last none
    pdb segtype_generic_XG24.pdb
}
coordpdb segtype_generic_XG24.pdb XG24
############################## Segment XG24 ends ###############################
# Segment XG24 ends
# Segment XG25 begins
##################### Segment XG25 begins as image of XG25 #####################
set XG25 [atomselect $m2 "serial 36569 to 36595"]
$XG25 set segname XG25
$XG25 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG25 writepdb segtype_generic_XG25.pdb
segment XG25 {
    first none
    last none
    pdb segtype_generic_XG25.pdb
}
coordpdb segtype_generic_XG25.pdb XG25
############################## Segment XG25 ends ###############################
# Segment XG25 ends
# Segment XG26 begins
##################### Segment XG26 begins as image of XG26 #####################
set XG26 [atomselect $m2 "serial 36596 to 36622"]
$XG26 set segname XG26
$XG26 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG26 writepdb segtype_generic_XG26.pdb
segment XG26 {
    first none
    last none
    pdb segtype_generic_XG26.pdb
}
coordpdb segtype_generic_XG26.pdb XG26
############################## Segment XG26 ends ###############################
# Segment XG26 ends
# Segment XG27 begins
##################### Segment XG27 begins as image of XG27 #####################
set XG27 [atomselect $m2 "serial 36623 to 36644"]
$XG27 set segname XG27
$XG27 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG27 writepdb segtype_generic_XG27.pdb
segment XG27 {
    first none
    last none
    pdb segtype_generic_XG27.pdb
}
coordpdb segtype_generic_XG27.pdb XG27
############################## Segment XG27 ends ###############################
# Segment XG27 ends
# Segment XG28 begins
##################### Segment XG28 begins as image of XG28 #####################
set XG28 [atomselect $m2 "serial 36645 to 36671"]
$XG28 set segname XG28
$XG28 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG28 writepdb segtype_generic_XG28.pdb
segment XG28 {
    first none
    last none
    pdb segtype_generic_XG28.pdb
}
coordpdb segtype_generic_XG28.pdb XG28
############################## Segment XG28 ends ###############################
# Segment XG28 ends
# Segment XG29 begins
##################### Segment XG29 begins as image of XG29 #####################
set XG29 [atomselect $m2 "serial 36672 to 36698"]
$XG29 set segname XG29
$XG29 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG29 writepdb segtype_generic_XG29.pdb
segment XG29 {
    first none
    last none
    pdb segtype_generic_XG29.pdb
}
coordpdb segtype_generic_XG29.pdb XG29
############################## Segment XG29 ends ###############################
# Segment XG29 ends
# Segment XG3 begins
###################### Segment XG3 begins as image of XG3 ######################
set XG3 [atomselect $m2 "serial 36029 to 36055"]
$XG3 set segname XG3
$XG3 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG3 writepdb segtype_generic_XG3.pdb
segment XG3 {
    first none
    last none
    pdb segtype_generic_XG3.pdb
}
coordpdb segtype_generic_XG3.pdb XG3
############################### Segment XG3 ends ###############################
# Segment XG3 ends
# Segment XG30 begins
##################### Segment XG30 begins as image of XG30 #####################
set XG30 [atomselect $m2 "serial 36699 to 36720"]
$XG30 set segname XG30
$XG30 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG30 writepdb segtype_generic_XG30.pdb
segment XG30 {
    first none
    last none
    pdb segtype_generic_XG30.pdb
}
coordpdb segtype_generic_XG30.pdb XG30
############################## Segment XG30 ends ###############################
# Segment XG30 ends
# Segment XG31 begins
##################### Segment XG31 begins as image of XG31 #####################
set XG31 [atomselect $m2 "serial 36721 to 36747"]
$XG31 set segname XG31
$XG31 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG31 writepdb segtype_generic_XG31.pdb
segment XG31 {
    first none
    last none
    pdb segtype_generic_XG31.pdb
}
coordpdb segtype_generic_XG31.pdb XG31
############################## Segment XG31 ends ###############################
# Segment XG31 ends
# Segment XG32 begins
##################### Segment XG32 begins as image of XG32 #####################
set XG32 [atomselect $m2 "serial 36748 to 36775"]
$XG32 set segname XG32
$XG32 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG32 writepdb segtype_generic_XG32.pdb
segment XG32 {
    first none
    last none
    pdb segtype_generic_XG32.pdb
}
coordpdb segtype_generic_XG32.pdb XG32
############################## Segment XG32 ends ###############################
# Segment XG32 ends
# Segment XG33 begins
##################### Segment XG33 begins as image of XG33 #####################
set XG33 [atomselect $m2 "serial 36776 to 36802"]
$XG33 set segname XG33
$XG33 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG33 writepdb segtype_generic_XG33.pdb
segment XG33 {
    first none
    last none
    pdb segtype_generic_XG33.pdb
}
coordpdb segtype_generic_XG33.pdb XG33
############################## Segment XG33 ends ###############################
# Segment XG33 ends
# Segment XG34 begins
##################### Segment XG34 begins as image of XG34 #####################
set XG34 [atomselect $m2 "serial 36803 to 36829"]
$XG34 set segname XG34
$XG34 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG34 writepdb segtype_generic_XG34.pdb
segment XG34 {
    first none
    last none
    pdb segtype_generic_XG34.pdb
}
coordpdb segtype_generic_XG34.pdb XG34
############################## Segment XG34 ends ###############################
# Segment XG34 ends
# Segment XG35 begins
##################### Segment XG35 begins as image of XG35 #####################
set XG35 [atomselect $m2 "serial 36830 to 36851"]
$XG35 set segname XG35
$XG35 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG35 writepdb segtype_generic_XG35.pdb
segment XG35 {
    first none
    last none
    pdb segtype_generic_XG35.pdb
}
coordpdb segtype_generic_XG35.pdb XG35
############################## Segment XG35 ends ###############################
# Segment XG35 ends
# Segment XG36 begins
##################### Segment XG36 begins as image of XG36 #####################
set XG36 [atomselect $m2 "serial 42119 to 42146"]
$XG36 set segname XG36
$XG36 set resid [list 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601]
$XG36 writepdb segtype_generic_XG36.pdb
segment XG36 {
    first none
    last none
    pdb segtype_generic_XG36.pdb
}
coordpdb segtype_generic_XG36.pdb XG36
############################## Segment XG36 ends ###############################
# Segment XG36 ends
# Segment XG37 begins
##################### Segment XG37 begins as image of XG37 #####################
set XG37 [atomselect $m2 "serial 42147 to 42174"]
$XG37 set segname XG37
$XG37 set resid [list 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602]
$XG37 writepdb segtype_generic_XG37.pdb
segment XG37 {
    first none
    last none
    pdb segtype_generic_XG37.pdb
}
coordpdb segtype_generic_XG37.pdb XG37
############################## Segment XG37 ends ###############################
# Segment XG37 ends
# Segment XG38 begins
##################### Segment XG38 begins as image of XG38 #####################
set XG38 [atomselect $m2 "serial 42175 to 42202"]
$XG38 set segname XG38
$XG38 set resid [list 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603]
$XG38 writepdb segtype_generic_XG38.pdb
segment XG38 {
    first none
    last none
    pdb segtype_generic_XG38.pdb
}
coordpdb segtype_generic_XG38.pdb XG38
############################## Segment XG38 ends ###############################
# Segment XG38 ends
# Segment XG39 begins
##################### Segment XG39 begins as image of XG39 #####################
set XG39 [atomselect $m2 "serial 42203 to 42230"]
$XG39 set segname XG39
$XG39 set resid [list 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604]
$XG39 writepdb segtype_generic_XG39.pdb
segment XG39 {
    first none
    last none
    pdb segtype_generic_XG39.pdb
}
coordpdb segtype_generic_XG39.pdb XG39
############################## Segment XG39 ends ###############################
# Segment XG39 ends
# Segment XG4 begins
###################### Segment XG4 begins as image of XG4 ######################
set XG4 [atomselect $m2 "serial 36056 to 36082"]
$XG4 set segname XG4
$XG4 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG4 writepdb segtype_generic_XG4.pdb
segment XG4 {
    first none
    last none
    pdb segtype_generic_XG4.pdb
}
coordpdb segtype_generic_XG4.pdb XG4
############################### Segment XG4 ends ###############################
# Segment XG4 ends
# Segment XG5 begins
###################### Segment XG5 begins as image of XG5 ######################
set XG5 [atomselect $m2 "serial 36083 to 36104"]
$XG5 set segname XG5
$XG5 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG5 writepdb segtype_generic_XG5.pdb
segment XG5 {
    first none
    last none
    pdb segtype_generic_XG5.pdb
}
coordpdb segtype_generic_XG5.pdb XG5
############################### Segment XG5 ends ###############################
# Segment XG5 ends
# Segment XG6 begins
###################### Segment XG6 begins as image of XG6 ######################
set XG6 [atomselect $m2 "serial 36105 to 36131"]
$XG6 set segname XG6
$XG6 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XG6 writepdb segtype_generic_XG6.pdb
segment XG6 {
    first none
    last none
    pdb segtype_generic_XG6.pdb
}
coordpdb segtype_generic_XG6.pdb XG6
############################### Segment XG6 ends ###############################
# Segment XG6 ends
# Segment XG7 begins
###################### Segment XG7 begins as image of XG7 ######################
set XG7 [atomselect $m2 "serial 36132 to 36158"]
$XG7 set segname XG7
$XG7 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XG7 writepdb segtype_generic_XG7.pdb
segment XG7 {
    first none
    last none
    pdb segtype_generic_XG7.pdb
}
coordpdb segtype_generic_XG7.pdb XG7
############################### Segment XG7 ends ###############################
# Segment XG7 ends
# Segment XG8 begins
###################### Segment XG8 begins as image of XG8 ######################
set XG8 [atomselect $m2 "serial 36159 to 36178"]
$XG8 set segname XG8
$XG8 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XG8 writepdb segtype_generic_XG8.pdb
segment XG8 {
    first none
    last none
    pdb segtype_generic_XG8.pdb
}
coordpdb segtype_generic_XG8.pdb XG8
############################### Segment XG8 ends ###############################
# Segment XG8 ends
# Segment XG9 begins
###################### Segment XG9 begins as image of XG9 ######################
set XG9 [atomselect $m2 "serial 36179 to 36199"]
$XG9 set segname XG9
$XG9 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XG9 writepdb segtype_generic_XG9.pdb
segment XG9 {
    first none
    last none
    pdb segtype_generic_XG9.pdb
}
coordpdb segtype_generic_XG9.pdb XG9
############################### Segment XG9 ends ###############################
# Segment XG9 ends
# Segment XI1 begins
###################### Segment XI1 begins as image of XI1 ######################
set XI1 [atomselect $m2 "serial 36852 to 36878"]
$XI1 set segname XI1
$XI1 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI1 writepdb segtype_generic_XI1.pdb
segment XI1 {
    first none
    last none
    pdb segtype_generic_XI1.pdb
}
coordpdb segtype_generic_XI1.pdb XI1
############################### Segment XI1 ends ###############################
# Segment XI1 ends
# Segment XI10 begins
##################### Segment XI10 begins as image of XI10 #####################
set XI10 [atomselect $m2 "serial 37078 to 37098"]
$XI10 set segname XI10
$XI10 set resid [list 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5]
$XI10 writepdb segtype_generic_XI10.pdb
segment XI10 {
    first none
    last none
    pdb segtype_generic_XI10.pdb
}
coordpdb segtype_generic_XI10.pdb XI10
############################## Segment XI10 ends ###############################
# Segment XI10 ends
# Segment XI11 begins
##################### Segment XI11 begins as image of XI11 #####################
set XI11 [atomselect $m2 "serial 37099 to 37120"]
$XI11 set segname XI11
$XI11 set resid [list 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6]
$XI11 writepdb segtype_generic_XI11.pdb
segment XI11 {
    first none
    last none
    pdb segtype_generic_XI11.pdb
}
coordpdb segtype_generic_XI11.pdb XI11
############################## Segment XI11 ends ###############################
# Segment XI11 ends
# Segment XI12 begins
##################### Segment XI12 begins as image of XI12 #####################
set XI12 [atomselect $m2 "serial 37121 to 37142"]
$XI12 set segname XI12
$XI12 set resid [list 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7]
$XI12 writepdb segtype_generic_XI12.pdb
segment XI12 {
    first none
    last none
    pdb segtype_generic_XI12.pdb
}
coordpdb segtype_generic_XI12.pdb XI12
############################## Segment XI12 ends ###############################
# Segment XI12 ends
# Segment XI13 begins
##################### Segment XI13 begins as image of XI13 #####################
set XI13 [atomselect $m2 "serial 37143 to 37169"]
$XI13 set segname XI13
$XI13 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI13 writepdb segtype_generic_XI13.pdb
segment XI13 {
    first none
    last none
    pdb segtype_generic_XI13.pdb
}
coordpdb segtype_generic_XI13.pdb XI13
############################## Segment XI13 ends ###############################
# Segment XI13 ends
# Segment XI14 begins
##################### Segment XI14 begins as image of XI14 #####################
set XI14 [atomselect $m2 "serial 37170 to 37196"]
$XI14 set segname XI14
$XI14 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI14 writepdb segtype_generic_XI14.pdb
segment XI14 {
    first none
    last none
    pdb segtype_generic_XI14.pdb
}
coordpdb segtype_generic_XI14.pdb XI14
############################## Segment XI14 ends ###############################
# Segment XI14 ends
# Segment XI15 begins
##################### Segment XI15 begins as image of XI15 #####################
set XI15 [atomselect $m2 "serial 37197 to 37217"]
$XI15 set segname XI15
$XI15 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI15 writepdb segtype_generic_XI15.pdb
segment XI15 {
    first none
    last none
    pdb segtype_generic_XI15.pdb
}
coordpdb segtype_generic_XI15.pdb XI15
############################## Segment XI15 ends ###############################
# Segment XI15 ends
# Segment XI16 begins
##################### Segment XI16 begins as image of XI16 #####################
set XI16 [atomselect $m2 "serial 37218 to 37239"]
$XI16 set segname XI16
$XI16 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XI16 writepdb segtype_generic_XI16.pdb
segment XI16 {
    first none
    last none
    pdb segtype_generic_XI16.pdb
}
coordpdb segtype_generic_XI16.pdb XI16
############################## Segment XI16 ends ###############################
# Segment XI16 ends
# Segment XI17 begins
##################### Segment XI17 begins as image of XI17 #####################
set XI17 [atomselect $m2 "serial 37240 to 37266"]
$XI17 set segname XI17
$XI17 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI17 writepdb segtype_generic_XI17.pdb
segment XI17 {
    first none
    last none
    pdb segtype_generic_XI17.pdb
}
coordpdb segtype_generic_XI17.pdb XI17
############################## Segment XI17 ends ###############################
# Segment XI17 ends
# Segment XI18 begins
##################### Segment XI18 begins as image of XI18 #####################
set XI18 [atomselect $m2 "serial 37267 to 37294"]
$XI18 set segname XI18
$XI18 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI18 writepdb segtype_generic_XI18.pdb
segment XI18 {
    first none
    last none
    pdb segtype_generic_XI18.pdb
}
coordpdb segtype_generic_XI18.pdb XI18
############################## Segment XI18 ends ###############################
# Segment XI18 ends
# Segment XI19 begins
##################### Segment XI19 begins as image of XI19 #####################
set XI19 [atomselect $m2 "serial 37295 to 37321"]
$XI19 set segname XI19
$XI19 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI19 writepdb segtype_generic_XI19.pdb
segment XI19 {
    first none
    last none
    pdb segtype_generic_XI19.pdb
}
coordpdb segtype_generic_XI19.pdb XI19
############################## Segment XI19 ends ###############################
# Segment XI19 ends
# Segment XI2 begins
###################### Segment XI2 begins as image of XI2 ######################
set XI2 [atomselect $m2 "serial 36879 to 36906"]
$XI2 set segname XI2
$XI2 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI2 writepdb segtype_generic_XI2.pdb
segment XI2 {
    first none
    last none
    pdb segtype_generic_XI2.pdb
}
coordpdb segtype_generic_XI2.pdb XI2
############################### Segment XI2 ends ###############################
# Segment XI2 ends
# Segment XI20 begins
##################### Segment XI20 begins as image of XI20 #####################
set XI20 [atomselect $m2 "serial 37322 to 37349"]
$XI20 set segname XI20
$XI20 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI20 writepdb segtype_generic_XI20.pdb
segment XI20 {
    first none
    last none
    pdb segtype_generic_XI20.pdb
}
coordpdb segtype_generic_XI20.pdb XI20
############################## Segment XI20 ends ###############################
# Segment XI20 ends
# Segment XI21 begins
##################### Segment XI21 begins as image of XI21 #####################
set XI21 [atomselect $m2 "serial 37350 to 37376"]
$XI21 set segname XI21
$XI21 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI21 writepdb segtype_generic_XI21.pdb
segment XI21 {
    first none
    last none
    pdb segtype_generic_XI21.pdb
}
coordpdb segtype_generic_XI21.pdb XI21
############################## Segment XI21 ends ###############################
# Segment XI21 ends
# Segment XI22 begins
##################### Segment XI22 begins as image of XI22 #####################
set XI22 [atomselect $m2 "serial 37377 to 37403"]
$XI22 set segname XI22
$XI22 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI22 writepdb segtype_generic_XI22.pdb
segment XI22 {
    first none
    last none
    pdb segtype_generic_XI22.pdb
}
coordpdb segtype_generic_XI22.pdb XI22
############################## Segment XI22 ends ###############################
# Segment XI22 ends
# Segment XI23 begins
##################### Segment XI23 begins as image of XI23 #####################
set XI23 [atomselect $m2 "serial 37404 to 37424"]
$XI23 set segname XI23
$XI23 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI23 writepdb segtype_generic_XI23.pdb
segment XI23 {
    first none
    last none
    pdb segtype_generic_XI23.pdb
}
coordpdb segtype_generic_XI23.pdb XI23
############################## Segment XI23 ends ###############################
# Segment XI23 ends
# Segment XI24 begins
##################### Segment XI24 begins as image of XI24 #####################
set XI24 [atomselect $m2 "serial 37425 to 37446"]
$XI24 set segname XI24
$XI24 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XI24 writepdb segtype_generic_XI24.pdb
segment XI24 {
    first none
    last none
    pdb segtype_generic_XI24.pdb
}
coordpdb segtype_generic_XI24.pdb XI24
############################## Segment XI24 ends ###############################
# Segment XI24 ends
# Segment XI25 begins
##################### Segment XI25 begins as image of XI25 #####################
set XI25 [atomselect $m2 "serial 37447 to 37473"]
$XI25 set segname XI25
$XI25 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI25 writepdb segtype_generic_XI25.pdb
segment XI25 {
    first none
    last none
    pdb segtype_generic_XI25.pdb
}
coordpdb segtype_generic_XI25.pdb XI25
############################## Segment XI25 ends ###############################
# Segment XI25 ends
# Segment XI26 begins
##################### Segment XI26 begins as image of XI26 #####################
set XI26 [atomselect $m2 "serial 37474 to 37500"]
$XI26 set segname XI26
$XI26 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI26 writepdb segtype_generic_XI26.pdb
segment XI26 {
    first none
    last none
    pdb segtype_generic_XI26.pdb
}
coordpdb segtype_generic_XI26.pdb XI26
############################## Segment XI26 ends ###############################
# Segment XI26 ends
# Segment XI27 begins
##################### Segment XI27 begins as image of XI27 #####################
set XI27 [atomselect $m2 "serial 37501 to 37522"]
$XI27 set segname XI27
$XI27 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI27 writepdb segtype_generic_XI27.pdb
segment XI27 {
    first none
    last none
    pdb segtype_generic_XI27.pdb
}
coordpdb segtype_generic_XI27.pdb XI27
############################## Segment XI27 ends ###############################
# Segment XI27 ends
# Segment XI28 begins
##################### Segment XI28 begins as image of XI28 #####################
set XI28 [atomselect $m2 "serial 37523 to 37549"]
$XI28 set segname XI28
$XI28 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI28 writepdb segtype_generic_XI28.pdb
segment XI28 {
    first none
    last none
    pdb segtype_generic_XI28.pdb
}
coordpdb segtype_generic_XI28.pdb XI28
############################## Segment XI28 ends ###############################
# Segment XI28 ends
# Segment XI29 begins
##################### Segment XI29 begins as image of XI29 #####################
set XI29 [atomselect $m2 "serial 37550 to 37576"]
$XI29 set segname XI29
$XI29 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI29 writepdb segtype_generic_XI29.pdb
segment XI29 {
    first none
    last none
    pdb segtype_generic_XI29.pdb
}
coordpdb segtype_generic_XI29.pdb XI29
############################## Segment XI29 ends ###############################
# Segment XI29 ends
# Segment XI3 begins
###################### Segment XI3 begins as image of XI3 ######################
set XI3 [atomselect $m2 "serial 36907 to 36933"]
$XI3 set segname XI3
$XI3 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI3 writepdb segtype_generic_XI3.pdb
segment XI3 {
    first none
    last none
    pdb segtype_generic_XI3.pdb
}
coordpdb segtype_generic_XI3.pdb XI3
############################### Segment XI3 ends ###############################
# Segment XI3 ends
# Segment XI30 begins
##################### Segment XI30 begins as image of XI30 #####################
set XI30 [atomselect $m2 "serial 37577 to 37598"]
$XI30 set segname XI30
$XI30 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI30 writepdb segtype_generic_XI30.pdb
segment XI30 {
    first none
    last none
    pdb segtype_generic_XI30.pdb
}
coordpdb segtype_generic_XI30.pdb XI30
############################## Segment XI30 ends ###############################
# Segment XI30 ends
# Segment XI31 begins
##################### Segment XI31 begins as image of XI31 #####################
set XI31 [atomselect $m2 "serial 37599 to 37625"]
$XI31 set segname XI31
$XI31 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI31 writepdb segtype_generic_XI31.pdb
segment XI31 {
    first none
    last none
    pdb segtype_generic_XI31.pdb
}
coordpdb segtype_generic_XI31.pdb XI31
############################## Segment XI31 ends ###############################
# Segment XI31 ends
# Segment XI32 begins
##################### Segment XI32 begins as image of XI32 #####################
set XI32 [atomselect $m2 "serial 37626 to 37653"]
$XI32 set segname XI32
$XI32 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI32 writepdb segtype_generic_XI32.pdb
segment XI32 {
    first none
    last none
    pdb segtype_generic_XI32.pdb
}
coordpdb segtype_generic_XI32.pdb XI32
############################## Segment XI32 ends ###############################
# Segment XI32 ends
# Segment XI33 begins
##################### Segment XI33 begins as image of XI33 #####################
set XI33 [atomselect $m2 "serial 37654 to 37680"]
$XI33 set segname XI33
$XI33 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI33 writepdb segtype_generic_XI33.pdb
segment XI33 {
    first none
    last none
    pdb segtype_generic_XI33.pdb
}
coordpdb segtype_generic_XI33.pdb XI33
############################## Segment XI33 ends ###############################
# Segment XI33 ends
# Segment XI34 begins
##################### Segment XI34 begins as image of XI34 #####################
set XI34 [atomselect $m2 "serial 37681 to 37707"]
$XI34 set segname XI34
$XI34 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI34 writepdb segtype_generic_XI34.pdb
segment XI34 {
    first none
    last none
    pdb segtype_generic_XI34.pdb
}
coordpdb segtype_generic_XI34.pdb XI34
############################## Segment XI34 ends ###############################
# Segment XI34 ends
# Segment XI35 begins
##################### Segment XI35 begins as image of XI35 #####################
set XI35 [atomselect $m2 "serial 37708 to 37729"]
$XI35 set segname XI35
$XI35 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI35 writepdb segtype_generic_XI35.pdb
segment XI35 {
    first none
    last none
    pdb segtype_generic_XI35.pdb
}
coordpdb segtype_generic_XI35.pdb XI35
############################## Segment XI35 ends ###############################
# Segment XI35 ends
# Segment XI36 begins
##################### Segment XI36 begins as image of XI36 #####################
set XI36 [atomselect $m2 "serial 48986 to 49013"]
$XI36 set segname XI36
$XI36 set resid [list 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601]
$XI36 writepdb segtype_generic_XI36.pdb
segment XI36 {
    first none
    last none
    pdb segtype_generic_XI36.pdb
}
coordpdb segtype_generic_XI36.pdb XI36
############################## Segment XI36 ends ###############################
# Segment XI36 ends
# Segment XI37 begins
##################### Segment XI37 begins as image of XI37 #####################
set XI37 [atomselect $m2 "serial 49014 to 49041"]
$XI37 set segname XI37
$XI37 set resid [list 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602]
$XI37 writepdb segtype_generic_XI37.pdb
segment XI37 {
    first none
    last none
    pdb segtype_generic_XI37.pdb
}
coordpdb segtype_generic_XI37.pdb XI37
############################## Segment XI37 ends ###############################
# Segment XI37 ends
# Segment XI38 begins
##################### Segment XI38 begins as image of XI38 #####################
set XI38 [atomselect $m2 "serial 49042 to 49069"]
$XI38 set segname XI38
$XI38 set resid [list 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603]
$XI38 writepdb segtype_generic_XI38.pdb
segment XI38 {
    first none
    last none
    pdb segtype_generic_XI38.pdb
}
coordpdb segtype_generic_XI38.pdb XI38
############################## Segment XI38 ends ###############################
# Segment XI38 ends
# Segment XI39 begins
##################### Segment XI39 begins as image of XI39 #####################
set XI39 [atomselect $m2 "serial 49070 to 49097"]
$XI39 set segname XI39
$XI39 set resid [list 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604]
$XI39 writepdb segtype_generic_XI39.pdb
segment XI39 {
    first none
    last none
    pdb segtype_generic_XI39.pdb
}
coordpdb segtype_generic_XI39.pdb XI39
############################## Segment XI39 ends ###############################
# Segment XI39 ends
# Segment XI4 begins
###################### Segment XI4 begins as image of XI4 ######################
set XI4 [atomselect $m2 "serial 36934 to 36960"]
$XI4 set segname XI4
$XI4 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI4 writepdb segtype_generic_XI4.pdb
segment XI4 {
    first none
    last none
    pdb segtype_generic_XI4.pdb
}
coordpdb segtype_generic_XI4.pdb XI4
############################### Segment XI4 ends ###############################
# Segment XI4 ends
# Segment XI5 begins
###################### Segment XI5 begins as image of XI5 ######################
set XI5 [atomselect $m2 "serial 36961 to 36982"]
$XI5 set segname XI5
$XI5 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI5 writepdb segtype_generic_XI5.pdb
segment XI5 {
    first none
    last none
    pdb segtype_generic_XI5.pdb
}
coordpdb segtype_generic_XI5.pdb XI5
############################### Segment XI5 ends ###############################
# Segment XI5 ends
# Segment XI6 begins
###################### Segment XI6 begins as image of XI6 ######################
set XI6 [atomselect $m2 "serial 36983 to 37009"]
$XI6 set segname XI6
$XI6 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XI6 writepdb segtype_generic_XI6.pdb
segment XI6 {
    first none
    last none
    pdb segtype_generic_XI6.pdb
}
coordpdb segtype_generic_XI6.pdb XI6
############################### Segment XI6 ends ###############################
# Segment XI6 ends
# Segment XI7 begins
###################### Segment XI7 begins as image of XI7 ######################
set XI7 [atomselect $m2 "serial 37010 to 37036"]
$XI7 set segname XI7
$XI7 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XI7 writepdb segtype_generic_XI7.pdb
segment XI7 {
    first none
    last none
    pdb segtype_generic_XI7.pdb
}
coordpdb segtype_generic_XI7.pdb XI7
############################### Segment XI7 ends ###############################
# Segment XI7 ends
# Segment XI8 begins
###################### Segment XI8 begins as image of XI8 ######################
set XI8 [atomselect $m2 "serial 37037 to 37056"]
$XI8 set segname XI8
$XI8 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XI8 writepdb segtype_generic_XI8.pdb
segment XI8 {
    first none
    last none
    pdb segtype_generic_XI8.pdb
}
coordpdb segtype_generic_XI8.pdb XI8
############################### Segment XI8 ends ###############################
# Segment XI8 ends
# Segment XI9 begins
###################### Segment XI9 begins as image of XI9 ######################
set XI9 [atomselect $m2 "serial 37057 to 37077"]
$XI9 set segname XI9
$XI9 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XI9 writepdb segtype_generic_XI9.pdb
segment XI9 {
    first none
    last none
    pdb segtype_generic_XI9.pdb
}
coordpdb segtype_generic_XI9.pdb XI9
############################### Segment XI9 ends ###############################
# Segment XI9 ends
# Segment XJ1 begins
###################### Segment XJ1 begins as image of XJ1 ######################
set XJ1 [atomselect $m2 "serial 37730 to 37756"]
$XJ1 set segname XJ1
$XJ1 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ1 writepdb segtype_generic_XJ1.pdb
segment XJ1 {
    first none
    last none
    pdb segtype_generic_XJ1.pdb
}
coordpdb segtype_generic_XJ1.pdb XJ1
############################### Segment XJ1 ends ###############################
# Segment XJ1 ends
# Segment XJ10 begins
##################### Segment XJ10 begins as image of XJ10 #####################
set XJ10 [atomselect $m2 "serial 37956 to 37976"]
$XJ10 set segname XJ10
$XJ10 set resid [list 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5]
$XJ10 writepdb segtype_generic_XJ10.pdb
segment XJ10 {
    first none
    last none
    pdb segtype_generic_XJ10.pdb
}
coordpdb segtype_generic_XJ10.pdb XJ10
############################## Segment XJ10 ends ###############################
# Segment XJ10 ends
# Segment XJ11 begins
##################### Segment XJ11 begins as image of XJ11 #####################
set XJ11 [atomselect $m2 "serial 37977 to 37998"]
$XJ11 set segname XJ11
$XJ11 set resid [list 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6]
$XJ11 writepdb segtype_generic_XJ11.pdb
segment XJ11 {
    first none
    last none
    pdb segtype_generic_XJ11.pdb
}
coordpdb segtype_generic_XJ11.pdb XJ11
############################## Segment XJ11 ends ###############################
# Segment XJ11 ends
# Segment XJ12 begins
##################### Segment XJ12 begins as image of XJ12 #####################
set XJ12 [atomselect $m2 "serial 37999 to 38020"]
$XJ12 set segname XJ12
$XJ12 set resid [list 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7 7]
$XJ12 writepdb segtype_generic_XJ12.pdb
segment XJ12 {
    first none
    last none
    pdb segtype_generic_XJ12.pdb
}
coordpdb segtype_generic_XJ12.pdb XJ12
############################## Segment XJ12 ends ###############################
# Segment XJ12 ends
# Segment XJ13 begins
##################### Segment XJ13 begins as image of XJ13 #####################
set XJ13 [atomselect $m2 "serial 38021 to 38047"]
$XJ13 set segname XJ13
$XJ13 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ13 writepdb segtype_generic_XJ13.pdb
segment XJ13 {
    first none
    last none
    pdb segtype_generic_XJ13.pdb
}
coordpdb segtype_generic_XJ13.pdb XJ13
############################## Segment XJ13 ends ###############################
# Segment XJ13 ends
# Segment XJ14 begins
##################### Segment XJ14 begins as image of XJ14 #####################
set XJ14 [atomselect $m2 "serial 38048 to 38074"]
$XJ14 set segname XJ14
$XJ14 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ14 writepdb segtype_generic_XJ14.pdb
segment XJ14 {
    first none
    last none
    pdb segtype_generic_XJ14.pdb
}
coordpdb segtype_generic_XJ14.pdb XJ14
############################## Segment XJ14 ends ###############################
# Segment XJ14 ends
# Segment XJ15 begins
##################### Segment XJ15 begins as image of XJ15 #####################
set XJ15 [atomselect $m2 "serial 38075 to 38095"]
$XJ15 set segname XJ15
$XJ15 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ15 writepdb segtype_generic_XJ15.pdb
segment XJ15 {
    first none
    last none
    pdb segtype_generic_XJ15.pdb
}
coordpdb segtype_generic_XJ15.pdb XJ15
############################## Segment XJ15 ends ###############################
# Segment XJ15 ends
# Segment XJ16 begins
##################### Segment XJ16 begins as image of XJ16 #####################
set XJ16 [atomselect $m2 "serial 38096 to 38117"]
$XJ16 set segname XJ16
$XJ16 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XJ16 writepdb segtype_generic_XJ16.pdb
segment XJ16 {
    first none
    last none
    pdb segtype_generic_XJ16.pdb
}
coordpdb segtype_generic_XJ16.pdb XJ16
############################## Segment XJ16 ends ###############################
# Segment XJ16 ends
# Segment XJ17 begins
##################### Segment XJ17 begins as image of XJ17 #####################
set XJ17 [atomselect $m2 "serial 38118 to 38144"]
$XJ17 set segname XJ17
$XJ17 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ17 writepdb segtype_generic_XJ17.pdb
segment XJ17 {
    first none
    last none
    pdb segtype_generic_XJ17.pdb
}
coordpdb segtype_generic_XJ17.pdb XJ17
############################## Segment XJ17 ends ###############################
# Segment XJ17 ends
# Segment XJ18 begins
##################### Segment XJ18 begins as image of XJ18 #####################
set XJ18 [atomselect $m2 "serial 38145 to 38172"]
$XJ18 set segname XJ18
$XJ18 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ18 writepdb segtype_generic_XJ18.pdb
segment XJ18 {
    first none
    last none
    pdb segtype_generic_XJ18.pdb
}
coordpdb segtype_generic_XJ18.pdb XJ18
############################## Segment XJ18 ends ###############################
# Segment XJ18 ends
# Segment XJ19 begins
##################### Segment XJ19 begins as image of XJ19 #####################
set XJ19 [atomselect $m2 "serial 38173 to 38199"]
$XJ19 set segname XJ19
$XJ19 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ19 writepdb segtype_generic_XJ19.pdb
segment XJ19 {
    first none
    last none
    pdb segtype_generic_XJ19.pdb
}
coordpdb segtype_generic_XJ19.pdb XJ19
############################## Segment XJ19 ends ###############################
# Segment XJ19 ends
# Segment XJ2 begins
###################### Segment XJ2 begins as image of XJ2 ######################
set XJ2 [atomselect $m2 "serial 37757 to 37784"]
$XJ2 set segname XJ2
$XJ2 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ2 writepdb segtype_generic_XJ2.pdb
segment XJ2 {
    first none
    last none
    pdb segtype_generic_XJ2.pdb
}
coordpdb segtype_generic_XJ2.pdb XJ2
############################### Segment XJ2 ends ###############################
# Segment XJ2 ends
# Segment XJ20 begins
##################### Segment XJ20 begins as image of XJ20 #####################
set XJ20 [atomselect $m2 "serial 38200 to 38227"]
$XJ20 set segname XJ20
$XJ20 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ20 writepdb segtype_generic_XJ20.pdb
segment XJ20 {
    first none
    last none
    pdb segtype_generic_XJ20.pdb
}
coordpdb segtype_generic_XJ20.pdb XJ20
############################## Segment XJ20 ends ###############################
# Segment XJ20 ends
# Segment XJ21 begins
##################### Segment XJ21 begins as image of XJ21 #####################
set XJ21 [atomselect $m2 "serial 38228 to 38254"]
$XJ21 set segname XJ21
$XJ21 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ21 writepdb segtype_generic_XJ21.pdb
segment XJ21 {
    first none
    last none
    pdb segtype_generic_XJ21.pdb
}
coordpdb segtype_generic_XJ21.pdb XJ21
############################## Segment XJ21 ends ###############################
# Segment XJ21 ends
# Segment XJ22 begins
##################### Segment XJ22 begins as image of XJ22 #####################
set XJ22 [atomselect $m2 "serial 38255 to 38281"]
$XJ22 set segname XJ22
$XJ22 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ22 writepdb segtype_generic_XJ22.pdb
segment XJ22 {
    first none
    last none
    pdb segtype_generic_XJ22.pdb
}
coordpdb segtype_generic_XJ22.pdb XJ22
############################## Segment XJ22 ends ###############################
# Segment XJ22 ends
# Segment XJ23 begins
##################### Segment XJ23 begins as image of XJ23 #####################
set XJ23 [atomselect $m2 "serial 38282 to 38302"]
$XJ23 set segname XJ23
$XJ23 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ23 writepdb segtype_generic_XJ23.pdb
segment XJ23 {
    first none
    last none
    pdb segtype_generic_XJ23.pdb
}
coordpdb segtype_generic_XJ23.pdb XJ23
############################## Segment XJ23 ends ###############################
# Segment XJ23 ends
# Segment XJ24 begins
##################### Segment XJ24 begins as image of XJ24 #####################
set XJ24 [atomselect $m2 "serial 38303 to 38324"]
$XJ24 set segname XJ24
$XJ24 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XJ24 writepdb segtype_generic_XJ24.pdb
segment XJ24 {
    first none
    last none
    pdb segtype_generic_XJ24.pdb
}
coordpdb segtype_generic_XJ24.pdb XJ24
############################## Segment XJ24 ends ###############################
# Segment XJ24 ends
# Segment XJ25 begins
##################### Segment XJ25 begins as image of XJ25 #####################
set XJ25 [atomselect $m2 "serial 38325 to 38351"]
$XJ25 set segname XJ25
$XJ25 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ25 writepdb segtype_generic_XJ25.pdb
segment XJ25 {
    first none
    last none
    pdb segtype_generic_XJ25.pdb
}
coordpdb segtype_generic_XJ25.pdb XJ25
############################## Segment XJ25 ends ###############################
# Segment XJ25 ends
# Segment XJ26 begins
##################### Segment XJ26 begins as image of XJ26 #####################
set XJ26 [atomselect $m2 "serial 38352 to 38378"]
$XJ26 set segname XJ26
$XJ26 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ26 writepdb segtype_generic_XJ26.pdb
segment XJ26 {
    first none
    last none
    pdb segtype_generic_XJ26.pdb
}
coordpdb segtype_generic_XJ26.pdb XJ26
############################## Segment XJ26 ends ###############################
# Segment XJ26 ends
# Segment XJ27 begins
##################### Segment XJ27 begins as image of XJ27 #####################
set XJ27 [atomselect $m2 "serial 38379 to 38400"]
$XJ27 set segname XJ27
$XJ27 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ27 writepdb segtype_generic_XJ27.pdb
segment XJ27 {
    first none
    last none
    pdb segtype_generic_XJ27.pdb
}
coordpdb segtype_generic_XJ27.pdb XJ27
############################## Segment XJ27 ends ###############################
# Segment XJ27 ends
# Segment XJ28 begins
##################### Segment XJ28 begins as image of XJ28 #####################
set XJ28 [atomselect $m2 "serial 38401 to 38427"]
$XJ28 set segname XJ28
$XJ28 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ28 writepdb segtype_generic_XJ28.pdb
segment XJ28 {
    first none
    last none
    pdb segtype_generic_XJ28.pdb
}
coordpdb segtype_generic_XJ28.pdb XJ28
############################## Segment XJ28 ends ###############################
# Segment XJ28 ends
# Segment XJ29 begins
##################### Segment XJ29 begins as image of XJ29 #####################
set XJ29 [atomselect $m2 "serial 38428 to 38454"]
$XJ29 set segname XJ29
$XJ29 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ29 writepdb segtype_generic_XJ29.pdb
segment XJ29 {
    first none
    last none
    pdb segtype_generic_XJ29.pdb
}
coordpdb segtype_generic_XJ29.pdb XJ29
############################## Segment XJ29 ends ###############################
# Segment XJ29 ends
# Segment XJ3 begins
###################### Segment XJ3 begins as image of XJ3 ######################
set XJ3 [atomselect $m2 "serial 37785 to 37811"]
$XJ3 set segname XJ3
$XJ3 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ3 writepdb segtype_generic_XJ3.pdb
segment XJ3 {
    first none
    last none
    pdb segtype_generic_XJ3.pdb
}
coordpdb segtype_generic_XJ3.pdb XJ3
############################### Segment XJ3 ends ###############################
# Segment XJ3 ends
# Segment XJ30 begins
##################### Segment XJ30 begins as image of XJ30 #####################
set XJ30 [atomselect $m2 "serial 38455 to 38476"]
$XJ30 set segname XJ30
$XJ30 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ30 writepdb segtype_generic_XJ30.pdb
segment XJ30 {
    first none
    last none
    pdb segtype_generic_XJ30.pdb
}
coordpdb segtype_generic_XJ30.pdb XJ30
############################## Segment XJ30 ends ###############################
# Segment XJ30 ends
# Segment XJ31 begins
##################### Segment XJ31 begins as image of XJ31 #####################
set XJ31 [atomselect $m2 "serial 38477 to 38503"]
$XJ31 set segname XJ31
$XJ31 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ31 writepdb segtype_generic_XJ31.pdb
segment XJ31 {
    first none
    last none
    pdb segtype_generic_XJ31.pdb
}
coordpdb segtype_generic_XJ31.pdb XJ31
############################## Segment XJ31 ends ###############################
# Segment XJ31 ends
# Segment XJ32 begins
##################### Segment XJ32 begins as image of XJ32 #####################
set XJ32 [atomselect $m2 "serial 38504 to 38531"]
$XJ32 set segname XJ32
$XJ32 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ32 writepdb segtype_generic_XJ32.pdb
segment XJ32 {
    first none
    last none
    pdb segtype_generic_XJ32.pdb
}
coordpdb segtype_generic_XJ32.pdb XJ32
############################## Segment XJ32 ends ###############################
# Segment XJ32 ends
# Segment XJ33 begins
##################### Segment XJ33 begins as image of XJ33 #####################
set XJ33 [atomselect $m2 "serial 38532 to 38558"]
$XJ33 set segname XJ33
$XJ33 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ33 writepdb segtype_generic_XJ33.pdb
segment XJ33 {
    first none
    last none
    pdb segtype_generic_XJ33.pdb
}
coordpdb segtype_generic_XJ33.pdb XJ33
############################## Segment XJ33 ends ###############################
# Segment XJ33 ends
# Segment XJ34 begins
##################### Segment XJ34 begins as image of XJ34 #####################
set XJ34 [atomselect $m2 "serial 38559 to 38585"]
$XJ34 set segname XJ34
$XJ34 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ34 writepdb segtype_generic_XJ34.pdb
segment XJ34 {
    first none
    last none
    pdb segtype_generic_XJ34.pdb
}
coordpdb segtype_generic_XJ34.pdb XJ34
############################## Segment XJ34 ends ###############################
# Segment XJ34 ends
# Segment XJ35 begins
##################### Segment XJ35 begins as image of XJ35 #####################
set XJ35 [atomselect $m2 "serial 38586 to 38607"]
$XJ35 set segname XJ35
$XJ35 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ35 writepdb segtype_generic_XJ35.pdb
segment XJ35 {
    first none
    last none
    pdb segtype_generic_XJ35.pdb
}
coordpdb segtype_generic_XJ35.pdb XJ35
############################## Segment XJ35 ends ###############################
# Segment XJ35 ends
# Segment XJ36 begins
##################### Segment XJ36 begins as image of XJ36 #####################
set XJ36 [atomselect $m2 "serial 55853 to 55880"]
$XJ36 set segname XJ36
$XJ36 set resid [list 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601]
$XJ36 writepdb segtype_generic_XJ36.pdb
segment XJ36 {
    first none
    last none
    pdb segtype_generic_XJ36.pdb
}
coordpdb segtype_generic_XJ36.pdb XJ36
############################## Segment XJ36 ends ###############################
# Segment XJ36 ends
# Segment XJ37 begins
##################### Segment XJ37 begins as image of XJ37 #####################
set XJ37 [atomselect $m2 "serial 55881 to 55908"]
$XJ37 set segname XJ37
$XJ37 set resid [list 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602 602]
$XJ37 writepdb segtype_generic_XJ37.pdb
segment XJ37 {
    first none
    last none
    pdb segtype_generic_XJ37.pdb
}
coordpdb segtype_generic_XJ37.pdb XJ37
############################## Segment XJ37 ends ###############################
# Segment XJ37 ends
# Segment XJ38 begins
##################### Segment XJ38 begins as image of XJ38 #####################
set XJ38 [atomselect $m2 "serial 55909 to 55936"]
$XJ38 set segname XJ38
$XJ38 set resid [list 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603 603]
$XJ38 writepdb segtype_generic_XJ38.pdb
segment XJ38 {
    first none
    last none
    pdb segtype_generic_XJ38.pdb
}
coordpdb segtype_generic_XJ38.pdb XJ38
############################## Segment XJ38 ends ###############################
# Segment XJ38 ends
# Segment XJ39 begins
##################### Segment XJ39 begins as image of XJ39 #####################
set XJ39 [atomselect $m2 "serial 55937 to 55964"]
$XJ39 set segname XJ39
$XJ39 set resid [list 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604 604]
$XJ39 writepdb segtype_generic_XJ39.pdb
segment XJ39 {
    first none
    last none
    pdb segtype_generic_XJ39.pdb
}
coordpdb segtype_generic_XJ39.pdb XJ39
############################## Segment XJ39 ends ###############################
# Segment XJ39 ends
# Segment XJ4 begins
###################### Segment XJ4 begins as image of XJ4 ######################
set XJ4 [atomselect $m2 "serial 37812 to 37838"]
$XJ4 set segname XJ4
$XJ4 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ4 writepdb segtype_generic_XJ4.pdb
segment XJ4 {
    first none
    last none
    pdb segtype_generic_XJ4.pdb
}
coordpdb segtype_generic_XJ4.pdb XJ4
############################### Segment XJ4 ends ###############################
# Segment XJ4 ends
# Segment XJ5 begins
###################### Segment XJ5 begins as image of XJ5 ######################
set XJ5 [atomselect $m2 "serial 37839 to 37860"]
$XJ5 set segname XJ5
$XJ5 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ5 writepdb segtype_generic_XJ5.pdb
segment XJ5 {
    first none
    last none
    pdb segtype_generic_XJ5.pdb
}
coordpdb segtype_generic_XJ5.pdb XJ5
############################### Segment XJ5 ends ###############################
# Segment XJ5 ends
# Segment XJ6 begins
###################### Segment XJ6 begins as image of XJ6 ######################
set XJ6 [atomselect $m2 "serial 37861 to 37887"]
$XJ6 set segname XJ6
$XJ6 set resid [list 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
$XJ6 writepdb segtype_generic_XJ6.pdb
segment XJ6 {
    first none
    last none
    pdb segtype_generic_XJ6.pdb
}
coordpdb segtype_generic_XJ6.pdb XJ6
############################### Segment XJ6 ends ###############################
# Segment XJ6 ends
# Segment XJ7 begins
###################### Segment XJ7 begins as image of XJ7 ######################
set XJ7 [atomselect $m2 "serial 37888 to 37914"]
$XJ7 set segname XJ7
$XJ7 set resid [list 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]
$XJ7 writepdb segtype_generic_XJ7.pdb
segment XJ7 {
    first none
    last none
    pdb segtype_generic_XJ7.pdb
}
coordpdb segtype_generic_XJ7.pdb XJ7
############################### Segment XJ7 ends ###############################
# Segment XJ7 ends
# Segment XJ8 begins
###################### Segment XJ8 begins as image of XJ8 ######################
set XJ8 [atomselect $m2 "serial 37915 to 37934"]
$XJ8 set segname XJ8
$XJ8 set resid [list 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]
$XJ8 writepdb segtype_generic_XJ8.pdb
segment XJ8 {
    first none
    last none
    pdb segtype_generic_XJ8.pdb
}
coordpdb segtype_generic_XJ8.pdb XJ8
############################### Segment XJ8 ends ###############################
# Segment XJ8 ends
# Segment XJ9 begins
###################### Segment XJ9 begins as image of XJ9 ######################
set XJ9 [atomselect $m2 "serial 37935 to 37955"]
$XJ9 set segname XJ9
$XJ9 set resid [list 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4]
$XJ9 writepdb segtype_generic_XJ9.pdb
segment XJ9 {
    first none
    last none
    pdb segtype_generic_XJ9.pdb
}
coordpdb segtype_generic_XJ9.pdb XJ9
############################### Segment XJ9 ends ###############################
# Segment XJ9 ends
############################# DISU patches follow ##############################
patch DISU A:598 A:604
patch DISU G:54 G:74
patch DISU G:119 G:205
patch DISU G:126 G:196
patch DISU G:218 G:247
patch DISU G:228 G:239
patch DISU G:296 G:331
patch DISU G:378 G:445
patch DISU G:385 G:418
patch DISU C:16 C:84
patch DISU C:130 C:159
patch DISU B:598 B:604
patch DISU I:54 I:74
patch DISU I:119 I:205
patch DISU I:126 I:196
patch DISU I:218 I:247
patch DISU I:228 I:239
patch DISU I:296 I:331
patch DISU I:378 I:445
patch DISU I:385 I:418
patch DISU E:16 E:84
patch DISU E:130 E:159
patch DISU D:598 D:604
patch DISU J:54 J:74
patch DISU J:119 J:205
patch DISU J:126 J:196
patch DISU J:218 J:247
patch DISU J:228 J:239
patch DISU J:296 J:331
patch DISU J:378 J:445
patch DISU J:385 J:418
patch DISU F:16 F:84
patch DISU F:130 F:159
patch DISU H:22 H:92
patch DISU H:140 H:196
patch DISU K:22 K:92
patch DISU K:140 K:196
patch DISU M:22 M:92
patch DISU M:140 M:196
patch DISU L:23 L:88
patch DISU L:134 L:194
patch DISU N:23 N:88
patch DISU N:134 N:194
patch DISU O:23 O:88
patch DISU O:134 O:194
############################# LINK patches follow ##############################
############################### Transform 1 ends ###############################
guesscoord
regenerate angles dihedrals
writepsf cmap 00-01-00_psfgen-build.psf
writepdb 00-01-00_psfgen-build.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using Pestifer! #########################
