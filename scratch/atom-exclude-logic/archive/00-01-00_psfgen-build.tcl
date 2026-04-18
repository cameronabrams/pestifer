# pestifer.scripters: 00-01-00_psfgen-build.tcl
####################### Created Sat Apr 18 15:18:05 2026 #######################
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
mol new 5u3o.pdb waitfor all
set m1 [molinfo top get id]
set nf [molinfo $m1 get numframes]
if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $m1 }
mol top $m1
############################## Transform 1 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
# A.U. chain H: Image chain H
# A.U. chain L: Image chain L
# A.U. chain A: Image chain A
############################### Segments follow ################################
# Segment H begins
############################### Segment H begins ###############################
set H00 [atomselect $m1 "serial 1 to 3563"]
$H00 set segname H
$H00 writepdb segtype_polymer_H_1_to_216.pdb
segment H {
    pdb segtype_polymer_H_1_to_216.pdb
}
################################ End segment H #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_H_1_to_216.pdb H
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment H ends ################################
# Segment H ends
# Segment L begins
############################### Segment L begins ###############################
set L00 [atomselect $m1 "serial 3564 to 6877"]
$L00 set segname L
############ Atom with serial 6878 in PDB needs serial 6877 for VMD ############
$L00 writepdb segtype_polymer_L_1_to_214.pdb
segment L {
    pdb segtype_polymer_L_1_to_214.pdb
}
################################ End segment L #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_L_1_to_214.pdb L
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment L ends ################################
# Segment L ends
# Segment A begins
############################### Segment A begins ###############################
set A01 [atomselect $m1 "serial 6878 to 7170"]
$A01 set segname A
############ Atom with serial 7172 in PDB needs serial 7170 for VMD ############
$A01 writepdb segtype_polymer_A_669_to_685.pdb
segment A {
    pdb segtype_polymer_A_669_to_685.pdb
}
################################ End segment A #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_A_669_to_685.pdb A
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment A ends ################################
# Segment A ends
############################# DISU patches follow ##############################
patch DISU H:22 H:92
patch DISU H:216 L:214
patch DISU L:23 L:88
patch DISU L:134 L:194
############################# LINK patches follow ##############################
############################### Transform 1 ends ###############################
guesscoord
regenerate angles dihedrals
writepsf cmap 00-01-00_psfgen-build.psf
writepdb 00-01-00_psfgen-build.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using Pestifer! #########################
