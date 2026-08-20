# pestifer.scripters: 00-01-000_psfgen-build.tcl
######################## pestifer 3.19.0  seed 27021972 ########################
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
topology toppar_all36_carb_imlab.str
topology toppar_all36_prot_modify_res.str
topology toppar_all36_moreions.str
topology top_all35_ethers.rtf
pdbalias atom ILE CD1 CD
pdbalias atom MET SE SD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5AC C10 C
pdbalias atom ANE5AC C11 CT
pdbalias atom ANE5AC N5 N
pdbalias atom ANE5AC O1A O11
pdbalias atom ANE5AC O1B O12
pdbalias atom ANE5AC O10 O
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
pdbalias atom CLA CL CLA
pdbalias atom ACT CH3 C1
pdbalias atom ACT O O1
pdbalias atom ACT OXT O2
pdbalias atom ACT C C2
pdbalias atom CAL CA CAL
pdbalias atom ILE CD1 CD
pdbalias atom MET SE SD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5AC C10 C
pdbalias atom ANE5AC C11 CT
pdbalias atom ANE5AC N5 N
pdbalias atom ANE5AC O1A O11
pdbalias atom ANE5AC O1B O12
pdbalias atom ANE5AC O10 O
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
pdbalias atom CLA CL CLA
pdbalias atom ACT CH3 C1
pdbalias atom ACT O O1
pdbalias atom ACT OXT O2
pdbalias atom ACT C C2
pdbalias atom CAL CA CAL
pdbalias atom ILE CD1 CD
pdbalias atom MET SE SD
pdbalias atom BGLCNA C7 C
pdbalias atom BGLCNA O7 O
pdbalias atom BGLCNA C8 CT
pdbalias atom BGLCNA N2 N
pdbalias atom ANE5AC C10 C
pdbalias atom ANE5AC C11 CT
pdbalias atom ANE5AC N5 N
pdbalias atom ANE5AC O1A O11
pdbalias atom ANE5AC O1B O12
pdbalias atom ANE5AC O10 O
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
pdbalias atom CLA CL CLA
pdbalias atom ACT CH3 C1
pdbalias atom ACT O O1
pdbalias atom ACT OXT O2
pdbalias atom ACT C C2
pdbalias atom CAL CA CAL
pdbalias residue HIS HSD
pdbalias residue MSE MET
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
pdbalias residue ACT ACET
pdbalias residue CA CAL
pdbalias residue HIS HSD
pdbalias residue MSE MET
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
pdbalias residue ACT ACET
pdbalias residue CA CAL
pdbalias residue HIS HSD
pdbalias residue MSE MET
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
pdbalias residue ACT ACET
pdbalias residue CA CAL
mol new my_6pti.pdb waitfor all autobonds off
set m27 [molinfo top get id]
set nf [molinfo $m27 get numframes]
if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $m27 }
mol top $m27
############################# Transform 50 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
############################### Segments follow ################################
# Segment A begins
############################### Segment A begins ###############################
segment A {
    pdb segtype_polymer_A_1_to_57.pdb
    mutate 24 ALA
}
################################ End segment A #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_A_1_to_57.pdb A
####################### Intra-segmental terminal patches #######################
################################ Segment A ends ################################
# Segment A ends
# Segment B begins
######################## Segment B begins as image of B ########################
segment B {
    first none
    last none
    pdb segtype_generic_B.pdb
}
coordpdb segtype_generic_B.pdb B
################################ Segment B ends ################################
# Segment B ends
# Segment C begins
######################## Segment C begins as image of C ########################
segment C {
    first none
    last none
    pdb segtype_generic_C.pdb
}
coordpdb segtype_generic_C.pdb C
################################ Segment C ends ################################
# Segment C ends
# Segment ION begins
###################### Segment ION begins as image of ION ######################
segment ION {
    first none
    last none
    pdb segtype_generic_ION.pdb
}
coordpdb segtype_generic_ION.pdb ION
############################### Segment ION ends ###############################
# Segment ION ends
# Segment WT1 begins
###################### Segment WT1 begins as image of WT1 ######################
segment WT1 {
    first none
    last none
    pdb segtype_generic_WT1.pdb
}
coordpdb segtype_generic_WT1.pdb WT1
############################### Segment WT1 ends ###############################
# Segment WT1 ends
############################# DISU patches follow ##############################
patch DISU A:5 A:55
patch DISU A:14 A:38
patch DISU A:30 A:51
############################# LINK patches follow ##############################
############################## Transform 50 ends ###############################
guesscoord
regenerate angles dihedrals
writepsf cmap 00-01-000_psfgen-build.psf
writepdb 00-01-000_psfgen-build.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using Pestifer! #########################
