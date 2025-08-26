# pestifer.scripters: 00-01-00_psfgen-build.tcl
####################### Created Mon Aug 25 22:56:58 2025 #######################
package require PestiferCRot
namespace import PestiferCRot::*
package require psfgen
psfcontext mixedcase
topology top_all36_prot.rtf
topology top_all36_cgenff.rtf
topology top_all36_lipid.rtf
topology top_all36_carb.rtf
topology top_all36_na.rtf
topology top_all36_prot.rtf
topology top_all35_ethers.rtf
topology top_all36_cgenff.rtf
topology top_all36_lipid.rtf
topology top_all36_carb.rtf
topology top_all36_na.rtf
topology toppar_water_ions.str
topology toppar_all36_carb_glycopeptide.str
topology toppar_all36_prot_modify_res.str
topology toppar_water_ions.str
topology toppar_all36_carb_glycopeptide.str
topology toppar_all36_prot_modify_res.str
topology toppar_all36_moreions.str
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
mol new 6pti.pdb waitfor all
set m1 [molinfo top get id]
set nf [molinfo $m1 get numframes]
if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $m1 }
mol top $m1
############################## Transform 1 begins ##############################
############### The following mappings of A.U. asym ids is used: ###############
# A.U. chain A: Image chain A
# A.U. chain B: Image chain B
# A.U. chain C: Image chain C
############################### Segments follow ################################
# Segment A begins
############################### Segment A begins ###############################
set A00 [atomselect $m1 "serial 1 to 458"]
$A00 set segname A
$A00 writepdb segtype_polymer_A_1_to_57.pdb
segment A {
    pdb segtype_polymer_A_1_to_57.pdb
}
################################ End segment A #################################
###################### Coordinate-specification commands #######################
# Subsegment [0] is a resolved run
coordpdb segtype_polymer_A_1_to_57.pdb A
####################### Intra-segmental terminal patches #######################
############## Restoring A.U. state for all resolved subsegments ###############
################################ Segment A ends ################################
# Segment A ends
# Segment B begins
######################## Segment B begins as image of B ########################
set B [atomselect $m1 "serial 459 to 463"]
$B set segname B
$B set resid [list 100 100 100 100 100]
$B writepdb segtype_generic_B.pdb
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
set C [atomselect $m1 "serial 464 to 536"]
$C set segname C
$C set resid [list 103 105 108 110 111 112 113 115 121 122 126 129 130 134 138 143 144 149 151 152 153 156 157 158 159 160 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247]
$C writepdb segtype_generic_C.pdb
segment C {
    first none
    last none
    pdb segtype_generic_C.pdb
}
coordpdb segtype_generic_C.pdb C
################################ Segment C ends ################################
# Segment C ends
############################# DISU patches follow ##############################
patch DISU A:5 A:55
patch DISU A:14 A:38
patch DISU A:30 A:51
############################# LINK patches follow ##############################
############################### Transform 1 ends ###############################
guesscoord
regenerate angles dihedrals
writepsf cmap 00-01-00_psfgen-build.psf
writepdb 00-01-00_psfgen-build.pdb
exit
###################### END PESTIFER.SCRIPTERS VMD SCRIPT #######################
################### Thank you for using pestifer.scripters! ####################
