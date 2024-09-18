######################## pestifer: pestifer-script.tcl #########################
####################### Created Wed Sep 18 12:09:56 2024 #######################
package require psfgen
psfcontext mixedcase
topology top_all36_prot.rtf
topology top_all35_ethers.rtf
topology top_all36_cgenff.rtf
topology top_all36_lipid.rtf
topology top_all36_carb.rtf
topology top_all36_na.rtf
topology toppar_water_ions.str
topology toppar_all36_carb_glycopeptide.str
topology toppar_all36_prot_modify_res.str
topology toppar_all36_moreions.str
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
pdbalias residue NAG BGLCNA
pdbalias residue BGLC BGLCNA
pdbalias residue FUC AFUC
pdbalias residue GAL BGAL
pdbalias residue ANE5 ANE5AC
pdbalias residue SIA ANE5AC
pdbalias residue EIC LIN
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias residue CL CLA
################################### TESTING ####################################
guesscoord
regenerate angles dihedrals
writepsf cmap test-state.psf
writepdb test-state.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using pestifer! #########################
