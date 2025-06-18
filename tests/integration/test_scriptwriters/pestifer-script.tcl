# pestifer: pestifer-script.tcl
####################### Created Tue Jun 17 16:52:30 2025 #######################
package require psfgen
psfcontext mixedcase
topology top_all36_lipid.rtf
topology top_all36_prot.rtf
topology toppar_all36_prot_modify_res.str
topology toppar_all36_carb_glycopeptide.str
topology toppar_all36_moreions.str
topology top_all36_cgenff.rtf
topology toppar_water_ions.str
topology top_all36_carb.rtf
topology top_all36_na.rtf
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
pdbalias residue NDG AGLCNA
pdbalias residue BGLC BGLCNA
pdbalias residue FUC AFUC
pdbalias residue FUL BFUC
pdbalias residue GAL BGAL
pdbalias residue ANE5 ANE5AC
pdbalias residue SIA ANE5AC
pdbalias residue EIC LIN
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias residue CL CLA
pdbalias residue C6DH C6DHPC
pdbalias residue C7DH C7DHPC
################################### TESTING ####################################
guesscoord
regenerate angles dihedrals
writepsf cmap test-state.psf
writepdb test-state.pdb
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using pestifer! #########################
