#### BEGIN PESTIFER PSFGEN ####
#### BEGIN HEADER ####
source /home/cfa/Git/pestifer/pestifer/Resources/tcl/modules/src/loopmc.tcl
source /home/cfa/Git/pestifer/pestifer/Resources/tcl/vmdrc.tcl
package require psfgen
psfcontxt mixedcase
topology /home/cfa/charmm/toppar/top_all36_prot.rtf
topology /home/cfa/charmm/toppar/top_all35_ethers.rtf
topology /home/cfa/charmm/toppar/top_all36_cgenff.rtf
topology /home/cfa/charmm/toppar/top_all36_lipid.rtf
topology /home/cfa/charmm/toppar/top_all36_na.rtf
topology /home/cfa/charmm/toppar/stream/carb/toppar_all36_carb_glycopeptide.str
topology /home/cfa/Git/pestifer/pestifer/Resources/charmm/top_all36_carb.rtf
topology /home/cfa/Git/pestifer/pestifer/Resources/charmm/toppar_water_ions.str
residue HIS HSD
atom ILE CD1 CD
residue NAG BGNA
atom BGNA C7 C
atom BGNA O7 O
atom BGNA C8 CT
atom BGNA N2 N
residue SIA ANE5
atom ANE5 C10 C
atom ANE5 C11 CT
atom ANE5 N5 N
atom ANE5 O1A O11
atom ANE5 O1B O12
atom ANE5 O10 O
atom VCG C01 C1
atom VCG C01 C1
atom VCG C02 C2
atom VCG C03 C3
atom VCG C04 C4
atom VCG C05 C5
atom VCG C06 C6
atom VCG C07 C7
atom VCG C08 C8
atom VCG C09 C9
residue EIC LIN
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
#### END HEADER ####
mol new 1gc1 waitfor all
set m2 [molinfo top get id]
mol top $m2
#### BEGIN SEGMENTS ####
#### END SEGMENTS ####
exit
#### END PESTIFER PSFGEN ####
