# pestifer.scripters: 00-01-01_psfgen-crotations.tcl
####################### Created Thu Apr 16 16:18:36 2026 #######################
package require PestiferCRot
namespace import PestiferCRot::*
mol new 00-01-00_psfgen-build.psf
mol addfile 00-01-00_psfgen-build.pdb waitfor all
set mCM [molinfo top get id]
set r1 [[atomselect $mCM "chain A and resid 682 and name CA"] get residue]
set r2 [[atomselect $mCM "chain A and resid 710 and name CA"] get residue]
set rterm [[atomselect $mCM "chain A and resid 710 and name CA"] get residue]
fold_alpha $r1 $r2 $rterm $mCM
set r1 [[atomselect $mCM "chain B and resid 682 and name CA"] get residue]
set r2 [[atomselect $mCM "chain B and resid 710 and name CA"] get residue]
set rterm [[atomselect $mCM "chain B and resid 710 and name CA"] get residue]
fold_alpha $r1 $r2 $rterm $mCM
set r1 [[atomselect $mCM "chain D and resid 682 and name CA"] get residue]
set r2 [[atomselect $mCM "chain D and resid 710 and name CA"] get residue]
set rterm [[atomselect $mCM "chain D and resid 710 and name CA"] get residue]
fold_alpha $r1 $r2 $rterm $mCM
set r1 [[atomselect $mCM "chain A and resid 683 and name CA"] get residue]
set r2 [[atomselect $mCM "chain A and resid 710 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -10.0
set r1 [[atomselect $mCM "chain B and resid 683 and name CA"] get residue]
set r2 [[atomselect $mCM "chain B and resid 710 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -10.0
set r1 [[atomselect $mCM "chain D and resid 683 and name CA"] get residue]
set r2 [[atomselect $mCM "chain D and resid 710 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -10.0
set TMP [atomselect $mCM all]
$TMP writepdb 00-01-01_psfgen-crotations.pdb
$TMP delete
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using Pestifer! #########################
