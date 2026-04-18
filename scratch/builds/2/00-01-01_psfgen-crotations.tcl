# pestifer.scripters: 00-01-01_psfgen-crotations.tcl
####################### Created Fri Apr 17 12:31:13 2026 #######################
package require PestiferCRot
namespace import PestiferCRot::*
mol new 00-01-00_psfgen-build.psf
mol addfile 00-01-00_psfgen-build.pdb waitfor all
set mCM [molinfo top get id]
set r1 [[atomselect $mCM "chain B and resid 35 and name CA"] get residue]
set r2 [[atomselect $mCM "chain B and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 psi C -180.0
set r1 [[atomselect $mCM "chain B and resid 35 and name CA"] get residue]
set r2 [[atomselect $mCM "chain B and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -60.0
set r1 [[atomselect $mCM "chain B and resid 39 and name CA"] get residue]
set r2 [[atomselect $mCM "chain B and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -30.0
set r1 [[atomselect $mCM "chain D and resid 35 and name CA"] get residue]
set r2 [[atomselect $mCM "chain D and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 psi C -180.0
set r1 [[atomselect $mCM "chain D and resid 35 and name CA"] get residue]
set r2 [[atomselect $mCM "chain D and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -60.0
set r1 [[atomselect $mCM "chain D and resid 39 and name CA"] get residue]
set r2 [[atomselect $mCM "chain D and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -30.0
set r1 [[atomselect $mCM "chain F and resid 35 and name CA"] get residue]
set r2 [[atomselect $mCM "chain F and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 psi C -180.0
set r1 [[atomselect $mCM "chain F and resid 35 and name CA"] get residue]
set r2 [[atomselect $mCM "chain F and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -60.0
set r1 [[atomselect $mCM "chain F and resid 39 and name CA"] get residue]
set r2 [[atomselect $mCM "chain F and resid 57 and name CA"] get residue]
brot $mCM $r1 $r2 phi C -30.0
set TMP [atomselect $mCM all]
$TMP writepdb 00-01-01_psfgen-crotations.pdb
$TMP delete
exit
########################### END PESTIFER VMD SCRIPT ############################
######################## Thank you for using Pestifer! #########################
