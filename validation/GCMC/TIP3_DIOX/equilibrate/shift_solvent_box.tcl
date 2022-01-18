set SolventNamePDB LIQ_BOX_NPT_EQ
set SolventNamePSF LIQ_BOX

set output_pdb_psf_file_name LIQ_BOX_NPT_EQ_SHIFTED
#***************************************************************************


# load solvent molecule
set solvent [mol new $SolventNamePSF.psf waitfor all]
mol addfile $SolventNamePDB.pdb mol $solvent waitfor all

set all [atomselect top "all"]


set cell [pbc get -now]

puts "cell $cell"
set cellList [split $cell]
set xAx [lindex $cellList 0]
regsub {\{} $xAx {} xAx 
set yAx [lindex $cellList 1]
set zAx [lindex $cellList 2]
puts "xAx $xAx"
puts "yAx $yAx"
puts "zAx $zAx"

set boxDim [list $xAx $yAx $zAx]
set halfX [expr $xAx/2]
set halfY [expr $yAx/2]
set halfZ [expr $zAx/2]
set trueOrigin [list $halfX $halfY $halfZ]

set cen [measure center [atomselect top all]]
puts "CENTER: $cen" 
set x1 [lindex $cen 0]
set y1 [lindex $cen 1]
set z1 [lindex $cen 2]
puts "x1 $x1"
puts "y1 $y1"
puts "z1 $z1"

set max 0

set geoCenter [list $x1 $y1 $z1]

set transformationVector [vecsub $trueOrigin $geoCenter]
$all moveby $transformationVector

puts "transformationVector $transformationVector"

$all writepsf ../common/$output_pdb_psf_file_name.psf
$all writepdb ../common/$output_pdb_psf_file_name.pdb


