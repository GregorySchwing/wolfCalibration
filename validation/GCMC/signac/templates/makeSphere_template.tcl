set r R_ARG
set padding PADDING_ARG
# Name it this so I can reuse the current signac script
set output_pdb_psf_file_name OUTPUT

set r2 [expr $r*$r]
set r_box [expr $r+$padding]
set 2_r_box [expr $r_box+$r_box]
set ratio [expr (sqrt(5)+1)*0.5]
set listOfBoxHalfAx [list $r_box $r_box $r_box]
set listOf2BoxHalfAx [list $2_r_box $2_r_box $2_r_box]

set numShellAtoms [expr round($M_PI * pow($r,2) * 16 / 9 / sqrt(3))]

package require psfgen
topology GAS_TOPOLOGY
segment GAS {
    auto none
    first NONE
    last NONE
    for {set res 1} {$res < $numShellAtoms} {incr res} {
        residue $res NE1
    }
}
writepsf sphere.psf
writepdb sphere.pdb
resetpsf
 
mol load psf sphere.psf pdb sphere.pdb
set all [atomselect top all]

set max [$all num]
for {set i 0} {$i<$max} {incr i} {
    set atom [atomselect top "index $i"]
    set theta [expr 2 * $M_PI * $i / $ratio]
    set phi [expr acos(1 - 2*($i+0.5)/$max)]
    $atom set x [expr cos($theta)*sin($phi)*$r]
    $atom set y [expr sin($theta)*sin($phi)*$r]
    $atom set z [expr cos($phi)*$r]
    $atom delete
}

$all writepdb sphere.pdb
$all delete
mol delete top

package require solvate
solvate sphere.psf sphere.pdb -minmax [list [vecscale $r_box {-1 -1 -1}] [vecscale $r_box {1 1 1}]]

set all [atomselect top "all"]

package require pbctools
pbc set $listOf2BoxHalfAx -all
$all moveby $listOfBoxHalfAx
$all writepsf $output_pdb_psf_file_name.psf
$all writepdb $output_pdb_psf_file_name.pdb

