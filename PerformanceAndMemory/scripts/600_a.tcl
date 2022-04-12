package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 600	
set ymax 600	
set zmax 600	
solvate -o 600_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
