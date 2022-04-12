package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 700	
set ymax 700	
set zmax 700	
solvate -o 700_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
