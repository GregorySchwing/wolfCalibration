package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 35	
set ymax 35	
set zmax 35	
solvate -o 35_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
