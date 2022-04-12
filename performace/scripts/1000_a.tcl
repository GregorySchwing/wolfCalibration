package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 1000	
set ymax 1000	
set zmax 1000	
solvate -o 1000_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
