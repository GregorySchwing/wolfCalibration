package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 25	
set ymax 25	
set zmax 25	
solvate -o 25_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
