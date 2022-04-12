package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 200	
set ymax 200	
set zmax 200	
solvate -o 200_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
