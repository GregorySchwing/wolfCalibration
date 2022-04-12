package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 300	
set ymax 300	
set zmax 300	
solvate -o 300_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
