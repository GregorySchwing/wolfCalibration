package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 50	
set ymax 50	
set zmax 50	
solvate -o 50_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
