package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 900	
set ymax 900	
set zmax 900	
solvate -o 900_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
