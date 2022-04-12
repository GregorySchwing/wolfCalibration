package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 400	
set ymax 400	
set zmax 400	
solvate -o 400_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
