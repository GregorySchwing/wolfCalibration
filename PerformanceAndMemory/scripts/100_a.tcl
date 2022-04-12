package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 100	
set ymax 100	
set zmax 100	
solvate -o 100_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
