package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 800	
set ymax 800	
set zmax 800	
solvate -o 800_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
