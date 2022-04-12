package require solvate 
set xmin 0	
set ymin 0	
set zmin 0	
set xmax 500	
set ymax 500	
set zmax 500	
solvate -o 500_a -minmax [list [list $xmin $ymin $zmin] [list $xmax $ymax $zmax]] 
