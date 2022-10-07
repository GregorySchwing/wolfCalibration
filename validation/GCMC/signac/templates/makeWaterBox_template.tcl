set r R_ARG
# Name it this so I can reuse the current signac script
set output_pdb_psf_file_name OUTPUT

package require solvate
solvate -o $output_pdb_psf_file_name -minmax [list [vecscale $r {-1 -1 -1}] [vecscale $r {1 1 1}]]
