#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH --qos=primary
#
# Number of cores
#SBATCH -c 4 -w, 
#SBATCH --nodelist=wsu205
#SBATCH --mem=200G
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=168:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err

mkdir -p 600_a
cd 600_a
vmd < 600_a.tcl 
cd ..
# Production MD
tar czfv 600.tar.gz ./600_a
# End of submit file
