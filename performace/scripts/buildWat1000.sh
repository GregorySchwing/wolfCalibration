#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH --qos=primary
#
# Number of cores
#SBATCH -c 4 -w, 
#SBATCH --mem=200G
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=168:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err

vmd < 1000_a.tcl 
cd ..
# Production MD
tar czfv 1000.tar.gz ./1000_a
# End of submit file
