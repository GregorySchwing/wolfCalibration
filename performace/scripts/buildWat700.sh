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

mkdir -p 700_a
cd 700_a
vmd < 700_a.tcl 
cd ..
# Production MD
tar czfv 700.tar.gz ./700_a
# End of submit file
