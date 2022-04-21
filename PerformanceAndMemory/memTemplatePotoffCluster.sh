#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH -A greg
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each.
#SBATCH -n 4
# Number of cores
#SBATCH -c 1 -w, --nodelist=potoff8
#SBATCH --mem=XXXG
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=168:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err


# Clear the environment from any previously loaded modules
now="$(date)"
printf "Launch date and time %s\n" "$now"
/home/scratch/GOMCBIN/GOMC_CPU_NVT +p4 NVT_water.conf > log.txt
now="$(date)"
printf "End date and time %s\n" "$now"
