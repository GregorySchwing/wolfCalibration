#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH -A greg
#
# Number of cores
#SBATCH -c 6 -w, --nodelist=potoff33
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=168:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err


# Clear the environment from any previously loaded modules
#/home6/greg/GOMC/bin/GOMC_CPU_NPT NPT_Prod_water_ethanol_fe.conf > log 
/wsu/home/go/go24/go2432/GOMC2/GOMC/bin/GOMC_CPU_NPT NPT_Prod_water_ethanol_fe.conf > log
# End of submit file
