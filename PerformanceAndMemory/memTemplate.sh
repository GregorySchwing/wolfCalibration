#!/bin/bash
# Project/Account
#SBATCH --qos=primary
#
# Number of cores
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each.
#SBATCH -n 8
#SBATCH --mem=XXXG
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=14-0:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err

module swap gnu7 intel/2019
now="$(date)"
printf "Launch date and time %s\n" "$now"
/wsu/home/go/go24/go2432/wolfCalibration/PerformanceAndMemory/GOMC_CPU_NVT NVT_water.conf > log.txt
now="$(date)"
printf "End date and time %s\n" "$now"
