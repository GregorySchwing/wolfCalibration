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
#SBATCH --time=24:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err


valgrind --tool=massif /wsu/home/go/go24/go2432/GOMC/bin/GOMC_CPU_NVT NVT_water.conf > log
