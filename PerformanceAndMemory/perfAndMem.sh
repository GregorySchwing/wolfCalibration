#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH --qos=primary
#
# Number of cores
#SBATCH -c 8 -w, 
#SBATCH --nodelist=wsu205
#SBATCH --mem=1500G
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=168:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err


# Clear the environment from any previously loaded modules
#module purge > /dev/null 2>&1
#cd /home6/greg/fz/martini/2M/charmm-gui-2758929829/gromacs
#cd /wsu/home/go/go24/go2432/GOMC
#./metamake.sh
#cp ./bin/GOMC_CPU_NVT /wsu/home/go/go24/go2432/big/15m/15mil
cd /wsu/home/go/go24/go2432/big/15m/15mil
valgrind --tool=massif ./GOMC_CPU_NVT in.conf 
