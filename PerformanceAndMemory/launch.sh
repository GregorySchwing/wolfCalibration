#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH -c 1 -w, --nodelist=potoff32
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each. 
#SBATCH -n 1
#SBATCH --mem=5G
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=4-0:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err
source /home6/greg/anaconda3/etc
conda activate mamba
python3 MemoryTest.py
