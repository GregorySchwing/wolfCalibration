#!/bin/bash

#SBATCH --job-name Simple

#SBATCH -q requeue

#SBATCH -N 1

#SBATCH -n 1

#SBATCH --mem=5G

#SBATCH --constraint=avx2

#SBATCH --mail-type=ALL

#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err

#SBATCH -t 2-0:0:0
eval "$(conda shell.bash hook)"
source /wsu/home/go/go24/go2432/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38-new
python project.py status
echo $HOSTNAME
