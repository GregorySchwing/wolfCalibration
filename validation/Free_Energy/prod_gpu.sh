#!/bin/bash
# Job name
#SBATCH --job-name NVT_GPU
# Submit to the gpu QoS, the requeue QoS can also be used for gpu's 
#SBATCH -q gpu
#SBATCH --gres=gpu:1
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each.
#SBATCH -n 1
# Request memory
#SBATCH --mem=8G
# Create an output file that will be output_<jobid>.out
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit
#SBATCH -t 1-0:0:0
module swap gnu7 intel/2019
/wsu/home/go/go24/go2432/testMPandMEMC/GOMC/bin/GOMC_GPU_NVT NVT_Prod_water_ethanol_fe.conf > log
