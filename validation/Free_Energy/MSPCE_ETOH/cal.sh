#!/bin/bash
# Job name
#SBATCH --job-name NVT
# Submit to the gpu QoS, the requeue QoS can also be used for gpu's 
#SBATCH -q secondary
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each.
#SBATCH -n 8
# Request memory
#SBATCH --mem=8G
# Mail when the job begins, ends, fails, requeues
#SBATCH --mail-type=ALL
# Where to send email alerts
#SBATCH --mail-user=go2432@wayne.edu
# Create an output file that will be output_<jobid>.out
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit
#SBATCH -t 1-0:0:0
/wsu/home/go/go24/go2432/GOMC2/GOMC/bin/GOMC_CPU_NPT +p8 NPT_Prod_water_ethanol_fe.conf > log
