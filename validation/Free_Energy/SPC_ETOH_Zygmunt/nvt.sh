#!/bin/bash
# Job name
#SBATCH --job-name SPC_ZYG_NVT_EQ
# Submit to the gpu QoS, the requeue QoS can also be used for gpu's 
#SBATCH -q primary
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each.
#SBATCH -n 8
# Request memory
#SBATCH --mem=8G
# Mail when the job begins, ends, fails, requeues
# Where to send email alerts
# Create an output file that will be output_<jobid>.out
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit
#SBATCH -t 1-12:0:0
/wsu/home/go/go24/go2432/GOMC/bin/GOMC_CPU_NVT +p8 NVT_Eq_water_ethanol_fe.conf > log














