#!/bin/bash
# Job name
#SBATCH --job-name Simple
# Submit to the secondary QoS
#SBATCH -q secondary
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each. 
#SBATCH -n 1
# Request memory
#SBATCH --mem=1G
# Mail when the job begins, ends, fails, requeues 
#SBATCH --mail-type=ALL
# Where to send email alerts
#SBATCH --mail-user=go2432@wayne.edu
# Create an output file that will be output_<jobid>.out 
#SBATCH -o output_%j.out
# Create an error file that will be error_<jobid>.out
#SBATCH -e errors_%j.err
# Set maximum time limit 
#SBATCH -t 0-2:00:0
for i in {0..22..1}
  do 
     cp npt_prod.sh
     cd Prod_Ewald/state_$i
     sbatch npt_prod.sh
     cd ../..
 done
