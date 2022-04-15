#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH --qos=secondary
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each. 
#SBATCH -n 1
#SBATCH --mem=5G
#
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=1:00:00
#
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err
mkdir -p ../systems

declare -a SMALL_SYS=(25 35 50)


for r in "${SMALL_SYS[@]}"
 do
    mkdir -p ../systems/${r}_a
    cp buildWat${r}.sh ../systems/${r}_a
    cp ${r}_a.tcl ../systems/${r}_a
    cd ../systems/${r}_a
    sbatch buildWat${r}.sh
    cd ../../scripts
 done

for r in {1..10..1}
 do
    mkdir -p ../systems/${r}00_a
    cp buildWat${r}00.sh ../systems/${r}00_a
    cp ${r}00_a.tcl ../systems/${r}00_a
    cd ../systems/${r}00_a
    sbatch buildWat${r}00.sh
    cd ../../scripts
 done

