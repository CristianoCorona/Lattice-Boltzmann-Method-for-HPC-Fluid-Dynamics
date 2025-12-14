#!/bin/bash
# Example PBS job script for LBM simulation
# Name
#PBS -N LBM_Sim_1000Re                  
# Standard CPU queue
#PBS -q cpu                             
# 1 node, 28  cores
#PBS -l select=1:ncpus=28:mem=4gb
# Max time (2 hours)
#PBS -l walltime=02:00:00

source /software/spack/share/spack/setup-env.sh
spack load gcc@15.2.0

#move to the dir where the script has been launched
cd $PBS_O_WORKDIR

# Debug info
echo "Job running on host: $(hostname)"
echo "Working directory: $(pwd)"

export OMP_NUM_THREADS=32

./LBM_cluster