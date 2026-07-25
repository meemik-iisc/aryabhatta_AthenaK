#!/bin/bash

#SBATCH --job-name="tCGM_tMdot_tvw"
#SBATCH -p normal
#SBATCH -t 01-00:00:00  # dd-hh:mm:ss
#SBATCH -n 16
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err.log
#SBATCH --export=ALL


# Unload all currently loaded modules
module purge

# Load the desired openmpi module
module load openmpi/4.1.1

# Run your MPI application
mpirun -np 16 ./athena -i /scratch/meemik/athenak/inputs/hydro/wind_outflow_tCGM_tvw/tCGM1_tvw2_10kpcBox_10Myr.athinput -d tCGM1_tvw2_10MyrOut
