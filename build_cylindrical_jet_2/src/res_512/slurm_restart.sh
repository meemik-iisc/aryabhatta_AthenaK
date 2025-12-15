#!/bin/bash

#SBATCH --job-name="res_128_cylindrical_jet"
#SBATCH --mail-type=END
#SBATCH --mail-user=meemikroy@iisc.ac.in
#SBATCH -p normal
#SBATCH -t 01-00:00:00  # dd-hh:mm:ss
#SBATCH -n 128
#SBATCH --output=%x-%j.log

# Unload all currently loaded modules
module purge

# Load the desired openmpi module
module load openmpi/4.1.1


# Run your MPI application
mpirun -np 128 ./athena -r /scratch/meemik/athenak/build_cylindrical_jet_2/src/cylindricalRestartOut/rst/jet_amb.00004.rst -d cylindricalRestartOut