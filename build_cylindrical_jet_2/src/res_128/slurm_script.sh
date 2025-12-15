#!/bin/bash

#SBATCH --job-name="res_128_cylindrical_jet"
#SBATCH --mail-type=END
#SBATCH --mail-user=meemikroy@iisc.ac.in
#SBATCH -p normal
#SBATCH -t 01-00:00:00  # dd-hh:mm:ss
#SBATCH -n 32
#SBATCH --output=%x-%j.log

# Unload all currently loaded modules
module purge

# Load the desired openmpi module
module load openmpi/4.1.1


# Set MCA parameters for InfiniBand communication
export OMPI_MCA_btl=self,vader,openib
export OMPI_MCA_btl_openib_if_include=ibp23s0


# Run your MPI application
mpirun -np 128 ../athena -i /scratch/meemik/athenak/inputs/hydro/cylindrical_jet_2/res_128.athinput -d cylindricalTestOut
