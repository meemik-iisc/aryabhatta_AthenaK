#!/bin/bash

#SBATCH --job-name="1Gyr_adia_bondi_accr"
#SBATCH --mail-type=END
#SBATCH --mail-user=meemikroy@iisc.ac.in
#SBATCH -p normal
#SBATCH -t 01-00:00:00  # dd-hh:mm:ss
#SBATCH -n 64
#SBATCH --output=%x-%j.log
#SBATCH --export=ALL


# Unload all currently loaded modules
module purge

# Load the desired openmpi module
module load openmpi/4.1.1

# Run your MPI application
mpirun -np 64 ./athena -i /scratch/meemik/athenak/inputs/hydro/bondi/bondi_adia_1Gyr.athinput -d adia1GyrBondiAccrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

