#!/bin/bash

#SBATCH --job-name="dens1e3mpcc_01pcbubble_10kyr_wind_rad_bondi_accr"
#SBATCH --mail-type=END
#SBATCH --mail-user=meemikroy@iisc.ac.in
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
mpirun -np 16 ../athena -i /scratch/meemik/athenak/inputs/hydro/wind_rad_bondi/bubble_dens1e3mpcc_box10pc_wind_rad_bondi_10kyr.athinput -d bubbledens1e3mpcc01pcWindRadBondi10kyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

