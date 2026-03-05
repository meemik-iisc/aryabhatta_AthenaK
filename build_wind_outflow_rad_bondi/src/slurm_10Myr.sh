#!/bin/bash

#SBATCH --job-name="10Myr_wind_1e3kms_outflow_rad_bondi"
#SBATCH -p normal
#SBATCH -t 01-00:00:00  # dd-hh:mm:ss
#SBATCH -n 64
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err.log
#SBATCH --export=ALL


# Unload all currently loaded modules
module purge

# Load the desired openmpi module
module load openmpi/4.1.1

# Run your MPI application
mpirun -np 64 ./athena -i /scratch/meemik/athenak/inputs/hydro/wind_outflow_rad_bondi/wind_outflow_rad_bondi_high_dens_vw_1e4kms_10kpc_box_10Myr.athinput -d densWindOutflowBox1e4kms10kpcRadBondi10MyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

