#!/bin/bash

#SBATCH --job-name="time_cgm_dens1e3_wind_1e2kms_outflow_rad_bondi"
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
mpirun -np 16 ./athena -i /scratch/meemik/athenak/inputs/hydro/wind_outflow_time_cgm/time_cgm_dens_1e3_wind_outflow_rad_bondi_vw_1e2kms_10kpc_box_10Myr.athinput -d timeCGMdens1e3windOutflowBox1e2kms10kpcRadBondi10MyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

