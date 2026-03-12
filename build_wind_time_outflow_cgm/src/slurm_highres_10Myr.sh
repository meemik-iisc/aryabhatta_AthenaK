#!/bin/bash

#SBATCH --job-name="high_res_time_outflow_cgm"
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
mpirun -np 64 ./athena -i /scratch/meemik/athenak/inputs/hydro/time_outflow_cgm_wind_rad_bondi/high_res_time_5Myr_cgm_dens_4e3_wind_outflow_rad_bondi_vw_5e3kms_10kpc_box_10Myr.athinput -d highRestime5MyrCGMdens4e3windOutflowBox5e3kms10kpcRadBondi10MyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

