#!/bin/bash

#SBATCH --job-name="time_dependent_mdot_04_cgm_rho01mpcc_vw_1e4kms_outflow_rad_bondi"
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
mpirun -np 16 ./athena -i /scratch/meemik/athenak/inputs/hydro/wind_outflow_time_dependent_cgm/time_dependent_mdot04_cgm_rho0_1mpcc_vw_1e4kms_10kpc_box_10Myr.athinput -d timeCGMMdot04rho01mpccVw1e4kms10kpcRadBondi10MyrOut
