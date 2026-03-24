#!/bin/bash

#SBATCH --job-name="ton_1Myr_mdot4_vw1e3_dens4e3_outflow_cgm"
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
mpirun -np 16 ./athena -i /scratch/meemik/athenak/inputs/hydro/time_outflow_cgm_wind_rad_bondi/ton_1Myr_tduty_05_mdot_4_cgm_dens_4e3_wind_outflow_rad_bondi_vw_1e3kms_20kpc_box_10Myr.athinput -d ton1Myrtduty05CGMMdot4dens4e3windOutflowBox1e3kms20kpcRadBondi10MyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

