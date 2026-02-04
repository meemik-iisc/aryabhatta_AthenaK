#!/bin/bash

#SBATCH --job-name="100Myr_dens_1_mpcc_10kpc_box_wind_rad_bondi_accr"
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
mpirun -np $SLURM_NTASKS ../athena -i /scratch/meemik/athenak/inputs/hydro/wind_rad_bondi/dens_1_mp_10kpc_box_wind_rad_bondi_100Myr.athinput -d adiaDens1mpcc10kpcBoxWindRadBondi100MyrAccrOut
# mpirun -np 144 ../athena -r /scratch/meemik/athenak/build_wind_rad_bondi/src/dens1mpcc/dens1mpcc60kpcBoxWindRadBondi100MyrAccrOut/rst/wind_rad_bondi_accr.00008.rst -d dens1mpcc60kpcBoxWindRadBondi100MyrAccrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

