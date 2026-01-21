#!/bin/bash

#SBATCH --job-name="60kpc_box_100Myr_0_1_mpcc_wind_rad_bondi_accr"
#SBATCH -p normal
#SBATCH -t 03-00:00:00  # dd-hh:mm:ss
#SBATCH -n 128
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err.log
#SBATCH --export=ALL


# Unload all currently loaded modules
module purge

# Load the desired openmpi module
module load openmpi/4.1.1

# Print job info for debugging
echo "=== Job started on $(date) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Number of tasks: $SLURM_NTASKS"
echo "Nodes allocated: $SLURM_JOB_NODELIST"
echo "Working directory: $PWD"
echo "Launching Athena with $SLURM_NTASKS MPI processes..."
echo ""

# Run your MPI application
# mpirun -np $SLURM_NTASKS ../athena -i /scratch/meemik/athenak/inputs/hydro/wind_rad_bondi/dens_0_1_mp_30kpc_box_wind_rad_bondi_100Myr.athinput -d box30kpcdens0_1mpccWindRadBondi100MyrOut
mpirun -np $SLURM_NTASKS ../athena -i /scratch/meemik/athenak/inputs/hydro/wind_rad_bondi/dens_0_1_mp_60kpc_box_wind_rad_bondi_100Myr.athinput -d dens0_1mpcc60kpcBoxWindRadBondi100MyrOut

echo ""
echo "=== Job finished on $(date) ==="

# mpirun -np 148 ../athena -r /scratch/meemik/athenak/build_wind_rad_bondi/src/dens0_1mpcc/dens0_1mpcc60kpcBoxWindRadBondi100MyrOut/rst/wind_rad_bondi_accr.00007.rst -d dens0_1mpcc60kpcBoxWindRadBondi100MyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

