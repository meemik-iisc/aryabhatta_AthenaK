#!/bin/bash

#SBATCH --job-name="pOut_uCGM_vw3_10MyrOut"
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
mpirun -np 16 ./athena -i /scratch/meemik/athenak/inputs/hydro/pOut_uCGM/pOut_vw3_ton1Myr_tduty_05_uCGM_10kyrBox_10MyrOut.athinput -d pOuttOn1MyrtDuty05uCGMvw310kpcBox10MyrOut
# mpirun -np 4 --mca orte_base_help_aggregate 0 --mca orte_debug_daemons 1 ../athena -i /scratch/meemik/athenak/inputs/hydro/conical_jet/res_128.athinput -d conicalJetAmbHalfOut

