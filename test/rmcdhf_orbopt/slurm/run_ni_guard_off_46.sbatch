#!/bin/zsh
#SBATCH -J rmcdhf-ni-guard-off
#SBATCH -N 1
#SBATCH --ntasks-per-node=46
#SBATCH -p batch
#SBATCH --output=%j_%x.log
#SBATCH --error=%j_%x.log
export GRASP_REJECT_NODE_CHANGE=0
export GRASP_RESULT_TAG=ni-guard-off-46mpi
exec zsh ${GRASP_WORKSPACE:?Set GRASP_WORKSPACE}/grasp-atomic-suite/test/rmcdhf_orbopt/slurm/run_ni_as1_as2_chain_diagnostic.sbatch
