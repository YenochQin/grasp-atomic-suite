#!/bin/zsh
#SBATCH -J rmcdhf-ni-guard-on
#SBATCH -N 1
#SBATCH --ntasks-per-node=46
#SBATCH -p batch
#SBATCH --output=%j_%x.log
#SBATCH --error=%j_%x.log
export GRASP_REJECT_NODE_CHANGE=1
export GRASP_RESULT_TAG=ni-guard-on-46mpi
exec zsh ${GRASP_WORKSPACE:?Set GRASP_WORKSPACE}/grasp-atomic-suite/test/rmcdhf_orbopt/slurm/run_ni_as1_as2_chain_diagnostic.sbatch
