#!/bin/zsh
#SBATCH -J rmcdhf-ni-expected-nodes
#SBATCH -N 1
#SBATCH --ntasks-per-node=46
#SBATCH -p batch
#SBATCH --output=%j_%x.log
#SBATCH --error=%j_%x.log

# Regression for b47ac87: require each candidate to match the orbital's
# NNODEP value.  Job 589 proved that this raw/post-damping criterion rejects
# recoverable intermediate orbitals; keep this script only to reproduce that
# failure boundary.
export GRASP_REJECT_NODE_CHANGE=1
unset GRASP_DEFER_ORTHY
exec zsh ${GRASP_WORKSPACE:?Set GRASP_WORKSPACE}/grasp-atomic-suite/test/rmcdhf_orbopt/slurm/run_ni_as1_as2_chain_diagnostic.sbatch
