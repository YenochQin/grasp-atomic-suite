#!/bin/zsh
#SBATCH -J summarize-b8
#SBATCH -N 1
#SBATCH -p batch
#SBATCH --output=%j_%x.log
#SBATCH --error=%j_%x.log
set -euo pipefail
root=${GRASP_WORKSPACE:?Set GRASP_WORKSPACE to the parent workspace directory}
src=$root/data/rmcdhf_test_data/results/b8-state-sets-521
cmp=$root/grasp-atomic-suite/test/rmcdhf_orbopt/compare_sum.py
for d in "$src"/B8_full "$src"/B8_target_only; do
  python3 "$cmp" "$d/rmcdhf.sum" "$d/archived.sum" --allow-energy-differences --allow-radial-grid-difference --allow-level-differences > "$d/archive_comparison.csv"
done
printf 'status,complete\n' > "$src/summary_status.csv"
