# Retained RMCDHF regression utilities

This directory contains shared readers and comparison utilities for the original
RMCDHF programs. The four `rmcdhf_orbopt*` programs and their
dedicated runners, SLURM templates, fixtures and policy tests were removed on
2026-11-08. Historical descriptions are in
[the archive](../../docs/archive/orbopt/README.md); they are not current run instructions.

- The rejected strict/diagnostic variants and their dedicated tests were removed
  on 2026-11-11; see [retirement scope](../../docs/implemented/rejected_rmcdhf_variants.md).
- `compare_rwfn.py`, `compare_sum.py`, `crosscheck_rwfn.py` and
  `compare_fine_structure.py` remain generic readers/comparison utilities.

CMake registers the original-program tests under `test/f2018`. The fixed-sweep
comparison driver owns its stdin fixture and has no dependency on retired solvers.
Use gfortran 15 with `FC` and `OMPI_FC` set to
`/usr/bin/gfortran-15`. Run registered tests with:

```fish
ctest --test-dir build-program-names-release \
  -R 'f2018_first_sweep_driver|f2018_ci_damping_rmcdhf' --output-on-failure
```

These small software fixtures do not establish physical or active-space convergence.
