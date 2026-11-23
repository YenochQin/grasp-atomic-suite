# Retained RMCDHF regression utilities

This directory contains shared readers and comparison utilities for the original
RMCDHF programs. The four `rmcdhf_orbopt*` programs and their
dedicated runners, SLURM templates, fixtures and policy tests were removed on
2026-11-08. See the [historical results](../../docs/archive/orbopt/evidence-0024-orbopt-as2-diagnostics.md)
for trials of the retired programs.

- The rejected strict/diagnostic variants and their dedicated tests were removed
  on 2026-11-11; see [retirement scope](../../docs/implemented/adr-0002-retire-strict-and-diagnostic-rmcdhf.md).
- `compare_rwfn.py`, `compare_sum.py`, `crosscheck_rwfn.py` and
  `compare_fine_structure.py` remain generic readers/comparison utilities.

CMake registers the original-program tests under `test/f2018`. The fixed-sweep
comparison driver owns its stdin fixture and has no dependency on retired solvers.
Configure and compile using the [Suite build instructions](../../README.md#build).
Use CTest to run the registered regressions.

```fish
ctest --test-dir build-gf16 \
  -R 'f2018_first_sweep_driver|f2018_ci_damping_rmcdhf' --output-on-failure
```

These small software fixtures do not establish physical or active-space convergence.
