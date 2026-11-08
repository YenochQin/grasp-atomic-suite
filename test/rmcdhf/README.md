# Retained RMCDHF regression utilities

This directory contains shared checks for the original RMCDHF programs and the
independent strict/diagnostic solvers. The four `rmcdhf_orbopt*` programs and their
dedicated runners, SLURM templates, fixtures and policy tests were removed on
2026-11-08. Historical descriptions are in
[the archive](../../docs/archive/orbopt/README.md); they are not current run instructions.

- `test_start_trace.py`, `test_p0_trace.py` and `p0_trace_mpi.py` check the retained
  diagnostic implementations.
- `strict_mem_restart.py` uses the original `rmcdhf_mem_mpi` as the one-sweep
  comparison baseline. Its strict-convergence assertions remain unchanged.
- `compare_rwfn.py`, `compare_sum.py`, `crosscheck_rwfn.py` and
  `compare_fine_structure.py` remain generic readers/comparison utilities.

CMake registers the applicable tests when their programs and preparation tools
are available. Use gfortran 15 with `FC` and `OMPI_FC` set to
`/usr/bin/gfortran-15`. Run registered tests with:

```sh
ctest --test-dir build-gf15-interfaces \
  -R 'p0_trace|start_trace|strict_mem_restart' --output-on-failure
```

These small software fixtures do not establish physical or active-space convergence.
