# F2018 compatibility regressions

These are software checks. They do not establish atomic-model convergence or
physical accuracy. Use gfortran 15 and set `FC` and `OMPI_FC` to
`/usr/bin/gfortran-15` before configuration and validation.

The regular CMake/CTest flow registers:

- `f2018_nmtejj_parity` and `f2018_nmtejj_csfg_parity`: 22,701 integer lookup
  inputs per implementation, compared with a frozen stdout SHA-256 from the
  d1260ba GNU 15.3.0 legacy libraries. The fixture records the original source
  hashes. These cover hits, misses, all three table branches and the upper guard.
- `f2018_csfg_reallocation`: preserve values across shrinking/growing arrays
  using LONG extents, and preserve the restart counter's 64-bit width.
- `f2018_storage_bytes`: scalar, empty, multidimensional, strided and packed
  intrinsic payloads. Arrays must be allocated/associated when passed.
- `f2018_system_paths`: real production POSIX/MPI helper interfaces in a
  disposable directory; no MPI launcher is required for the successful helper
  path. Tests include quoted paths, failure states, short buffers, shell wait
  status, exit status and buffered file output. Run on Linux/POSIX.

The separate baseline comparison requires both a preserved legacy build and the
candidate build. It generates a synthetic finite-nucleus Z=10 `1s^2` input,
runs preparation followed by **one SCF iteration**, and compares binary
wavefunction/mixing/MCP outputs for the two retained original serial RMCDHF variants. Match compiler,
math library, options and thread counts; the script does not certify the supplied
binary directories' provenance automatically.

```sh
export FC=/usr/bin/gfortran-15
export OMPI_FC=/usr/bin/gfortran-15
../graspkit-tools/.venv/bin/python test/f2018/compare_first_sweep.py \
  --reference-bindir build-gf15-standard-import/bin \
  --candidate-bindir build-gf15-f2018-batch/bin \
  --output build-gf15-f2018-batch/audit/first-sweep-parity
```

The output directory must not already exist. Inputs, stdout/stderr and per-file
hashes remain there for diagnosis. Success is byte equality; tolerances are not
relaxed if a difference occurs. This check does not cover CI damping, multi-root
SCF convergence, full RCI or transition properties.

The interface stage additionally registers `f2018_linear_algebra_contracts`,
`f2018_interface_rejections`, and `f2018_mpi_contracts_1/2/4`. These exercise
LP64 BLAS/LAPACK values, rejected invalid call patterns, and actual MPI payloads.
The orbopt programs are retired; their old binaries are not current test targets.
