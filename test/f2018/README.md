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

`f2018_davidson_contracts` and `f2018_davidson_mpi_1/2/4` call production
DINIT, INIEST/INIEST2, SPICMV2 and GDVD. They check stride/guard values, diagonal
initial estimates, sparse multiplication, and the lowest eigenpair of a fully
coupled 6-by-6 tridiagonal matrix against its analytic eigenvalue and residual
(1e-12). These protect the callback interfaces and the MGS zero-column boundary.
They do not cover the known early-convergence failure of the MPI solver on a
matrix whose starting vector spans an isolated 2-by-2 block; the current report
records that failed additional probe separately. The BLAS contract also checks
SLAMCH/DLAMCH return types. NP path cases live in `test/npol`, including the
legacy static-storage configuration that originally returned an empty path.

## Acceptance profiles and mandatory runner (2026-11-09)

Resolve the installed GNU 15 compiler on the current host; the Linux reference
uses `/usr/bin/gfortran-15`, while a local macOS installation may use a different
path. Set both `FC` and `OMPI_FC` to that exact executable. All profiles require
GNU 15; keep builds from other compilers separate.

```sh
export FC="$(command -v gfortran-15)"
export OMPI_FC="$FC"
./configure.sh --debug --profile f2018-check --serial \
  --build-dir build-gf15-check
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-check
```

Use `--profile f2018-release -DCMAKE_BUILD_TYPE=Release` in a fresh directory
for the O2 gate. Default `f2018-audit` retains static storage and F2018 syntax;
`legacy` retains static storage with legacy syntax. Those two diagnostic
profiles require `check.py --audit` and never establish strict acceptance.
The check/release all-target builds now pass on GNU 15.3/macOS. Each CTest run
has 64 passes and 4 Linux-only skips out of 68, so mandatory acceptance still
**fails**. Do not weaken flags or skip handling to obtain a green result.

`check.py` refreshes configuration, verifies `required.json`, compiler identity,
all compile commands and target coverage, builds all targets, and checks every
CTest result. Missing tests, failures, and skips fail the run. `--inventory-only`
does not run a build or tests and reports that limited status. Logs and source,
compiler, executable-input identities are retained under the build directory;
see [the current validation report](../../docs/reviews/f2018_strict_interfaces_npol_2026-11-09.md).
The required list covers current programs/tests, not the still-missing F00–F10
suite. Labels are `f2018_unit`, `f2018_serial`, `f2018_mpi`,
`f2018_properties`, and `f2018_io`; no real-case test is registered yet.

`configure.sh` also accepts `--no-mpi`, `--no-csfg`, `--no-openmp`,
`--prefix PATH`, and quoted `-D...` arguments. Targets write modules separately;
`gj90` no longer copies into source `bin/` during a build. Default installation
is `<build>/install`, or the explicit prefix. Validation does not install.

## Production CI damping branch

`ci_damping.F90` compiles six production MATRIX implementations with controlled
assembly/eigensolver/input seams and the real BLAS library. It checks FIRST,
nonzero two-root damping/Schmidt orthogonalization, a zero-damping control,
nonzero block offsets and finite normalized orthogonal outputs. Two serial
cases plus four MPI variants at 1/2/4 ranks give 14 tests. This is not an
end-to-end SCF/RCI or CSFG validation. The checked fixture also exposed and now
protects the MPI IENDC allocation's required zero lower bound.

## Complete first-sweep comparison

The comparator requires exactly seven binary outputs and the frozen MCP record
shape, decodes complete finite wavefunction/mixing records, and checks program
completion and the expected single-iteration limit. All stdout/stderr and native
scientific text are compared. Only the anchored date/time lines and the seconds
line immediately following `Wall time:` are normalized. Warnings, scientific
values and stop states remain significant. Every failure leaves a structured
`comparison.json` and the available input/output logs.

Both solvers remain mandatory by default. On macOS, the upstream memory solver
requires missing Linux `/proc/meminfo`; use an explicit `--solvers rmcdhf` for a
**partial** comparison, or run both on Linux. No automatic skip is performed.
Reference and candidate binary directories must be different. Their compiler,
math library and flags still need independent verification; executable hashes
are identities, not an automatic proof of matching build provenance.
