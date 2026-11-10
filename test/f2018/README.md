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

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
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
`f2018_davidson_early` and `f2018_davidson_early_mpi_1/2/4` additionally exercise
an isolated 2-by-2 block with analytic lowest eigenvalue 1. A directly converged
residual must be accepted before preconditioning/MGS, even in the first loops;
the energy/coefficient three-loop guards remain unchanged. Both cases check
finite values, eigenvalues, explicit residuals and normalization at 1e-12.
`f2018_sparse_empty_rank` calls the production RANGULAR_MPI SETSDA with occupied
and empty partitions and checks its MCP records. Only this diagnostic test
compilation poisons local integers (`-finit-integer=-999999`) to catch the former
uninitialized empty-rank index deterministically. The BLAS contract also checks
SLAMCH/DLAMCH return types. NP path cases live in `test/npol`, including the
legacy static-storage configuration that originally returned an empty path.

`f2018_strict_orbital_consistency` compiles the real strict IMPROV/CONSIS/DAMPOR
with a controlled radial proposal and real QUAD/RINT. An analytic unit-normalized
orbital tests scale factors 2 and 1/4, a real P/Q shape change, and a small change
below the unchanged ACCY. The shape case remains unconverged even when damping
makes the saved step smaller than ACCY. This catches both comparing before
normalization and comparing after damping. Strict SCNSTY now compares normalized,
undamped proposals; its formula, threshold and two-pass stopping gate are unchanged.

The first prioritized warning-repair batch additionally registers:

- `f2018_weights_*`: seven production GETALDWT variants, with equal, standard,
  custom and retry inputs. A test-only integer poison exposes the former undefined
  read extent; no production default initialization is added.
- `f2018_rk_count_*`: three production GENINTRK bodies with a counted SLATER seam,
  simulated rank layouts 1–4, empty partitions, repeated calls and integer-boundary
  arithmetic. These are cache-partition checks, not full RCI or large-count MPI tests.
- `f2018_csfg_dacon` and `f2018_csfg_debug_mpi_1/2/4`: the actual DACON update has
  the same analytic values for all diagnostic modes, and its two CSFG-owned flags
  are broadcast correctly without indexing beyond the shared 30-entry array.
- `f2018_interp_contracts`: finite polynomial values, first node, short/zero tables,
  singleton extrapolation with an explicit unavailable-accuracy message, exact
  singleton lookup, empty-table rejection, and 71 regular points. ACCY is unchanged.
- `f2018_dgemvmpi_1/2/4`: actual N/T/C products, guards, an empty rank/dimension,
  and rejection of unsupported strides before data access/communication.

The existing Davidson tests also check preservation of earlier projected columns
and another root's convergence flag. A removed unreachable branch is distinguished
from a demonstrated numerical bug in the [repair report](../../docs/reviews/f2018_warning_repairs_2026-11-09.md).

## Acceptance profiles and mandatory runner (2026-11-09)

Resolve the installed GNU 15 compiler on the current host; the Linux reference
uses `/usr/bin/gfortran-15`, while a local macOS installation may use a different
path. Set both `FC` and `OMPI_FC` to that exact executable. All profiles require
GNU 15; keep builds from other compilers separate.

```fish
set -gx FC (command -v gfortran-15)
set -gx OMPI_FC $FC
./configure.sh --debug --profile f2018-check --serial \
  --build-dir build-gf15-check
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-check
```

Current production defaults to `f2018-release` with O2 and target OpenMP only.
Use `f2018-check`/Debug for standard/interface/runtime diagnostics and
`f2018-audit`/Release for optimized standard/interface/warning diagnostics.
All use automatic local storage; the legacy profile, implicit-SAVE switches and
the former `check.py --audit` compatibility mode are retired. CSFG is isolated
and off by default. See [current build policy](../../docs/implemented/build_profiles.md).
The strict normalization stage passed both all-target profiles with 74/74
tests; its two restart cases converge in 11 cycles with unchanged inputs,
assertions and thresholds. That record, the earlier Linux 71/73 and macOS 64/68
results remain historical evidence. The matrix/warning stage added one Python
infrastructure test; the repair batch adds 18 regression entries, giving 93 tests
in the historical full configuration. That 32-configuration run passed, totaling 1,960 CTest executions with no failures
or skips and no warning-budget increases. Current registered-scope acceptance
is not full P0–P9 certification.

`check.py` refreshes configuration, verifies `required.json`, compiler identity,
all compile commands and target coverage, builds all targets, verifies the
per-file warning budget in Check/Audit, and checks every CTest result. Production
reports warning collection as not performed rather than declaring warnings clear. Missing diagnostic warning records,
new/increased warnings, missing tests, failures, and skips fail the run. `--inventory-only`
does not run a build or tests and reports that limited status. Logs, source,
compiler and object identities are retained under the build directory;
see [the current validation report](../../docs/reviews/f2018_warning_repairs_2026-11-09.md).
The required list covers current programs/tests, not the still-missing F00–F10
suite. Labels are `f2018_unit`, `f2018_serial`, `f2018_mpi`,
`f2018_properties`, and `f2018_io`; no real-case test is registered yet.

`configure.sh` also accepts `--no-mpi`, `--no-csfg`, `--no-openmp`,
`--prefix PATH`, and quoted `-D...` arguments. Targets write modules separately;
`gj90` no longer copies into source `bin/` during a build. Default installation
is `<build>/install`, or the explicit prefix. Validation does not install.

## Complete configuration matrix and warning baseline

`matrix.py` enumerates Check, Audit and production Release over
`GRASP_BUILD_SERIAL_APPS`, `GRASP_ENABLE_MPI`, and `GRASP_JJ2LSJ_OPENMP`, with
the angular switch following the applicable serial/OpenMP combination: **24 ordinary builds**.
`GRASP_BUILD_CSFG=OFF` in every active matrix case. Serial/MPI both disabled still
builds the shared libraries, jj2lsj and their applicable tests. A partial
`--case` selection never reports full-matrix acceptance.

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
../graspkit-tools/.venv/bin/python test/f2018/matrix.py build-gf15-matrix \
  --parallel-configs 2 --jobs 3
```

This configures, builds and tests every case, checks enabled/disabled target
groups, checks actual OpenMP flags and CSFG exclusion, and enforces the warning baseline for Check/Audit. Each
case retains its CTest XML and acceptance summary; `matrix-validate.json`
reports the full result. MPI execution requires local process communication.
The run is a correctness regression, not a performance benchmark.

Diagnostic profiles default to `GRASP_CAPTURE_WARNINGS=ON`; production uses OFF. A diagnostic compiler launcher
records GNU JSON diagnostics separately for every object, including zero-warning
compilations, and preserves any existing compiler launcher. Records contain the
command and source/object hashes. Parallel log interleaving cannot attach a
warning to the wrong file. `check.py` rejects missing/stale records; after enabling
capture on an older build, use a fresh directory or rebuild with `--clean-first`.
The launcher adds only diagnostic-format flags, preserving numerical options.

The portable [JSON baseline](../../docs/reference/f2018/gfortran15-warning-baseline.json)
has a [per-file index](../../docs/reference/f2018/gfortran15-warning-baseline.md).
Diagnostic identity includes target, source, position, warning option and message;
counts are separate for each configuration. A new diagnostic cannot be offset by
removing another, or by removing it from another target. Source/target coverage and
the exact compiler version are checked too. `-w` and `-Wno-*` cannot hide warnings
in strict acceptance.

`correctness-review` flags potential uninitialized values, bounds, narrowing,
truncation, interface mismatch and implicit static-storage risks. These are
**unresolved review items**, not a list of proven bugs or accepted exceptions.
A warning-regression pass says there is no increase relative to the recorded
baseline; it does not satisfy P3's separate correctness-warning clearance gate.
An uninitialized value is read before a defined value has been established by
assignment, input, initialization or a procedure result; a Fortran declaration
alone does not guarantee zero. `-Wmaybe-uninitialized` can also mean that the
compiler cannot prove the assignment path or that a branch is unreachable.
Diagnostic counts are not counts of variables or confirmed bugs: one issue can
appear at several positions or in several targets. Verify the dataflow and valid
inputs, and repair the actual contract rather than initializing everything to zero.
Linker diagnostics and compiler calls made inside disposable Python test fixtures
remain in their logs and are outside this CMake Fortran compilation baseline.

To explicitly establish a replacement candidate, first collect all builds without
claiming test acceptance. The recorder requires all 16 ordinary diagnostic cases and refuses to
overwrite an existing baseline:

```fish
../graspkit-tools/.venv/bin/python test/f2018/matrix.py build-gf15-matrix \
  --stage build --parallel-configs 2 --jobs 3
../graspkit-tools/.venv/bin/python test/f2018/warning_inventory.py record \
  build-gf15-matrix --output build-gf15-matrix/candidate-warnings.json
../graspkit-tools/.venv/bin/python test/f2018/matrix.py build-gf15-matrix \
  --baseline build-gf15-matrix/candidate-warnings.json
```

Review source/diagnostic differences before intentionally replacing the checked-in
baseline, including reductions after fixes; no validation command updates it
automatically. Rebuild the Markdown index with `warning_inventory.py report`.
The former `--audit` compatibility mode is removed; production validation never certifies strict
acceptance. The infrastructure regression `f2018_matrix_warnings_driver` checks
the full Cartesian matrix and negative parser, coverage, hash, version and budget
cases without running atomic calculations.

```fish
../graspkit-tools/.venv/bin/python test/f2018/warning_inventory.py report \
  --output docs/reference/f2018/gfortran15-warning-baseline.md
```

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
