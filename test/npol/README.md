# Nuclear-polarization code tests

The application consumers come from the fresh GRASP copies in `src/appl/rmcdhf*`,
not the existing orbopt sources in `src/appl/rmcdhf90*`.

All data are synthetic software fixtures. These tests never start an atomic SCF
solver, use preparation tools, load atomic datasets, or submit jobs.

Select only these tests with `ctest --test-dir <build> -L '^npol_code$' --output-on-failure`.
There are 9 tests in a serial configuration and 14 with MPI and Python available.
CMake prefers the existing Tools Python environment; no packages are installed.

| Test | Evidence |
| --- | --- |
| `npol_code_input` | Absent/disabled/single/dual-channel input; missing/unknown fields, versions, units, types, invalid charge/mass, missing strengths versus zero, NaN/Inf, cutoffs, string limits, I/O failure, malformed/trailing groups |
| `npol_code_kernel` | Independent closed values at origin and cutoff, arbitrary radius, switches, zero strength, extreme powers, subnormals, rejected overflow |
| `npol_code_compose` | Sign, repeated composition, fresh Coulomb capture, reinitialization, reset, failed initialization rollback |
| `npol_code_bounds` | N=1/NNNP, invalid lengths, short arrays, tail sentinels, bad radii, nonfinite and overflowing composition |
| `npol_code_consumers_<variant>` | Actual production YPOT/DACON/RINTI from each of the four variants; deterministic test-only derivative and weighted-quadrature stubs; diagonal/off-diagonal differences and kinetic invariance |
| `npol_code_upstream_copy` | SHA-256 of all 492 imported GRASP files after removing the reviewed NP-only changes; independent of the external source tree |
| `npol_code_entrypoints` | Four call/interface signatures and hook order; no implicit loading in NUCPOT or unsupported solvers; dynamic production initialization driver |
| `npol_code_mpi_broadcast_1/2` | Real MPI, distinct rank working directories, root-only input, every configuration field, full arrays, disabled/missing/zero behavior, root-only summary |
| `npol_code_mpi_error_2` | Root parse/build/output failure, worker invalid/inconsistent length, charge mismatch, composition failure; all ranks reach an error barrier, then nonzero abort |
| `npol_code_output` | Metadata read-back, replacement of longer output, disabled no-output, open/write failure, no published state on failure |

The two MPI consumer tests run at two ranks and sum actual production YPOT
contributions to verify that only one full external potential remains. DACON and
RINTI must still see the unscaled potential on every rank.

Arithmetic tolerances depend on machine epsilon and the bounded fixture scales;
they are not estimates of physical uncertainty. The off-diagonal-only fixture
uses a signed artificial perturbation, distinct from the attractive NP model.

The Python wrapper runs each process in a temporary directory with a 20-second
timeout. MPI tests have additional CTest timeouts and PROCESSORS metadata. Error
tests require both nonzero exit and an explicit marker emitted after all ranks
reach a barrier; a launcher error, crash, or timeout is a failure. On POSIX,
unreadable-file coverage requires a non-root test user. `/dev/full`, where
available, supplies an actual output-write failure through a temporary symlink;
the device itself is never replaced or deleted.

Build logs and CTest results live under the chosen build directories. Test-only
stubs are linked solely to the consumer executables, never to an application.
See the [implementation guide](../../docs/implemented/rmcdhf_nuclear_polarization.md)
for the exact supported scope and recorded validation.
