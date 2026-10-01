# grasp-atomic-suite

[中文](README_ZH.md) | English

`grasp-atomic-suite` is a Fortran research and development repository for relativistic atomic-structure and atomic-property calculations. It combines the orbital optimization programs migrated from `rmcdhf_test` with `rhfs90`, `ris4`, `gj90`, and their MPI support for hyperfine-structure, isotope-shift, and Landé `g_J` calculations.

This repository is not a minimal upstream GRASP mirror. It keeps the traditional GRASP numerical libraries and application layout, while adding CMake builds, MPI targets, an extracted `gj90` program, regression data, and development notes for `g_J` and MPI optimization work.

## Source Layout

```text
.
├── configure.sh              # out-of-source CMake configuration helper
├── CMakeLists.txt            # top-level CMake build
├── src/
│   ├── appl/
│   │   ├── gj90/             # standalone Landé g_J calculator
│   │   ├── rhfs90/           # relativistic hyperfine-structure program
│   │   ├── rmcdhf90*/        # four rmcdhf_orbopt orbital optimization variants
│   │   └── ris4/             # relativistic isotope-shift program
│   └── lib/
│       ├── libmod/           # global parameters, common state, shared modules
│       ├── lib9290/          # GRASP92-style I/O, constants, grids, utilities
│       ├── libdvd90/         # diagonalization-related routines
│       ├── libmcp90/         # MCP support routines
│       ├── librang90/        # angular algebra, Racah/tensor matrix elements
│       └── mpi90/            # MPI file, path, and parallel helper routines
├── docs/                     # theory notes, implementation traces, optimization reports
├── data/                     # sample inputs and reference outputs for gj90/RHFS checks
├── test/                     # orbital optimization and shared-library tests
├── scripts/                  # grid patcher and build cleanup helpers
├── bin/                      # installed executables
└── lib/                      # installed static libraries and Fortran modules
```

## Programs

| Program | Source directory | CMake targets | Purpose |
| --- | --- | --- | --- |
| `gj90` | `src/appl/gj90` | `gj90`, `gj90_mpi` | Standalone Landé `g_J` calculator extracted from the RHFS path. It reads `isodata`, `name.c`, `name.m/name.cm`, and `name.w`, then writes `name.gj/name.cgj`. |
| `rhfs90` | `src/appl/rhfs90` | `rhfs`, `rhfs_mpi` | Relativistic hyperfine-structure program. It computes hyperfine constants and related matrix elements, writing `name.h/name.ch` and `name.hoffd/name.choffd`. |
| `ris4` | `src/appl/ris4` | `ris4`, `ris4_mpi` | Relativistic isotope-shift program. It computes normal mass shift, specific mass shift, and field-shift electronic factors, writing `name.i/name.ci` and angular intermediate data. |
| Orbital optimization | `src/appl/rmcdhf90*` | `rmcdhf_orbopt`, `rmcdhf_orbopt_mpi`, `rmcdhf_orbopt_mem`, `rmcdhf_orbopt_mem_mpi` | Serial/MPI variants with conventional or in-memory MCP storage. Existing `rwfn.out`, `rmix.out`, `rmcdhf.sum`, and `rmcdhf.log` conventions are retained. |

The default CMake configuration attempts to build MPI versions. If no MPI Fortran toolchain is found, MPI targets are skipped. Serial application targets are disabled by default and must be enabled explicitly at CMake configuration time.

## Build

The recommended build flow is out-of-source CMake:

```sh
./configure.sh
cmake --build build --target install -j4
```

`configure.sh` passes `-DGRASP_ENABLE_MPI=ON` by default. If MPI Fortran is available, installation produces `bin/gj90_mpi`, `bin/rhfs_mpi`, `bin/ris4_mpi`, both `rmcdhf_orbopt*_mpi` programs, and the associated libraries.

To build all four orbital programs and all property programs together:

```sh
cmake -S . -B build-all -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-all --parallel 4
ctest --test-dir build-all --output-on-failure
cmake --install build-all
```

Debug build:

```sh
./configure.sh --debug
cmake --build build-debug --target install -j4
```

To disable MPI and build the serial applications, configure CMake directly:

```sh
mkdir build-serial
cmake -S . -B build-serial -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-serial --target install -j4
```

Remove an existing build directory before reconfiguring, or run:

```sh
./scripts/clean-build-artifacts.sh
```

## Dependencies And Configuration

The build requires a Fortran compiler, CMake, BLAS, and LAPACK. MPI targets additionally require an MPI Fortran environment such as `mpifort`.

CMake finds BLAS/LAPACK and adds compatibility flags for older Fortran sources when using GNU Fortran:

- `-fno-automatic`
- `-fallow-argument-mismatch`

## Inputs And Outputs

The property applications use GRASP-style state-name input. After you enter `name`, the program looks for matching files in the current working directory. Orbital programs use the upstream RMCDHF interactive input and files such as `rcsf.inp`, `rwfn.inp`, and MCP data. Preparation tools such as `rnucleus`, `rwfnestimate`, and `rangular_mpi` still come from an external GRASP installation.

Common input files:

- `isodata`: nuclear and isotope data.
- `name.c`: CSF list and coupling information.
- `name.w`: radial orbital wavefunctions.
- `name.m`: MCDHF/non-CI mixing-coefficient file.
- `name.cm`: CI mixing-coefficient file.

Common output files:

- `gj90`: `name.gj` or `name.cgj`.
- `rhfs90`: `name.h/name.ch` and `name.hoffd/name.choffd`.
- `ris4`: `name.i/name.ci`, and possibly angular intermediate files such as `name.IOB` and `name.ITB`.

`gj90_mpi` also supports command-line arguments:

```sh
mpiexec -n 2 ./bin/gj90_mpi test --nonci
mpiexec -n 2 ./bin/gj90_mpi test --ci
```

Without arguments, it falls back to the interactive input flow.

## Verification Data

The `data/` directory contains a `test` data set:

```text
data/isodata
data/test.c
data/test.m
data/test.w
data/test.gj
```

These files can be used to compare `gj90` Landé-factor output. CTest covers shared-library integration, sparse MPI buffers, orbital transaction/round logic, and the grid patcher:

```sh
ctest --test-dir build-all --output-on-failure
```

Orbital regression runners live in [test/rmcdhf_orbopt](test/rmcdhf_orbopt/README.md). Large calculation inputs/results remain outside the repository. SLURM templates require `GRASP_WORKSPACE` (the parent workspace); external input templates additionally require `GRASP_TEST_DATA_ROOT`. Set cluster module names and resource requests for your host before submitting jobs.

## Radial Grid And Migration

All shared limits and numerical defaults live in
[`suite_parameters_M.f90`](src/lib/libmod/suite_parameters_M.f90).
`NNN1=NNNP+10` and the orbital limit offsets are derived automatically.
All orbital and property programs now default to `N=NNNP` for finite nuclei
(currently 2990); point nuclei keep their separate defaults.
Input orbitals are interpolated onto the current calculation grid.
The shared initializer updates the default ACCY after interactive H changes,
while preserving subsequent explicit tolerance overrides.
Edit the configuration at the top of [scripts/patch_grasp_grid.py](scripts/README.md)
and run it without arguments. Select `SOURCE_LAYOUT="atomic-suite"` for this tree
or `"grasp2018"` for the full upstream tree; the original CLI remains available.
Preview is the default. Keep a separate source copy for each parameter set and
build it manually using the upstream README CMake workflow. Changes require a
full rebuild and numerical convergence checks.

See [docs/common_parameters.md](docs/common_parameters.md) for the parameter table,
initialization order, and validation. The suite layout now edits one central file.
For candidate values, their rationale, test datasets, and convergence criteria,
see [GRASP grid parameter tuning](docs/grasp_grid_parameter_tuning.md) (Chinese).

See [docs/rmcdhf_migration.md](docs/rmcdhf_migration.md) for provenance, shared-library decisions, and validation. New orbital development belongs here; the original `rmcdhf_test` checkout is retained for reference.

## Development Notes

- `docs/lande_g.md`: theory background for Landé `g_J`, MCDHF/RCI evaluation, and LS-coupling checks.
- `docs/RHFS_gJ_report.md`: implementation trace for `g_J` in `rhfs90`, including input loading, radial integrals, angular matrix elements, and ASF projection.
- `docs/gj90_blueprint.md`: standalone `gj90` design, minimal dependency chain, and compute-kernel extraction plan.
- `docs/RHFS_MPI_hfsgg_optimization.md`: scope, diagnosis, and results for `HFSGG_MPI` optimization.
- `docs/1.0.1-dev.1_optimization_summary.md`, `docs/1.1.1-dev.1_optimization_summary.md`: staged optimization summaries.

## Contributor Notes

- Preserve the local Fortran style of the file being edited; avoid broad formatting-only rewrites.
- Keep interface/helper units named as `*_I.f90`, and shared modules aligned with `*_C.f90` or the local directory convention.
- Numerical changes should include input data, output differences, or regression notes.
- MPI changes should be checked against the corresponding serial path for input, output, and numerical consistency.

## License

This repository is distributed under the [MIT License](LICENSE).
