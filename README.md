# grasp-atomic-suite

[中文](README_ZH.md) | English

`grasp-atomic-suite` is a Fortran research and development repository for relativistic atomic-structure and atomic-property calculations. It follows the GRASP92/GRASP2018-style source layout and currently focuses on `rhfs90`, `ris4`, `gj90`, and their MPI support for hyperfine-structure, isotope-shift, and Landé `g_J` calculations.

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
├── bin/                      # installed executables
└── lib/                      # installed static libraries and Fortran modules
```

## Programs

| Program | Source directory | CMake targets | Purpose |
| --- | --- | --- | --- |
| `gj90` | `src/appl/gj90` | `gj90`, `gj90_mpi` | Standalone Landé `g_J` calculator extracted from the RHFS path. It reads `isodata`, `name.c`, `name.m/name.cm`, and `name.w`, then writes `name.gj/name.cgj`. |
| `rhfs90` | `src/appl/rhfs90` | `rhfs`, `rhfs_mpi` | Relativistic hyperfine-structure program. It computes hyperfine constants and related matrix elements, writing `name.h/name.ch` and `name.hoffd/name.choffd`. |
| `ris4` | `src/appl/ris4` | `ris4`, `ris4_mpi` | Relativistic isotope-shift program. It computes normal mass shift, specific mass shift, and field-shift electronic factors, writing `name.i/name.ci` and angular intermediate data. |

The default CMake configuration attempts to build MPI versions. If no MPI Fortran toolchain is found, MPI targets are skipped. Serial application targets are disabled by default and must be enabled explicitly at CMake configuration time.

## Build

The recommended build flow is out-of-source CMake:

```sh
./configure.sh
cmake --build build --target install -j4
```

`configure.sh` passes `-DGRASP_ENABLE_MPI=ON` by default. If MPI Fortran is available, installation normally produces `bin/gj90_mpi`, `bin/rhfs_mpi`, `bin/ris4_mpi`, and the associated libraries.

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

The applications use GRASP-style state-name input. After you enter `name`, the program looks for matching files in the current working directory.

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

These files can be used to compare `gj90` Landé-factor output. The current `tests/` directory is empty; if CTest scripts are added later, run them from the build directory:

```sh
ctest
```

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
