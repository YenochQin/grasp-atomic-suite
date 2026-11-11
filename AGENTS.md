# Repository Guidelines

## Shell commands

All shell commands shown to the user, including replies, documentation, and copyable examples, must use fish shell syntax. Label shell code blocks as `fish`; do not output Bash commands or Bash-specific syntax. Use fish-compatible variable assignments, environment variables, control flow, and virtual-environment activation.

## Project Structure & Module Organization
Core source lives under `src/`. Application programs are in `src/appl/<program>/`, and shared libraries are in `src/lib/<library>/`. Build outputs are installed to `bin/` and `lib/`. Development notes and numerical reports live under `docs/`, and sample validation data lives under `data/`.

Use [docs/README.md](docs/README.md) as the documentation status index. Current
implementation guides belong in `docs/implemented/`, pending experiments in
`docs/plans/`, completed reviews in `docs/reviews/`, theory drafts in
`docs/reference/`, and superseded plans or historical reports in `docs/archive/`.
State the source revision and distinguish implemented code from completed
physical validation; update relative links and the index when moving documents.

The supported ordinary RMCDHF programs are `rmcdhf{,_mem,_mpi,_mem_mpi}`;
CSFG remains behind its separate opt-in boundary. The rejected `rmcdhf_orbopt*`,
`rmcdhf_mem_strict_mpi`, and `rmcdhf_mem_diag_mpi` programs are retired; historical
documents are not authorization to restore them. See [retirement scope](docs/implemented/rejected_rmcdhf_variants.md).
Retained programs share
`libmod`, `lib9290`, `libdvd90`, and `mpi90` with the property programs; do not
introduce a second copy of these libraries. Preserve suite MPI helpers.
Keep `src/lib/libmod/parameter_def_M.f90` as the legacy re-export interface.
Define shared capacities only in `suite_parameters_M.f90`, using upstream defaults:
`KEYORB=215`, `NNNP=590`, `NNN1=NNNP+10`, `NNNW=127`, `NNNWM1=NNNW-1`, and
`NNNWM2=NNNW-2`. The grid patcher edits this configuration module for capacity
changes; existing applications continue to `USE parameter_def`.
Each grid entry point registered in the patcher keeps its upstream point/finite-nucleus
settings and ACCY initialization order, including RDENSITY. Finite nuclei default to `N=NNNP=590`.
Upstream RMCDHF computes ACCY before interactive grid input; explicitly supply
ACCY when changing H in restart inputs. Grid experiments belong in separate copies.

## Build, Test, and Development Commands

Use **GNU Fortran 15** for builds and validation. Set Open MPI's `OMPI_FC`
explicitly to `/usr/bin/gfortran-15`. Initialize any compiler, MPI, and BLAS
modules first, then set both exported compiler variables as shown below. Verify that
`/usr/bin/gfortran-15 -dumpversion` reports major version 15 and
`mpifort --showme:command` reports `/usr/bin/gfortran-15`.
For direct CMake configuration, pass `-DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15`;
use a fresh build directory if an existing cache selected another compiler.

Prefer the out-of-source CMake flow:

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
./configure.sh
cmake --build build --target install -j4
```

Use `./configure.sh --debug` for a debug build in `build-debug/`. Use `./scripts/clean-build-artifacts.sh` to remove CMake build directories and installed build outputs.

Enable `-DGRASP_BUILD_SERIAL_APPS=ON` to build the original serial RMCDHF programs
alongside serial and MPI property programs. In this workspace, run Python
scripts/tests with `../graspkit-tools/.venv/bin/python`; do not create another
environment or run `uv sync`/`uv run` here.

For grid or capacity changes, use [common parameters](docs/implemented/common_parameters.md)
for the current contract and [scripts/README.md](scripts/README.md) for patch-tool
configuration, preview/apply, and restore commands. This partial tree needs the
`atomic-suite` layout; `grasp2018` requires a full upstream tree. Preserve omitted
settings and keep local config ignored. Use separate source copies per parameter
set; compilation is a separate step from applying the patch.

## Coding Style & Naming Conventions
Match the existing Fortran style in the file you are editing. The tree contains both older fixed-format-style indentation and newer free-form modules, so preserve local formatting instead of normalizing whole files. Keep filenames lowercase and aligned with existing suffix conventions such as `*_I.f90` for interface/helper units and `*_C.f90` for shared code modules. Ordinary program directories and targets use lowercase functional names like `rangular`, `rhfs`, and `gj`; use the build profile to identify the Fortran standard. For renamed paths and executable names, see [program naming](docs/implemented/program_naming.md).

## Testing Guidelines
CI uses the CMake workflow. Run tests covering the affected programs, using a
CTest name filter when appropriate:

```fish
ctest --test-dir build --output-on-failure
```

Broaden validation for shared parameters, libraries, numerical routines, or build
changes. Follow the common-parameter guide when changing array capacities or grid
defaults; a passing build does not establish physical convergence. For documentation
changes, check the relevant references and examples. Keep functional tests close to
the affected program or under `tests/`, and document required inputs.

## Commit & Pull Request Guidelines
Recent history favors short, imperative commit subjects, optionally with an issue reference, for example `Update README.md` or `Bug fix in librang (#112)`. Keep commits scoped to one change. Pull requests should explain the scientific or numerical impact, list the build/test commands you ran, link related issues, and include sample output or input-file changes when behavior or results change.

## Configuration Tips
Use CMake cache variables, compiler environment variables, or `CMakeLists.user` for local build customization. The required `/usr/bin/gfortran-15` path above is the repository toolchain contract; keep other machine-local settings and generated build directories out of commits.
