# Repository Guidelines

## Project Structure & Module Organization
Core source lives under `src/`. Application programs are in `src/appl/<program>/`, and shared libraries are in `src/lib/<library>/`. Build outputs are installed to `bin/` and `lib/`. Development notes and numerical reports live under `docs/`, and sample validation data lives under `data/`.

Use [docs/README.md](docs/README.md) as the documentation status index. Current
implementation guides belong in `docs/implemented/`, pending experiments in
`docs/plans/`, completed reviews in `docs/reviews/`, theory drafts in
`docs/reference/`, and superseded plans or historical reports in `docs/archive/`.
State the source revision and distinguish implemented code from completed
physical validation; update relative links and the index when moving documents.

The `rmcdhf90*` programs migrated from `rmcdhf_test` build as
`rmcdhf_orbopt`, `rmcdhf_orbopt_mem`, and their `_mpi` variants. They share
`libmod`, `lib9290`, `libdvd90`, and `mpi90` with the property programs; do not
introduce a second copy of these libraries. Preserve suite MPI helpers.
Shared capacities and numerical defaults live only in
`src/lib/libmod/suite_parameters_M.f90`; `parameter_def_M.f90` re-exports legacy
names and must not redefine them. `radial_grid_defaults_M.f90` initializes and
validates the legacy runtime state for all six entry points. All finite-nucleus
programs default to `N=NNNP` (2990); point-nucleus defaults remain separate.
Update ACCY after interactive grid input and before explicit ACCY overrides.

## Build, Test, and Development Commands
Prefer the out-of-source CMake flow:

```sh
./configure.sh
cmake --build build --target install -j4
```

Use `./configure.sh --debug` for a debug build in `build-debug/`. Use `./scripts/clean-build-artifacts.sh` to remove CMake build directories and installed build outputs.

Enable `-DGRASP_BUILD_SERIAL_APPS=ON` to build all four orbital programs
alongside serial and MPI property programs. In this workspace, run Python
scripts/tests with `../graspkit-tools/.venv/bin/python`; do not create another
environment or run `uv sync`/`uv run` here. On cluster hosts initialize the
compiler, MPI, and BLAS modules before configuration.

For grid or capacity changes, use [common parameters](docs/implemented/common_parameters.md)
for the current contract and [scripts/README.md](scripts/README.md) for patch-tool
configuration, preview/apply, and restore commands. This partial tree needs the
`atomic-suite` layout; `grasp2018` requires a full upstream tree. Preserve omitted
settings and keep local config ignored. Use separate source copies per parameter
set; compilation is a separate step from applying the patch.

## Coding Style & Naming Conventions
Match the existing Fortran style in the file you are editing. The tree contains both older fixed-format-style indentation and newer free-form modules, so preserve local formatting instead of normalizing whole files. Keep filenames lowercase and aligned with existing suffix conventions such as `*_I.f90` for interface/helper units and `*_C.f90` for shared code modules. Program directories and targets use lowercase names like `rhfs90`, `rdensity`, and `ris4`.

## Testing Guidelines
CI uses the CMake workflow. Run tests covering the affected programs, using a
CTest name filter when appropriate:

```sh
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
Use CMake cache variables, compiler environment variables, or `CMakeLists.user` for local build customization. Do not commit machine-local compiler paths or generated build directories.
