# Repository Guidelines

## Project Structure & Module Organization
Core source lives under `src/`. Application programs are in `src/appl/<program>/`, shared libraries are in `src/lib/<library>/`, and extra utilities are in `src/tool/`. Build outputs are installed to `bin/` and `lib/`. Example and regression-style runs live under `grasptest/<case-or-example>/script/`, with program-specific notes in files such as `src/appl/rdensity/README.md` and `src/appl/ris4/README.md`.

## Build, Test, and Development Commands
Prefer the out-of-source CMake flow:

```sh
./configure.sh
cd build && make -j4 install
```

Use `./configure.sh --debug` for a debug build in `build-debug/`. The legacy in-tree build still works:

```sh
make
make src/lib/libmod
make src/appl/rhfs90
```

Use `make clean` to remove object files and `make cleanall` to also clear installed binaries and archives.

## Coding Style & Naming Conventions
Match the existing Fortran style in the file you are editing. The tree contains both older fixed-format-style indentation and newer free-form modules, so preserve local formatting instead of normalizing whole files. Keep filenames lowercase and aligned with existing suffix conventions such as `*_I.f90` for interface/helper units and `*_C.f90` for shared code modules. Program directories and targets use lowercase names like `rhfs90`, `rdensity`, and `ris4`.

## Testing Guidelines
CI builds both CMake and legacy Makefile workflows. If your checkout includes `test/`, run:

```sh
cd build && ctest
```

For functional verification, run the scripts in `grasptest/`, for example `grasptest/example1/script/script_ex1` or `grasptest/case1/script/sh_case1`. Keep new tests close to the affected program and document required input data in that directory.

## Commit & Pull Request Guidelines
Recent history favors short, imperative commit subjects, optionally with an issue reference, for example `Update README.md` or `Bug fix in librang (#112)`. Keep commits scoped to one change. Pull requests should explain the scientific or numerical impact, list the build/test commands you ran, link related issues, and include sample output or input-file changes when behavior or results change.

## Configuration Tips
Do not edit the top-level `Makefile` for local compiler settings. Copy `Make.user.gfortran` or `Make.user.ifort` to `Make.user` and override `FC`, `FC_MPI`, `FC_FLAGS`, and linker paths there.
