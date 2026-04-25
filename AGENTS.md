# Repository Guidelines

## Project Structure & Module Organization
Core source lives under `src/`. Application programs are in `src/appl/<program>/`, and shared libraries are in `src/lib/<library>/`. Build outputs are installed to `bin/` and `lib/`. Development notes and numerical reports live under `docs/`, and sample validation data lives under `data/`.

## Build, Test, and Development Commands
Prefer the out-of-source CMake flow:

```sh
./configure.sh
cmake --build build --target install -j4
```

Use `./configure.sh --debug` for a debug build in `build-debug/`. Use `./scripts/clean-build-artifacts.sh` to remove CMake build directories and installed build outputs.

## Coding Style & Naming Conventions
Match the existing Fortran style in the file you are editing. The tree contains both older fixed-format-style indentation and newer free-form modules, so preserve local formatting instead of normalizing whole files. Keep filenames lowercase and aligned with existing suffix conventions such as `*_I.f90` for interface/helper units and `*_C.f90` for shared code modules. Program directories and targets use lowercase names like `rhfs90`, `rdensity`, and `ris4`.

## Testing Guidelines
CI uses the CMake workflow. If your checkout includes CTest tests, run:

```sh
cd build && ctest
```

For functional verification, keep new tests close to the affected program or under `tests/`, and document required input data in that directory.

## Commit & Pull Request Guidelines
Recent history favors short, imperative commit subjects, optionally with an issue reference, for example `Update README.md` or `Bug fix in librang (#112)`. Keep commits scoped to one change. Pull requests should explain the scientific or numerical impact, list the build/test commands you ran, link related issues, and include sample output or input-file changes when behavior or results change.

## Configuration Tips
Use CMake cache variables, compiler environment variables, or `CMakeLists.user` for local build customization. Do not commit machine-local compiler paths or generated build directories.
