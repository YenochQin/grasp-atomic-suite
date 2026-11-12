#!/usr/bin/env bash
set -euo pipefail
export GRASP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_directory="build"
# The default command configures every retained ordinary application, including
# rwfnestimate. Explicit options below may select a smaller or extended scope.
cmake_args=(-DGRASP_BUILD_SERIAL_APPS=ON -DGRASP_ENABLE_MPI=ON
            -DGRASP_JJ2LSJ_OPENMP=ON -DGRASP_RANGULAR_OPENMP=ON
            -DGRASP_BUILD_CSFG=OFF)

# Quote command arguments for the fish commands displayed to the user.
fish_quote() {
    local value=$1
    value=${value//\\/\\\\}
    value=${value//\'/\\\'}
    printf "'%s'" "$value"
}

cache_value() {
    sed -n "s/^${1}:[^=]*=//p" "$build_abspath/CMakeCache.txt"
}

show_next_steps() {
    local build_arg python_executable install_prefix rangular_openmp=OFF
    local build_scope="CUSTOM: compile only the selected program groups listed below"
    local compile_description="Compile the selected program groups"
    local rwfnestimate_status="not selected"
    build_arg=$(fish_quote "$build_abspath")
    python_executable=$(cache_value Python3_EXECUTABLE)
    if [[ -z "$python_executable" ]]; then
        python_executable=$(cache_value _Python3_EXECUTABLE)
    fi
    install_prefix=$(cache_value CMAKE_INSTALL_PREFIX)
    if [[ -f "$build_abspath/f2018-targets.txt" ]] &&
       grep -qx rangular_omp "$build_abspath/f2018-targets.txt"; then
        rangular_openmp=ON
    fi
    if grep -qx rwfnestimate "$build_abspath/f2018-targets.txt"; then
        rwfnestimate_status="included"
    fi
    if [[ $(cache_value GRASP_BUILD_SERIAL_APPS) == ON &&
          $(cache_value GRASP_ENABLE_MPI) == ON &&
          $(cache_value GRASP_JJ2LSJ_OPENMP) == ON &&
          "$rangular_openmp" == ON &&
          $(cache_value GRASP_BUILD_CSFG) == OFF ]]; then
        build_scope="ALL NON-CSFG PROGRAMS: serial + MPI + OpenMP"
        compile_description="Compile ALL non-CSFG programs, including rwfnestimate"
    fi
    cat <<EOF

Build configuration
  Source:          $GRASP
  Build directory: $build_abspath
  Profile:         $(cache_value GRASP_FORTRAN_PROFILE) ($(cache_value CMAKE_BUILD_TYPE))
  Fortran:         $(cache_value CMAKE_Fortran_COMPILER)
  Build scope:     $build_scope
  rwfnestimate:    $rwfnestimate_status
  Serial applications:            $(cache_value GRASP_BUILD_SERIAL_APPS)
  MPI applications:               $(cache_value GRASP_ENABLE_MPI)
  OpenMP: jj2lsj=$(cache_value GRASP_JJ2LSJ_OPENMP), rangular=$rangular_openmp
  CSFG extension:                 $(cache_value GRASP_BUILD_CSFG)

Output locations
  Executable programs: $build_abspath/bin
  Libraries:           $build_abspath/lib
  Optional install:    $install_prefix

Next commands (fish; usable from any working directory)

1. $compile_description:
   cmake --build $build_arg --parallel 4

2. Validate the build and run all registered tests:
   $(fish_quote "$python_executable") $(fish_quote "$GRASP/test/f2018/check.py") $build_arg --jobs 4
   This also performs an incremental build and checks compiler flags and targets.
   Result:   $build_abspath/audit/acceptance/summary.json
   Test log: $build_abspath/audit/acceptance/ctest.log

   To rerun only the tests after building:
   ctest --test-dir $build_arg --output-on-failure

3. Point generated Tools jobs to the built programs:
   set -gx SUITE_BINDIR $(fish_quote "$build_abspath/bin")
   For Tools, server_settings.suite_bindir can point to this same directory.

4. Optional: after validation passes, install to $install_prefix:
   cmake --install $build_arg
   Installed programs: $install_prefix/bin
   To use the installed programs with Tools:
   set -gx SUITE_BINDIR $(fish_quote "$install_prefix/bin")
EOF
}

show_help() {
    local python_option=""
    if [[ -x "$GRASP/../graspkit-tools/.venv/bin/python" ]]; then
        local shared_python
        shared_python="$(cd "$GRASP/../graspkit-tools/.venv/bin" && pwd)/python"
        python_option=" -DPython3_EXECUTABLE=$(fish_quote "$shared_python")"
    fi
    cat <<EOF
Configure ALL non-CSFG programs: serial, MPI and both OpenMP targets.
This includes rwfnestimate. Compilation is the next step shown after configuration.

Usage: ./configure.sh [options] [-DNAME=VALUE ...]
  --serial          Enable serial applications (already enabled by default)
  --no-serial       Explicitly omit ordinary serial applications, including rwfnestimate
  --no-mpi          Disable MPI programs (enabled by default)
  --no-openmp       Disable both jj2lsj and rangular OpenMP support
  --no-csfg         Disable the optional CSFG extension (already off by default)
  --debug           Use a Debug build; default directory becomes build-debug
  --profile NAME    f2018-release (production), f2018-check or f2018-audit
  --build-dir PATH  Build directory (default: build); compatible caches are updated
  --prefix PATH     Optional installation prefix (default: <build>/install)
  --pic             Enable position-independent code
  -h, --help        Show this help

All ordinary programs, with MPI and both OpenMP targets (fish):
  set -gx FC /usr/bin/gfortran-15
  set -gx OMPI_FC /usr/bin/gfortran-15
  ./configure.sh$python_option

GNU Fortran 15, CMake >= 3.31, BLAS/LAPACK and Python are required.
MPI programs also require a compatible MPI Fortran environment.
Load the compiler/MPI environment before configuring.
Use --debug with --profile f2018-check. Release/audit use Release builds.
The completion message reports the selected options and exact next commands.
Reusing a build directory applies these defaults again; repeat any opt-out options.
Use a new directory when changing the source tree or compiler.
EOF
}

while (($#)); do
    case "$1" in
        -h|--help) show_help; exit 0 ;;
        --debug) build_directory="build-debug"; cmake_args+=(-DCMAKE_BUILD_TYPE=Debug); shift ;;
        --pic) cmake_args+=(-DCMAKE_POSITION_INDEPENDENT_CODE=ON); shift ;;
        --no-mpi) cmake_args+=(-DGRASP_ENABLE_MPI=OFF); shift ;;
        --serial) cmake_args+=(-DGRASP_BUILD_SERIAL_APPS=ON); shift ;;
        --no-serial) cmake_args+=(-DGRASP_BUILD_SERIAL_APPS=OFF); shift ;;
        --no-csfg) cmake_args+=(-DGRASP_BUILD_CSFG=OFF); shift ;;
        --no-openmp) cmake_args+=(-DGRASP_JJ2LSJ_OPENMP=OFF -DGRASP_RANGULAR_OPENMP=OFF); shift ;;
        --profile|--prefix|--build-dir)
            if (($# < 2)) || [[ -z "$2" ]]; then echo "Missing value for $1; see ./configure.sh --help" >&2; exit 2; fi
            case "$1" in
                --profile) cmake_args+=("-DGRASP_FORTRAN_PROFILE=$2");;
                --prefix) cmake_args+=("-DCMAKE_INSTALL_PREFIX=$2");;
                --build-dir) build_directory="$2";;
            esac
            shift 2 ;;
        -D*) cmake_args+=("$1"); shift ;;
        *) echo "Unknown option: $1; see ./configure.sh --help" >&2; exit 2 ;;
    esac
done
if [[ -n ${FC:-} ]]; then
    cmake_args+=("-DCMAKE_Fortran_COMPILER=$FC")
    export OMPI_FC="${OMPI_FC:-$FC}"
fi
if [[ "$build_directory" = /* ]]; then
    build_abspath="$build_directory"
else
    build_abspath="$GRASP/$build_directory"
fi
if [[ -e "$build_abspath" ]]; then
    if [[ ! -f "$build_abspath/CMakeCache.txt" ]] ||
       [[ $(cache_value CMAKE_HOME_DIRECTORY) != "$GRASP" ]]; then
        echo "ERROR: Existing directory is not a configured build of this source tree: $build_abspath" >&2
        echo "Choose a new directory with --build-dir PATH." >&2
        exit 1
    fi
    # Changing compilers in an existing cache can mix incompatible objects and
    # modules. CMake retains the cached compiler when no override is supplied.
    requested_compiler=""
    for argument in "${cmake_args[@]}"; do
        case "$argument" in
            -DCMAKE_Fortran_COMPILER=*|-DCMAKE_Fortran_COMPILER:*=*)
                requested_compiler=${argument#*=} ;;
        esac
    done
    if [[ -n "$requested_compiler" ]]; then
        resolved_compiler=$(command -v "$requested_compiler" || true)
        if [[ -z "$resolved_compiler" ]] ||
           [[ ! "$resolved_compiler" -ef "$(cache_value CMAKE_Fortran_COMPILER)" ]]; then
            echo "ERROR: Requested Fortran compiler differs from the existing build: $requested_compiler" >&2
            echo "Choose a new directory with --build-dir PATH when changing compilers." >&2
            exit 1
        fi
    fi
    echo "Updating the existing build configuration: $build_abspath"
fi
cmake -S "$GRASP" -B "$build_abspath" "${cmake_args[@]}"
build_abspath=$(cd "$build_abspath" && pwd)
echo "Configuration complete. Continue with the build and validation commands below."
show_next_steps
