#!/usr/bin/env bash
set -euo pipefail
export GRASP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_directory="build"
cmake_args=(-DGRASP_ENABLE_MPI=ON)
while (($#)); do
    case "$1" in
        --debug) build_directory="build-debug"; cmake_args+=(-DCMAKE_BUILD_TYPE=Debug); shift ;;
        --pic) cmake_args+=(-DCMAKE_POSITION_INDEPENDENT_CODE=ON); shift ;;
        --no-mpi) cmake_args+=(-DGRASP_ENABLE_MPI=OFF); shift ;;
        --serial) cmake_args+=(-DGRASP_BUILD_SERIAL_APPS=ON); shift ;;
        --no-csfg) cmake_args+=(-DGRASP_BUILD_CSFG=OFF); shift ;;
        --no-openmp) cmake_args+=(-DGRASP_JJ2LSJ_OPENMP=OFF -DGRASP_RANGULAR_OPENMP=OFF); shift ;;
        --profile|--prefix|--build-dir)
            if (($# < 2)); then echo "Missing value for $1" >&2; exit 2; fi
            case "$1" in
                --profile) cmake_args+=("-DGRASP_FORTRAN_PROFILE=$2");;
                --prefix) cmake_args+=("-DCMAKE_INSTALL_PREFIX=$2");;
                --build-dir) build_directory="$2";;
            esac
            shift 2 ;;
        -D*) cmake_args+=("$1"); shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
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
    echo "ERROR: Build directory already exists: $build_abspath" >&2
    exit 1
fi
cmake -S "$GRASP" -B "$build_abspath" "${cmake_args[@]}"
cat <<EOF
Configured: $build_abspath
Build: cmake --build "$build_abspath" --parallel 4
Build/test validation: <shared-python> "$GRASP/test/f2018/check.py" "$build_abspath"
Production: f2018-release. Validation: f2018-check (Debug) or f2018-audit (Release).
Install only after validation: cmake --install "$build_abspath"
EOF
