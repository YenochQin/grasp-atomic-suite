#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 7 ]; then
    echo "usage: $0 <mpiexec> <numproc-flag> <nprocs> <gj90-mpi-bin> <data-dir> <case-name> <work-dir>" >&2
    exit 2
fi

mpiexec_bin=$1
numproc_flag=$2
nprocs=$3
gj90_mpi_bin=$4
data_dir=$5
case_name=$6
work_dir=$7

mkdir -p "${work_dir}"
rm -f "${work_dir}/${case_name}.gj" "${work_dir}/gj90_mpi.stdout"

cp "${data_dir}/isodata" "${work_dir}/"
cp "${data_dir}/${case_name}.c" "${work_dir}/"
cp "${data_dir}/${case_name}.m" "${work_dir}/"
cp "${data_dir}/${case_name}.w" "${work_dir}/"

(
    cd "${work_dir}"
    "${mpiexec_bin}" "${numproc_flag}" "${nprocs}" "${gj90_mpi_bin}" "${case_name}" > "gj90_mpi.stdout" 2>&1
)

if [ ! -f "${work_dir}/${case_name}.gj" ]; then
    echo "gj90_mpi did not produce ${case_name}.gj" >&2
    exit 1
fi

cmp -s "${work_dir}/${case_name}.gj" "${data_dir}/${case_name}.gj" || {
    echo "generated ${case_name}.gj differs from reference output" >&2
    diff -u "${data_dir}/${case_name}.gj" "${work_dir}/${case_name}.gj" >&2 || true
    exit 1
}
