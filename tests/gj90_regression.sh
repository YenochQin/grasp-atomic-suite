#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <gj90-bin> <data-dir> <case-name> <work-dir>" >&2
    exit 2
fi

gj90_bin=$1
data_dir=$2
case_name=$3
work_dir=$4

mkdir -p "${work_dir}"
rm -f "${work_dir}/${case_name}.gj" "${work_dir}/gj90.stdout" "${work_dir}/gj90.stdin"

cp "${data_dir}/isodata" "${work_dir}/"
cp "${data_dir}/${case_name}.c" "${work_dir}/"
cp "${data_dir}/${case_name}.m" "${work_dir}/"
cp "${data_dir}/${case_name}.w" "${work_dir}/"

cat > "${work_dir}/gj90.stdin" <<EOF
y
${case_name}
n
EOF

(
    cd "${work_dir}"
    "${gj90_bin}" < "gj90.stdin" > "gj90.stdout" 2>&1
)

if [ ! -f "${work_dir}/${case_name}.gj" ]; then
    echo "gj90 did not produce ${case_name}.gj" >&2
    exit 1
fi

cmp -s "${work_dir}/${case_name}.gj" "${data_dir}/${case_name}.gj" || {
    echo "generated ${case_name}.gj differs from reference output" >&2
    diff -u "${data_dir}/${case_name}.gj" "${work_dir}/${case_name}.gj" >&2 || true
    exit 1
}
