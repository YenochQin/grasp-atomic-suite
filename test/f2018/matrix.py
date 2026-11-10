"""Build ordinary-GRASP production and validation profiles; CSFG is isolated.

Use --stage build to collect a fresh matrix before explicitly creating/reviewing
a warning baseline. The default stage runs the complete CTest suite and checks
the reviewed warning budget for every configuration; partial selections never
claim full-matrix acceptance.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import product
import json
import os
from pathlib import Path
import subprocess
import sys

from check import run_logged
from warning_inventory import AXES, BASELINE, case_id, compare, digest, inventory, write_json

SOURCE = Path(__file__).resolve().parents[2]
OPTIONS = {"serial": "GRASP_BUILD_SERIAL_APPS", "mpi": "GRASP_ENABLE_MPI",
           "csfg": "GRASP_BUILD_CSFG", "openmp": "GRASP_JJ2LSJ_OPENMP"}


def cases() -> list[dict]:
    return [{"id": case_id(profile, features), "profile": profile, "features": features}
            for profile in ("f2018-check", "f2018-audit", "f2018-release")
            for serial, mpi, openmp in product((True, False), repeat=3)
            for angular_openmp in ((True, False) if serial else (False,))
            for features in [{"serial": serial, "mpi": mpi, "csfg": False, "openmp": openmp,
                              "rangular_openmp": angular_openmp}]]


def configure_command(case: dict, build: Path, compiler: Path, mpi_compiler: Path, mpiexec: Path) -> list[str]:
    args = ["cmake", "-S", str(SOURCE), "-B", str(build),
            f"-DCMAKE_Fortran_COMPILER={compiler}", f"-DPython3_EXECUTABLE={sys.executable}",
            f"-DGRASP_FORTRAN_PROFILE={case['profile']}",
            "-DGRASP_CAPTURE_WARNINGS=" + ("OFF" if case["profile"] == "f2018-release" else "ON"),
            "-DCMAKE_BUILD_TYPE=" + ("Debug" if case["profile"] == "f2018-check" else "Release")]
    args += [f"-D{OPTIONS[axis]}={'ON' if case['features'][axis] else 'OFF'}" for axis in AXES]
    args += [f"-DGRASP_RANGULAR_OPENMP={'ON' if case['features']['rangular_openmp'] else 'OFF'}"]
    if case["features"]["mpi"]:
        args += [f"-DMPI_Fortran_COMPILER={mpi_compiler}", f"-DMPIEXEC_EXECUTABLE={mpiexec}"]
    return args


def run_case(case: dict, args: argparse.Namespace, baseline: dict | None) -> dict:
    build = args.root / case["id"]
    build.mkdir(parents=True, exist_ok=True)
    audit = build / "audit"
    audit.mkdir(exist_ok=True)
    result = {**case, "status": "running", "strict_acceptance": False}
    env = os.environ | {"FC": str(args.compiler), "OMPI_FC": str(args.compiler), "LC_ALL": "C",
                        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_DYNAMIC": "FALSE"}
    try:
        run_logged(configure_command(case, build, args.compiler, args.mpi_compiler, args.mpiexec),
                   audit / "matrix-configure.log", env)
        features = json.loads((build / "f2018-features.json").read_text())
        if features.pop("profile") != case["profile"] or features != case["features"]:
            raise ValueError("Requested feature matrix differs from actual configured features")
        command = [sys.executable, str(SOURCE / "test/f2018/check.py"), str(build), "--jobs", str(args.jobs)]
        if args.stage == "build":
            command.append("--inventory-only")
        else:
            command += ["--warning-baseline", str(args.baseline.resolve())]
        run_logged(command, audit / "matrix-check.log", env)
        if args.stage == "build":
            run_logged(["cmake", "--build", str(build), "--parallel", str(args.jobs)], audit / "matrix-build.log", env)
        observation = None
        if case["profile"] != "f2018-release":
            observation = inventory(build)
            write_json(audit / "warnings.json", observation)
            result.update(compilation_count=observation["compilation_count"],
                          source_file_count=len(observation["files"]),
                          warning_count=sum(observation["counts"].values()),
                          category_counts=observation["category_counts"])
        else:
            result["warning_gate"] = "not_collected_in_production"
        summary = json.loads((audit / "acceptance/summary.json").read_text())
        result["test_count"] = summary["test_count"]
        result["acceptance_summary_sha256"] = digest(audit / "acceptance/summary.json")
        if args.stage == "validate":
            if not summary.get("scope_acceptance"):
                raise ValueError("Configuration did not pass its build/test acceptance")
            if observation is not None:
                result["warnings"] = compare(observation, baseline)
            result.update(status="passed", scope_acceptance=True, strict_acceptance=summary["strict_acceptance"])
        else:
            result["status"] = "built_not_tested"
    except (ValueError, KeyError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        result.update(status="failed", error=str(error))
    write_json(audit / "matrix-case.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--stage", choices=("build", "validate"), default="validate")
    parser.add_argument("--compiler", type=Path, default=Path("/usr/bin/gfortran-15"))
    parser.add_argument("--mpi-compiler", type=Path, default=Path("/usr/bin/mpifort"))
    parser.add_argument("--mpiexec", type=Path, default=Path("/usr/bin/mpiexec"))
    parser.add_argument("--jobs", type=int, default=3, help="Build parallelism per configuration")
    parser.add_argument("--parallel-configs", type=int, default=2)
    parser.add_argument("--case", action="append", choices=[case["id"] for case in cases()])
    parser.add_argument("--baseline", type=Path, default=BASELINE)
    args = parser.parse_args()
    if args.jobs < 1 or args.parallel_configs < 1:
        parser.error("Parallelism must be positive")
    args.root = args.root.resolve()
    args.root.mkdir(parents=True, exist_ok=True)
    # Keep the versioned compiler pathname: resolving its symlink can change
    # wrapper/toolchain identification on distributions with parallel GCCs.
    args.compiler = args.compiler.absolute()
    baseline = json.loads(args.baseline.read_text()) if args.stage == "validate" else None
    selected = [case for case in cases() if not args.case or case["id"] in args.case]
    outcomes = {}
    summary = {"schema_version": 1, "stage": args.stage, "full_matrix_acceptance": False,
               "expected_cases": [case["id"] for case in cases()], "cases": outcomes}
    write_json(args.root / f"matrix-{args.stage}.json", summary)
    with ThreadPoolExecutor(max_workers=args.parallel_configs) as pool:
        pending = {pool.submit(run_case, case, args, baseline): case["id"] for case in selected}
        for future in as_completed(pending):
            result = future.result()
            outcomes[result["id"]] = result
            print(f"{len(outcomes)}/{len(selected)} {result['id']}: {result['status']}"
                  + (f" — {result['error']}" if "error" in result else ""), flush=True)
            write_json(args.root / f"matrix-{args.stage}.json", summary)
    passed = all(result["status"] != "failed" for result in outcomes.values())
    summary["full_matrix_acceptance"] = (args.stage == "validate" and passed and len(outcomes) == len(cases()))
    summary["status"] = "passed" if passed else "failed"
    if baseline is not None:
        summary["baseline_sha256"] = digest(args.baseline)
    summary["cases"] = dict(sorted(outcomes.items()))
    write_json(args.root / f"matrix-{args.stage}.json", summary)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
