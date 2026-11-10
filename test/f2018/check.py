"""Fail-closed validation of the CURRENT target/test scope (not P9 certification)."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET

REQUIRED = Path(__file__).with_name("required.json")
SUPPORTED = {"f2018-check", "f2018-audit", "f2018-release"}
STRICT = {"f2018-check", "f2018-audit"}


def requirements(features: set[str]) -> tuple[set[str], set[str]]:
    targets: set[str] = set()
    tests: set[str] = set()
    for group in json.loads(REQUIRED.read_text())["groups"]:
        if set(group["requires"]) <= features:
            targets.update(group["targets"])
            tests.update(group["tests"])
    return targets, tests


def check_roster(features: set[str], targets: set[str], tests: set[str]) -> None:
    wanted_targets, wanted_tests = requirements(features)
    missing = {"targets": sorted(wanted_targets - targets), "tests": sorted(wanted_tests - tests)}
    if not tests or any(missing.values()):
        raise ValueError(f"Incomplete required scope: {missing}")
    # Disabled feature groups must not leave stale configured targets/tests behind.
    enabled_targets, enabled_tests = wanted_targets, wanted_tests
    for group in json.loads(REQUIRED.read_text())["groups"]:
        if not set(group["requires"]) <= features:
            extra = (set(group["targets"]) - enabled_targets) & targets
            extra_tests = (set(group["tests"]) - enabled_tests) & tests
            if extra or extra_tests:
                raise ValueError(f"Disabled group {group['name']} still present: {extra}, {extra_tests}")


def cache_values(build: Path) -> dict[str, str]:
    return dict(re.findall(r"^([^/#][^:\n]*):[^=\n]+=(.*)$", (build / "CMakeCache.txt").read_text(), re.M))


def check_commands(commands: list[dict], profile: str, compiler: Path) -> None:
    if profile not in SUPPORTED:
        raise ValueError(f"Unsupported/retired build profile: {profile}")
    if not commands:
        raise ValueError("Empty compile_commands.json")
    for entry in commands:
        if Path(entry["file"]).suffix.lower() not in {".f90", ".f", ".for", ".f95"}:
            continue
        args = entry.get("arguments") or shlex.split(entry["command"])
        if Path(args[0]).resolve() != compiler.resolve():
            raise ValueError(f"Foreign compiler for {entry['file']}: {args[0]}")
        required = set()
        if profile in STRICT:
            required |= {"-std=f2018", "-pedantic-errors", "-Wall", "-Wextra",
                         "-Werror=implicit-interface", "-Werror=implicit-procedure"}
        if profile == "f2018-check":
            required |= {"-O0", "-g", "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"}
        else:
            required.add("-O2")
        forbidden = [arg for arg in args if (
            arg in {"-fno-automatic", "-fallow-argument-mismatch", "-fall-intrinsics", "-ffast-math", "-Ofast", "-w"}
            or arg.startswith(("-fdefault-integer-", "-fdefault-real-", "-fdefault-double-", "-Wno-", "-fno-check"))
            or (arg.startswith("-std=") and arg != "-std=f2018")
            or (re.fullmatch(r"-O[0-9sgz]", arg) and arg not in required)
            or (profile == "f2018-check" and arg.startswith("-fcheck=") and arg != "-fcheck=all")
            or (profile == "f2018-check" and arg.startswith("-ffpe-trap=") and arg != "-ffpe-trap=invalid,zero,overflow")
        )]
        if required - set(args) or forbidden:
            raise ValueError(f"Invalid strict command for {entry['file']}: missing={required-set(args)}, forbidden={forbidden}")


def check_feature_commands(commands: list[dict], features: dict) -> None:
    openmp_targets = {"jj2lsj_conversion", "jj2lsj_thread_kernels", "jj2lsj_omp",
                      "test.jj2lsj_coefficients", "test.jj2lsj_ranking"}
    angular_targets = {"rangular_omp", "rangular_thread_state", "rangular_thread_kernels"}
    for entry in commands:
        if Path(entry["file"]).suffix.lower() not in {".f90", ".f", ".for", ".f95"}:
            continue
        if not features.get("csfg", False) and ("csfg" in entry["file"].lower() or "/src/tool/graspg/" in entry["file"]):
            raise ValueError(f"CSFG source leaked into the ordinary build: {entry['file']}")
        args = entry.get("arguments") or shlex.split(entry["command"])
        targets = set(re.findall(r"CMakeFiles/([^/]+)\.dir/", " ".join(args)))
        enabled = bool(targets & openmp_targets and features["openmp"] or
                       targets & angular_targets and features.get("rangular_openmp", False))
        if "-fopenmp" in args and not enabled:
            raise ValueError(f"OpenMP flag leaked into a disabled configuration: {entry['file']}")
        if enabled and "-fopenmp" not in args:
            raise ValueError(f"Missing OpenMP flag for {entry['file']}")


def check_junit(path: Path, expected: set[str]) -> dict[str, str]:
    cases = ET.parse(path).getroot().findall(".//testcase")
    status: dict[str, str] = {}
    for case in cases:
        name = case.attrib["name"]
        if name in status:
            raise ValueError(f"Duplicate test result: {name}")
        status[name] = "passed"
        for tag in ("skipped", "failure", "error"):
            if case.find(tag) is not None:
                status[name] = tag
    if not cases or set(status) != expected or any(value != "passed" for value in status.values()):
        raise ValueError(f"Incomplete/failed test run: missing={sorted(expected-set(status))}, "
                         f"unexpected={sorted(set(status)-expected)}, "
                         f"not_passed={ {k:v for k,v in status.items() if v != 'passed'} }")
    return status


def run_logged(command: list[str], log: Path, env: dict[str, str]) -> None:
    with log.open("w") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=env, check=False)
    if result.returncode:
        raise RuntimeError(f"Command exited {result.returncode}; see {log}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("--check-config", action="store_true")
    parser.add_argument("--features", default="")
    parser.add_argument("--inventory-only", action="store_true", help="Check registration and commands without running tests")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--warning-baseline", type=Path, help="Reviewed GNU 15 per-file warning budget")
    args = parser.parse_args()
    build = args.build.resolve()
    targets = set((build / "f2018-targets.txt").read_text().splitlines())
    if args.check_config:
        check_roster(set(filter(None, args.features.split(","))), targets,
                     set((build / "f2018-tests.txt").read_text().splitlines()))
        return 0
    audit = build / "audit" / "acceptance"
    audit.mkdir(parents=True, exist_ok=True)
    if not args.inventory_only:
        (audit / "ctest.xml").unlink(missing_ok=True)
    summary: dict = {"status": "running", "strict_acceptance": False}
    try:
        cache = cache_values(build)
        compiler = Path(cache["CMAKE_Fortran_COMPILER"])
        env = os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_DYNAMIC": "FALSE",
                            "FC": str(compiler), "OMPI_FC": str(compiler)}
        # Refresh CMake before inspecting inventories: cmake --build may otherwise
        # regenerate a changed tree after the preflight checked stale commands.
        run_logged(["cmake", "-S", cache["CMAKE_HOME_DIRECTORY"], "-B", str(build)], audit / "configure.log", env)
        targets = set((build / "f2018-targets.txt").read_text().splitlines())
        features = json.loads((build / "f2018-features.json").read_text())
        profile = features.pop("profile")
        if profile not in SUPPORTED:
            raise ValueError("Use f2018-release for production or f2018-check/audit for validation")
        cache = cache_values(build)
        version = subprocess.check_output([str(compiler), "-dumpfullversion"], text=True).strip()
        if version.split(".")[0] != "15":
            raise ValueError(f"GNU 15 required, found {version}")
        summary.update(profile=profile, compiler=str(compiler), compiler_version=version,
                       compiler_sha256=hashlib.sha256(compiler.read_bytes()).hexdigest(), features=features,
                       required_manifest_sha256=hashlib.sha256(REQUIRED.read_bytes()).hexdigest())
        source = Path(cache["CMAKE_HOME_DIRECTORY"])
        summary["source"] = {"directory": str(source)}
        revision = subprocess.run(["git", "-C", str(source), "rev-parse", "HEAD"], capture_output=True, text=True)
        if revision.returncode == 0:
            summary["source"]["revision"] = revision.stdout.strip()
            diff = subprocess.check_output(["git", "-C", str(source), "diff", "HEAD", "--", "."])
            summary["source"]["tracked_diff_sha256"] = hashlib.sha256(diff).hexdigest()
            names = subprocess.check_output(["git", "-C", str(source), "ls-files", "-co", "--exclude-standard", "-z"])
            summary["source"]["files"] = {name: hashlib.sha256((source/name).read_bytes()).hexdigest()
                for name in sorted(set(names.decode().split("\0"))) if name and (source/name).is_file()}
        summary["cmake"] = subprocess.check_output(["cmake", "--version"], text=True).splitlines()[0]
        summary["libraries"] = {key: value for key, value in cache.items()
            if key.startswith(("BLAS_", "LAPACK_", "MPI_")) and any(word in key for word in ("LIBRARY", "COMPILER", "VERSION"))}
        listing = subprocess.check_output(["ctest", "--test-dir", str(build), "--show-only=json-v1"], text=True, env=env)
        (audit / "registered-tests.json").write_text(listing)
        tests = {test["name"] for test in json.loads(listing)["tests"]}
        check_roster({k for k, v in features.items() if v}, targets, tests)
        commands = json.loads((build / "compile_commands.json").read_text())
        check_commands(commands, profile, compiler)
        check_feature_commands(commands, features)
        compiled_targets = set()
        for entry in commands:
            command = " ".join(entry.get("arguments", [])) or entry.get("command", "")
            compiled_targets.update(re.findall(r"CMakeFiles/([^/]+)\.dir/", command))
        missing_compilation = requirements({k for k,v in features.items() if v})[0] - compiled_targets
        if missing_compilation:
            raise ValueError(f"Required targets have no compilation entries: {sorted(missing_compilation)}")
        summary.update(test_count=len(tests), compile_count=len(commands))
        if args.inventory_only:
            summary["status"] = "inventory_only"
        else:
            run_logged(["cmake", "--build", str(build), "--parallel", str(args.jobs)], audit / "build.log", env)
            if profile in STRICT:
                from warning_inventory import BASELINE, compare, inventory, write_json
                observation = inventory(build)
                write_json(build / "audit/warnings.json", observation)
                baseline_path = args.warning_baseline or BASELINE
                summary["warning_baseline_sha256"] = hashlib.sha256(baseline_path.read_bytes()).hexdigest()
                summary["warning_gate"] = compare(observation, json.loads(baseline_path.read_text()))
            else:
                summary["warning_gate"] = {"status": "not_collected", "reason": "Production flags; validate warnings in f2018-audit"}
            junit = audit / "ctest.xml"
            junit.unlink(missing_ok=True)
            # Parse results even when CTest returns failure, preserving skip information.
            with (audit / "ctest.log").open("w") as log:
                result = subprocess.run(["ctest", "--test-dir", str(build), "--output-on-failure",
                                         "--output-junit", str(junit)], stdout=log, stderr=subprocess.STDOUT, env=env)
            summary["tests"] = check_junit(junit, tests)
            if result.returncode:
                raise RuntimeError(f"CTest exited {result.returncode}")
            summary.update(status="passed", scope_acceptance=True, strict_acceptance=profile in STRICT)
        print(f"{summary['status']}: {profile}, {len(tests)} registered tests; strict_acceptance={summary['strict_acceptance']}")
        return 0
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError, ET.ParseError) as error:
        summary.update(status="failed", error=str(error))
        print(str(error), file=sys.stderr)
        return 1
    finally:
        junit = audit / "ctest.xml"
        if junit.is_file() and not args.inventory_only:
            try:
                summary["test_outcomes"] = {case.attrib["name"]:
                    next((tag for tag in ("skipped", "failure", "error") if case.find(tag) is not None), "passed")
                    for case in ET.parse(junit).getroot().findall(".//testcase")}
            except ET.ParseError:
                pass
        (audit / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    sys.exit(main())
