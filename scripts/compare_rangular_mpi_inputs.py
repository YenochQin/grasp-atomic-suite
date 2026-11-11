#!/usr/bin/env python3
"""Compare text and descriptor input in the same rangular_mpi executable."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from typing import Any

from compare_rangular_versions import THREAD_ENVIRONMENT, file_identity, save_json, stop_process
from verify_csf_descriptors import verify_descriptors


def read_rank_timings(directory: Path, ranks: int) -> list[dict[str, float]]:
    profiles = []
    for rank_index in range(ranks):
        profile = directory / f"{rank_index:03d}" / f"timings.rank{rank_index:03d}.csv"
        with profile.open() as stream:
            rows = list(csv.DictReader(stream))
        if not rows or rows[-1]["stage"] != "program_before_finalize_including_interaction_input":
            raise ValueError(f"Incomplete rank timings: {profile}")
        totals: dict[str, float] = {}
        for row in rows:
            elapsed = float(row["wall_seconds"])
            if not math.isfinite(elapsed) or elapsed < 0:
                raise ValueError(f"Invalid rank timing: {profile}")
            totals[row["stage"]] = totals.get(row["stage"], 0.0) + elapsed
        profiles.append(totals)
    return profiles


def run_mode(arguments: argparse.Namespace, inputs: Path, mode: str, repetition: int) -> dict[str, Any]:
    destination = arguments.output / f"run-{repetition:02d}-{mode}"
    environment = os.environ | THREAD_ENVIRONMENT
    environment["OMP_DYNAMIC"] = "false"
    # CPATH still uses 128-character paths. Short scratch paths avoid constraining --output.
    directory = Path(tempfile.mkdtemp(prefix="rgmpi-", dir=arguments.scratch_root)).resolve()
    arguments.pending_directories.append((directory, destination))
    if len(os.fsencode(str(directory))) > 100:
        raise ValueError("MPI scratch path is too long; choose a shorter --scratch-root")
    environment["MPI_TMP"] = str(directory)
    command = [str(arguments.mpiexec), *arguments.mpiexec_arg, arguments.numproc_flag,
               str(arguments.ranks), str(arguments.executable),
               "--csf-descriptors" if mode == "descriptors" else "--csf-text",
               str(inputs / ("csfs.bin" if mode == "descriptors" else "rcsf.inp")),
               "--timings", "timings.csv"]
    execution: dict[str, Any] = {
        "mode": mode, "repetition": repetition, "command": command, "directory": str(destination),
        "scratch_directory": str(directory), "thread_environment": THREAD_ENVIRONMENT,
        "status": "running",
    }
    print(f"Running {mode}, repetition {repetition}", flush=True)
    (directory / "stdin.txt").write_text("y\n")
    process: subprocess.Popen[bytes] | None = None
    started = time.perf_counter()
    try:
        with ((directory / "stdin.txt").open("rb") as stdin,
              (directory / "stdout.log").open("wb") as stdout,
              (directory / "stderr.log").open("wb") as stderr):
            process = subprocess.Popen(command, cwd=directory, env=environment, stdin=stdin,
                                       stdout=stdout, stderr=stderr, start_new_session=os.name == "posix")
            execution["returncode"] = process.wait(timeout=arguments.timeout)
        execution["status"] = "exited" if execution["returncode"] == 0 else "failed"
    except BaseException:
        execution["status"] = "interrupted_or_failed"
        if process is not None:
            stop_process(process)
        raise
    finally:
        execution["wall_seconds"] = time.perf_counter() - started
        # Keep all output on scratch until every timed run has completed.
        save_json(directory / "execution.json", execution)
    if execution["returncode"] != 0:
        raise RuntimeError(f"MPI {mode} failed; see {destination / 'stderr.log'}")
    execution["rank_stages"] = read_rank_timings(directory, arguments.ranks)
    root_stages = execution["rank_stages"][0]
    execution["root_input_seconds"] = root_stages["input_header_and_broadcast"] + root_stages["csf_load_and_broadcast"]
    execution["output_sha256"] = {
        path.relative_to(directory).as_posix(): file_identity(path)["sha256"]
        for path in sorted(directory.glob("[0-9][0-9][0-9]/mcp*.*"))
    }
    for rank_index in range(arguments.ranks):
        if not (directory / f"{rank_index:03d}" / f"mcp{rank_index:03d}.30").is_file():
            raise ValueError(f"Missing MCP sparse structure for rank {rank_index}")
    save_json(directory / "execution.json", execution)
    print(f"Finished {mode}: wall={execution['wall_seconds']:.6f}s, "
          f"rank-0 input+broadcast={execution['root_input_seconds']:.6f}s", flush=True)
    return execution


def summarize(runs: list[dict[str, Any]]) -> dict[str, Any]:
    statistics_by_mode = {}
    for mode in ("text", "descriptors"):
        selected = [run for run in runs if run["mode"] == mode]
        statistics_by_mode[mode] = {}
        for quantity in ("wall_seconds", "root_input_seconds"):
            values = [run[quantity] for run in selected]
            statistics_by_mode[mode][quantity] = {
                "values": values, "median": statistics.median(values), "min": min(values), "max": max(values)
            }
    ratios = {}
    for quantity in ("wall_seconds", "root_input_seconds"):
        reference = statistics_by_mode["text"][quantity]["median"]
        candidate = statistics_by_mode["descriptors"][quantity]["median"]
        ratios[quantity] = {"text_over_descriptors": reference / candidate,
                            "median_reduction_fraction": 1 - candidate / reference}
    return {"by_mode": statistics_by_mode, "ratios": ratios}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--csf", type=Path, required=True)
    parser.add_argument("--descriptors", type=Path, required=True)
    parser.add_argument("--verifier", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, default=Path(shutil.which("mpirun") or "mpirun"))
    parser.add_argument("--numproc-flag", default="-np")
    parser.add_argument("--mpiexec-arg", action="append", default=[], help="Repeatable; use --mpiexec-arg=VALUE")
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--scratch-root", type=Path, default=Path(tempfile.gettempdir()))
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    for name in ("executable", "csf", "descriptors", "verifier", "mpiexec", "output", "scratch_root"):
        setattr(arguments, name, getattr(arguments, name).resolve())
    if not 1 <= arguments.ranks <= 999 or arguments.repeat < 1 or not arguments.timeout > 0:
        parser.error("Require 1..999 ranks, a positive repeat count and timeout")
    if arguments.output.exists():
        parser.error("Output directory must not already exist")
    arguments.output.mkdir(parents=True)
    arguments.pending_directories = []
    report: dict[str, Any] = {"status": "running", "ranks": arguments.ranks, "repeat": arguments.repeat, "runs": []}
    try:
        inputs = arguments.output / "inputs"
        inputs.mkdir()
        shutil.copyfile(arguments.csf, inputs / "rcsf.inp")
        shutil.copyfile(arguments.descriptors, inputs / "csfs.bin")
        report["identity"] = {"executable": file_identity(arguments.executable),
                              "verifier": file_identity(arguments.verifier),
                              "text": file_identity(inputs / "rcsf.inp"),
                              "descriptors": file_identity(inputs / "csfs.bin"),
                              "platform": platform.platform(), "scratch_root": str(arguments.scratch_root),
                              "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None}
        verified = verify_descriptors(arguments.verifier, inputs / "csfs.bin", inputs / "rcsf.inp")
        (arguments.output / "array-verification.log").write_text(verified.stdout + verified.stderr)
        report["array_verification"] = "passed"
        expected_outputs = None
        for repetition in range(1, arguments.repeat + 1):
            modes = ("text", "descriptors") if repetition % 2 else ("descriptors", "text")
            for mode in modes:
                execution = run_mode(arguments, inputs, mode, repetition)
                report["runs"].append(execution)
                if expected_outputs is None:
                    expected_outputs = execution["output_sha256"]
                execution["byte_identical"] = execution["output_sha256"] == expected_outputs
                save_json(arguments.output / "summary.json", report)
                if not execution["byte_identical"]:
                    raise ValueError("MCP file set/content differs from the first text run")
        report["statistics"] = summarize(report["runs"])
        report["status"] = "passed"
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        report["status"] = "failed"
        report["error"] = str(error)
        print(error, file=sys.stderr)
    finally:
        if report["status"] == "running":
            report["status"] = "interrupted"
        for scratch_directory, destination in arguments.pending_directories:
            try:
                shutil.move(str(scratch_directory), destination)
            except OSError as error:
                report["status"] = "failed"
                report["error"] = f"Cannot preserve {scratch_directory}: {error}"
        save_json(arguments.output / "summary.json", report)
    if report["status"] != "passed":
        return 1
    print(json.dumps(report["statistics"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
