#!/usr/bin/env python3
"""Compare native rangular inputs, with explicit offline support for historical binaries."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import signal
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import mcp_compare
from verify_csf_descriptors import SUCCESS_MARKER

THREAD_ENVIRONMENT = {
    name: "1"
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    )
}


def file_identity(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def save_json(path: Path, report: dict[str, Any]) -> None:
    temporary = path.with_suffix(".json.pending")
    temporary.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def stop_process(process: subprocess.Popen[bytes]) -> None:
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    elif process.poll() is None:
        process.kill()
    process.wait()


def wait_measured(
    process: subprocess.Popen[bytes], timeout_seconds: float | None
) -> tuple[int, int | None]:
    """Collect this child's OS high-water RSS, never a cumulative children maximum."""
    if not hasattr(os, "wait4") or sys.platform not in {"darwin", "linux"}:
        return process.wait(timeout=timeout_seconds), None
    deadline = None if timeout_seconds is None else time.monotonic() + timeout_seconds
    while True:
        pid, status, usage = os.wait4(process.pid, os.WNOHANG)
        if pid:
            process.returncode = os.waitstatus_to_exitcode(status)
            unit = 1 if sys.platform == "darwin" else 1024
            return process.returncode, int(usage.ru_maxrss * unit)
        if deadline is not None and time.monotonic() >= deadline:
            raise subprocess.TimeoutExpired(process.args, timeout_seconds)
        time.sleep(0.01)


def run_logged(
    command: list[str],
    directory: Path,
    stdin_path: Path | None,
    timeout_seconds: float | None,
    *,
    input_csf: Path | None = None,
    input_descriptors: Path | None = None,
) -> dict[str, Any]:
    directory.mkdir()
    if input_csf is not None:
        shutil.copyfile(input_csf, directory / "rcsf.inp")
    if input_descriptors is not None:
        shutil.copyfile(input_descriptors, directory / "rcsf_descriptors.parquet")
        shutil.copyfile(input_descriptors.with_name("csfs_header.toml"), directory / "rcsf_header.toml")
    thread_environment = THREAD_ENVIRONMENT.copy()
    execution: dict[str, Any] = {
        "command": command,
        "directory": str(directory),
        "status": "running",
        "thread_environment": thread_environment,
    }
    save_json(directory / "execution.json", execution)
    print(f"Running {directory.name}", flush=True)
    start_time = time.perf_counter()
    process: subprocess.Popen[bytes] | None = None
    try:
        with (
            (directory / "stdout.log").open("wb") as stdout,
            (directory / "stderr.log").open("wb") as stderr,
            (stdin_path or Path(os.devnull)).open("rb") as stdin,
        ):
            process = subprocess.Popen(
                command,
                cwd=directory,
                stdin=stdin,
                stdout=stdout,
                stderr=stderr,
                env=os.environ | thread_environment,
                start_new_session=os.name == "posix",
            )
            execution["returncode"], execution["peak_rss_bytes"] = wait_measured(
                process, timeout_seconds
            )
            execution["peak_rss_method"] = (
                "wait4.ru_maxrss" if execution["peak_rss_bytes"] is not None else "unavailable"
            )
        execution["status"] = "exited" if execution["returncode"] == 0 else "failed"
    except subprocess.TimeoutExpired:
        execution["status"] = "timeout"
        if process is not None:
            stop_process(process)
        raise
    except BaseException:
        execution["status"] = "interrupted_or_failed"
        if process is not None:
            stop_process(process)
        raise
    finally:
        execution["wall_seconds"] = time.perf_counter() - start_time
        save_json(directory / "execution.json", execution)
    if execution["returncode"] != 0:
        raise RuntimeError(
            f"{directory.name} exited {execution['returncode']}; see its stdout.log/stderr.log"
        )
    return execution


def reuse_statistics(path: Path) -> dict[str, Any]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        expected = [
            "block", "pairs", "hits", "misses", "bypasses", "not_stored", "entries",
            "probes", "cache_reserved_bytes", "cache_used_bytes", "cache_limit_bytes",
        ]
        if reader.fieldnames != expected:
            raise ValueError(f"invalid reuse statistics header: {path}")
        blocks = [{key: int(value) for key, value in row.items()} for row in reader]
    if not blocks:
        raise ValueError(f"empty reuse statistics: {path}")
    for block in blocks:
        if (
            any(value < 0 for value in block.values())
            or block["pairs"] != block["hits"] + block["misses"] + block["bypasses"]
            or block["entries"] + block["not_stored"] != block["misses"]
            or not block["cache_used_bytes"] <= block["cache_reserved_bytes"] <= block["cache_limit_bytes"]
        ):
            raise ValueError(f"inconsistent reuse statistics: {path}")
    totals = {key: sum(block[key] for block in blocks) for key in expected[1:8]}
    return {
        "blocks": blocks,
        "totals": totals,
        "hit_rate": totals["hits"] / totals["pairs"] if totals["pairs"] else 0.0,
        "peak_cache_reserved_bytes": max(block["cache_reserved_bytes"] for block in blocks),
        "peak_cache_used_bytes": max(block["cache_used_bytes"] for block in blocks),
        "memory_scope": "Cache array storage only; process peak RSS includes all allocations and runtime overhead.",
    }


def stage_timings(path: Path) -> dict[str, float]:
    with path.open() as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != [
            "stage",
            "block",
            "channel",
            "wall_seconds",
            "records",
            "storage",
        ]:
            raise ValueError(f"invalid timings header: {path}")
        rows = list(reader)
    if not rows or rows[-1]["stage"] != "program_total_including_interaction_input":
        raise ValueError(f"incomplete candidate timings: {path}")
    totals: dict[str, float] = {}
    for row in rows:
        elapsed = float(row["wall_seconds"])
        if not math.isfinite(elapsed) or elapsed < 0:
            raise ValueError(f"invalid timing value: {path}")
        totals[row["stage"]] = totals.get(row["stage"], 0.0) + elapsed
    return totals


def compare_outputs(
    reference: Path, candidate: Path, arguments: argparse.Namespace
) -> dict[str, Any]:
    identity = mcp_compare.input_identity(reference, candidate, explicit=False)
    reference_data = mcp_compare.Dataset(reference)
    candidate_data = mcp_compare.Dataset(candidate)
    if max(reference_data.header.nw, candidate_data.header.nw) >= arguments.keyorb:
        raise ValueError("orbital count must be smaller than --keyorb")
    comparison = mcp_compare.compare(
        reference_data,
        candidate_data,
        atol=arguments.atol,
        rtol=arguments.rtol,
        keyorb=arguments.keyorb,
    )
    comparison["csf_identity"] = identity
    reference_files = {
        path.name: file_identity(path)["sha256"] for path in reference_data.paths
    }
    candidate_files = {
        path.name: file_identity(path)["sha256"] for path in candidate_data.paths
    }
    comparison["byte_identical"] = reference_files == candidate_files
    comparison["reference_mcp_sha256"] = reference_files
    comparison["candidate_mcp_sha256"] = candidate_files
    return comparison


def execute_comparison(arguments: argparse.Namespace) -> int:
    output = arguments.output.resolve()
    if output.exists():
        raise ValueError(
            "--output must be a new directory; previous results are never overwritten"
        )
    paths = {
        "reference": arguments.reference.resolve(),
        "candidate": arguments.candidate.resolve(),
        "csf": arguments.csf.resolve(),
    }
    for option in ("stdin_file", "binary", "descriptors", "header"):
        value = getattr(arguments, option)
        if value is not None:
            paths[option] = value.resolve()
    if arguments.binary is not None:
        raise ValueError("GCSDESC1 is retired; supply --descriptors and --header")
    descriptor_mode = arguments.descriptors is not None
    if descriptor_mode:
        paths["verifier"] = (
            arguments.verifier or paths["candidate"].with_name("csf_descriptor_verify")
        ).resolve()
    for name, path in paths.items():
        if not path.is_file():
            raise ValueError(f"{name} is not a file: {path}")
        if name in {"reference", "candidate", "verifier"} and not os.access(
            path, os.X_OK
        ):
            raise ValueError(f"{name} is not executable: {path}")
    identities = {name: file_identity(path) for name, path in paths.items()}
    output.mkdir(parents=True)
    snapshots = output / "inputs"
    snapshots.mkdir()
    summary: dict[str, Any] = {
        "schema_version": 1,
        "status": "running",
        "inputs": identities,
        "same_executable_bytes": identities["reference"]["sha256"]
        == identities["candidate"]["sha256"],
        "platform": {"system": platform.system(), "machine": platform.machine()},
        "thread_environment": THREAD_ENVIRONMENT,
        "repetitions": arguments.repeat,
        "candidate_reuse_cache_mib": arguments.candidate_reuse_cache_mib,
        "reference_timings": arguments.reference_timings,
        "atol": arguments.atol,
        "rtol": arguments.rtol,
        "keyorb": arguments.keyorb,
        "runs": [],
        "comparisons": [],
        "timing_scope": "External process wall time only; input copies, verification, preparation and MCP comparison excluded.",
        "cache_policy": "No cache flush. Input copying/validation can warm the page cache; these are not cold-cache benchmarks.",
    }
    summary_path = output / "summary.json"
    save_json(summary_path, summary)
    try:
        for name in ("csf", "stdin_file", "binary", "descriptors", "header"):
            if name not in paths:
                continue
            snapshot = (
                snapshots
                / {
                    "csf": "rcsf.inp",
                    "stdin_file": "rangular.stdin",
                    "binary": "csfs.bin",
                    "descriptors": "csfs_descriptors.parquet",
                    "header": "csfs_header.toml",
                }[name]
            )
            shutil.copyfile(paths[name], snapshot)
            if file_identity(snapshot)["sha256"] != identities[name]["sha256"]:
                raise ValueError(f"{name} changed while taking its input snapshot")
        stdin_path = snapshots / "rangular.stdin"
        if "stdin_file" not in paths:
            stdin_path.write_bytes(b"y\n")
        summary["effective_stdin"] = file_identity(stdin_path)
        if descriptor_mode:
            summary["effective_binary"] = file_identity(snapshots / "csfs_descriptors.parquet")
            verification_directory = output / "array_verification"
            summary["array_verification"] = run_logged(
                [
                    str(paths["verifier"]),
                    str(snapshots / "csfs_descriptors.parquet"),
                    str(snapshots / "rcsf.inp"),
                ],
                verification_directory,
                None,
                arguments.timeout,
            )
            with (verification_directory / "stdout.log").open() as stream:
                if not any(line.strip() == SUCCESS_MARKER for line in stream):
                    raise ValueError(
                        "array verifier exited without its success marker (possible legacy STOP)"
                    )
            summary["array_verification"]["status"] = "verified"
        modes = (["reference_text" if arguments.reference_text_oracle else "reference_default", "candidate_default"]
                 if descriptor_mode else ["reference_text", "candidate_text"])
        if descriptor_mode:
            modes.append("candidate_descriptors")
        mismatch_found = False
        for trial_index in range(arguments.repeat):
            rotation = trial_index % len(modes)
            execution_order = modes[rotation:] + modes[:rotation]
            trial_directories: dict[str, Path] = {}
            for mode in execution_order:
                directory = output / f"trial_{trial_index + 1:03d}_{mode}"
                binary_name = "reference" if mode in ("reference_text", "reference_default") else "candidate"
                executable = paths[binary_name]
                if (
                    file_identity(executable)["sha256"]
                    != identities[binary_name]["sha256"]
                ):
                    raise ValueError(
                        f"{binary_name} executable changed during the comparison"
                    )
                command = [str(executable)]
                if mode == "candidate_descriptors":
                    command += ["--csf-descriptors", str(snapshots / "csfs_descriptors.parquet")]
                if mode not in ("reference_text", "reference_default") or arguments.reference_timings:
                    command += ["--timings", "timings.csv"]
                if mode not in ("reference_text", "reference_default"):
                    if arguments.candidate_reuse_cache_mib is not None:
                        command += [
                            "--reuse-cache-mib", str(arguments.candidate_reuse_cache_mib),
                            "--reuse-stats", "reuse.csv",
                        ]
                execution = run_logged(
                    command,
                    directory,
                    stdin_path,
                    arguments.timeout,
                    input_csf=None
                    if descriptor_mode and mode != "reference_text"
                    else snapshots / "rcsf.inp",
                    input_descriptors=(snapshots / "csfs_descriptors.parquet") if descriptor_mode and mode in ("reference_default", "candidate_default") else None,
                )
                if (
                    file_identity(executable)["sha256"]
                    != identities[binary_name]["sha256"]
                ):
                    raise ValueError(
                        f"{binary_name} executable changed during execution"
                    )
                if descriptor_mode and mode != "reference_text":
                    if (directory / "rcsf.inp").exists():
                        raise ValueError("descriptor run unexpectedly created rcsf.inp")
                    shutil.copyfile(snapshots / "rcsf.inp", directory / "rcsf.inp")
                    execution["text_provenance_added_after_execution"] = True
                if (
                    file_identity(directory / "rcsf.inp")["sha256"]
                    != identities["csf"]["sha256"]
                ):
                    raise ValueError("a run modified its CSF input")
                dataset = mcp_compare.Dataset(directory)
                if dataset.header.nprocs != 1 or any(
                    path.name.startswith("mcp0") for path in dataset.paths
                ):
                    raise ValueError(
                        "this runner requires ordinary unsharded mcp.* output (no MPI/CSFG)"
                    )
                if dataset.header.nw >= arguments.keyorb:
                    raise ValueError("orbital count must be smaller than --keyorb")
                execution["mcp_inspection"] = mcp_compare.inspect(
                    dataset, arguments.keyorb
                )
                if mode not in ("reference_text", "reference_default") or arguments.reference_timings:
                    execution["stage_wall_seconds"] = stage_timings(
                        directory / "timings.csv"
                    )
                if mode not in ("reference_text", "reference_default") and arguments.candidate_reuse_cache_mib is not None:
                    execution["angular_reuse"] = reuse_statistics(directory / "reuse.csv")
                execution.update(mode=mode, trial=trial_index + 1, status="validated")
                summary["runs"].append(execution)
                save_json(directory / "execution.json", execution)
                save_json(summary_path, summary)
                trial_directories[mode] = directory
            for mode in modes[1:]:
                comparison = compare_outputs(
                    trial_directories[modes[0]],
                    trial_directories[mode],
                    arguments,
                )
                report_path = output / f"comparison_{trial_index + 1:03d}_{mode}.json"
                save_json(report_path, comparison)
                mismatch_found |= comparison["status"] != "equal"
                summary["comparisons"].append(
                    {
                        "trial": trial_index + 1,
                        "mode": mode,
                        "status": comparison["status"],
                        "byte_identical": comparison["byte_identical"],
                        "report": str(report_path),
                    }
                )
            if trial_index:
                comparison = compare_outputs(
                    output / f"trial_001_{modes[0]}",
                    trial_directories[modes[0]],
                    arguments,
                )
                report_path = (
                    output / f"reference_repeatability_{trial_index + 1:03d}.json"
                )
                save_json(report_path, comparison)
                mismatch_found |= comparison["status"] != "equal"
                summary["comparisons"].append(
                    {
                        "trial": trial_index + 1,
                        "mode": "reference_repeatability",
                        "status": comparison["status"],
                        "report": str(report_path),
                    }
                )
            save_json(summary_path, summary)
        medians = {
            mode: statistics.median(
                run["wall_seconds"] for run in summary["runs"] if run["mode"] == mode
            )
            for mode in modes
        }
        summary["median_wall_seconds"] = medians
        summary["median_peak_rss_bytes"] = {
            mode: statistics.median(
                run["peak_rss_bytes"] for run in summary["runs"]
                if run["mode"] == mode and run["peak_rss_bytes"] is not None
            )
            for mode in modes
            if any(run["mode"] == mode and run["peak_rss_bytes"] is not None for run in summary["runs"])
        }
        summary["median_stage_wall_seconds"] = {
            mode: {
                stage: statistics.median(
                    run["stage_wall_seconds"][stage] for run in summary["runs"] if run["mode"] == mode
                )
                for stage in next(run for run in summary["runs"] if run["mode"] == mode).get("stage_wall_seconds", {})
            }
            for mode in modes
        }
        if not mismatch_found:
            summary["reference_over_candidate_wall_ratio"] = {
                mode: medians[modes[0]] / medians[mode]
                for mode in modes[1:]
                if medians[mode] > 0
            }
        summary["status"] = "different" if mismatch_found else "equal"
        return 1 if mismatch_found else 0
    except KeyboardInterrupt:
        summary["status"] = "interrupted"
        return 130
    except Exception as error:  # noqa: BLE001 -- CLI boundary preserves a failed report for any producer/parser error.
        # This is the CLI orchestration boundary: incomplete evidence must never
        # retain a misleading running/equal status, including format-library errors.
        summary.update(
            status="error", error=str(error), error_type=type(error).__name__
        )
        print(f"Comparison failed: {error}", file=sys.stderr)
        return 2
    finally:
        save_json(summary_path, summary)
        print(f"{summary['status']}: {summary_path}", flush=True)


def positive_integer(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def positive_seconds(value: str) -> float:
    number = mcp_compare.nonnegative(value)
    if number == 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="original serial rangular executable",
    )
    parser.add_argument(
        "--candidate", type=Path, required=True, help="new suite rangular executable"
    )
    parser.add_argument(
        "--csf",
        type=Path,
        required=True,
        help="original CSF text (copied; never modified)",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="new directory for all evidence"
    )
    parser.add_argument(
        "--stdin-file",
        type=Path,
        help="interactive answers; default is full interaction: y",
    )
    parser.add_argument("--reference-text-oracle", action="store_true",
                        help="Use original text only for an explicitly selected historical reference executable")
    parser.add_argument("--legacy-text-oracle", action="store_true",
                        help="Offline comparison of two historical text-only executables")
    descriptors = parser.add_mutually_exclusive_group()
    descriptors.add_argument("--binary", type=Path, help=argparse.SUPPRESS)
    descriptors.add_argument(
        "--descriptors", type=Path, help="native V3 descriptor Parquet; no conversion"
    )
    parser.add_argument(
        "--header", type=Path, help="matching header, required with --descriptors"
    )
    parser.add_argument(
        "--verifier",
        type=Path,
        help="default: csf_descriptor_verify beside --candidate",
    )
    parser.add_argument("--repeat", type=positive_integer, default=1)
    parser.add_argument(
        "--candidate-reuse-cache-mib", type=positive_integer,
        help="enable conservative serial angular reuse and collect reuse.csv for each candidate run",
    )
    parser.add_argument(
        "--reference-timings", action="store_true",
        help="also pass --timings to a suite reference (older upstream binaries may not support it)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_seconds,
        help="seconds per child process; default: unlimited",
    )
    parser.add_argument("--atol", type=mcp_compare.nonnegative, default=0.0)
    parser.add_argument("--rtol", type=mcp_compare.nonnegative, default=0.0)
    parser.add_argument(
        "--keyorb",
        type=int,
        default=215,
        help="common integral-label radix of both builds",
    )
    arguments = parser.parse_args()
    if not (arguments.descriptors or arguments.binary or arguments.legacy_text_oracle):
        parser.error("Current programs require --descriptors and --header")
    if arguments.reference_text_oracle and not arguments.descriptors:
        parser.error("--reference-text-oracle requires --descriptors and --header")
    if arguments.legacy_text_oracle and (arguments.descriptors or arguments.binary):
        parser.error("--legacy-text-oracle cannot be combined with descriptor inputs")
    if bool(arguments.descriptors) != bool(arguments.header):
        parser.error("--descriptors and --header must be provided together")
    if arguments.verifier and not (arguments.binary or arguments.descriptors):
        parser.error("--verifier requires --binary or --descriptors")
    if arguments.keyorb < 2:
        parser.error("--keyorb must be >= 2")
    if arguments.candidate_reuse_cache_mib is not None:
        if arguments.candidate_reuse_cache_mib > 4096:
            parser.error("--candidate-reuse-cache-mib must be <= 4096")
    try:
        return execute_comparison(arguments)
    except (OSError, ValueError) as error:
        print(f"Cannot start comparison: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
