"""Read-only, bounded Linux /proc sampling for an explicitly selected MPI job.

No MPI calls, signals to the job, subprocesses, or third-party packages are used.
Counters are raw cumulative values; timestamps and units are saved for analysis.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import signal
import sys
import time
from typing import Any

Record = dict[str, Any]
RANK_KEYS = ("OMPI_COMM_WORLD_RANK", "PMIX_RANK", "PMI_RANK", "SLURM_PROCID")
THREAD_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
STATUS_KEYS = {"Uid", "Threads", "VmRSS", "VmHWM", "VmSwap", "Cpus_allowed_list",
               "Mems_allowed_list", "voluntary_ctxt_switches", "nonvoluntary_ctxt_switches"}


def read_text(path: Path, limit: int = 65536) -> Record:
    try:
        with path.open("rb") as stream:
            data = stream.read(limit + 1)
        return {"text": data[:limit].decode(errors="replace"), "truncated": len(data) > limit}
    except OSError as exc:
        return {"error": type(exc).__name__, "errno": exc.errno}


def read_link(path: Path) -> Record:
    try:
        return {"path": os.readlink(path)}
    except OSError as exc:
        return {"error": type(exc).__name__, "errno": exc.errno}


def parse_stat(text: str) -> Record:
    # comm can itself contain spaces and ')' characters.
    left, right = text.index("("), text.rindex(")")
    fields = text[right + 1:].split()
    result: Record = {"pid": int(text[:left]), "comm": text[left + 1:right], "state": fields[0]}
    for name, index in (("ppid", 1), ("minflt", 7), ("majflt", 9),
                        ("utime_ticks", 11), ("stime_ticks", 12), ("threads", 17),
                        ("starttime_ticks", 19), ("vsize_bytes", 20), ("rss_pages", 21),
                        ("processor", 36), ("blkio_delay_ticks", 39)):
        result[name] = int(fields[index]) if len(fields) > index else None
    return result


def stat_at(directory: Path) -> Record | None:
    try:
        return parse_stat(read_text(directory / "stat", 8192)["text"])
    except (KeyError, ValueError, IndexError):
        return None


def identity(stat: Record) -> tuple[int, int]:
    return stat["pid"], stat["starttime_ticks"]


def selected_status(path: Path) -> Record:
    result = read_text(path)
    if "text" in result:
        result["text"] = "\n".join(line for line in result["text"].splitlines()
                                    if line.partition(":")[0] in STATUS_KEYS)
    return result


def selected_environment(path: Path) -> Record:
    result = read_text(path, 1024 * 1024)
    if "text" not in result:
        return result
    pairs = (item.partition("=") for item in result["text"].split("\0"))
    # Never save the full environment (which can contain credentials).
    values = {key: value for key, separator, value in pairs
              if separator and key in RANK_KEYS + THREAD_KEYS}
    ranks = {key: value for key, value in values.items() if key in RANK_KEYS}
    return {"rank_candidates": ranks, "rank": next(iter(ranks.values()), None)
            if len(set(ranks.values())) <= 1 else None,
            "thread_settings": {key: value for key, value in values.items() if key in THREAD_KEYS},
            "truncated": result["truncated"]}


class ProcSampler:
    def __init__(self, pids: list[int], *, proc: Path = Path("/proc"),
                 cgroup_root: Path = Path("/sys/fs/cgroup"), details_interval: float = 30,
                 max_threads: int = 1024, max_fds: int = 128) -> None:
        self.proc, self.cgroup_root = proc, cgroup_root
        self.details_interval, self.max_threads, self.max_fds = details_interval, max_threads, max_fds
        self.roots: dict[int, int] = {}
        self.known: dict[int, int] = {}
        self.last_details: dict[tuple[int, int], float] = {}
        self.last_stack: dict[tuple[int, int], float] = {}
        for pid in pids:
            directory = proc / str(pid)
            stat = stat_at(directory)
            if pid <= 1 or pid == os.getpid() or stat is None or stat["state"] in "ZXx":
                raise ValueError(f"PID {pid} is not an attachable live process")
            if directory.stat().st_uid != os.getuid():
                raise ValueError(f"PID {pid} is not owned by the current user")
            self.roots[pid] = stat["starttime_ticks"]
        self.known.update(self.roots)

    def process_tree(self) -> tuple[dict[int, Record], int]:
        table: dict[int, Record] = {}
        unreadable = 0
        for directory in self.proc.iterdir():
            if not directory.name.isdigit():
                continue
            try:
                if directory.stat().st_uid != os.getuid():
                    continue
            except OSError:
                continue
            stat = stat_at(directory)
            if stat is None:
                unreadable += 1
            elif stat["state"] not in "ZXx":
                table[stat["pid"]] = stat
        selected = {pid: stat for pid, stat in table.items()
                    if self.known.get(pid) == stat["starttime_ticks"]}
        while True:
            added = {pid: stat for pid, stat in table.items()
                     if pid not in selected and pid != os.getpid() and pid not in self.roots
                     and stat["ppid"] in selected
                     and stat["starttime_ticks"] >= selected[stat["ppid"]]["starttime_ticks"]}
            if not added:
                break
            selected.update(added)
        self.known = {pid: stat["starttime_ticks"] for pid, stat in selected.items()}
        return selected, unreadable

    def file_details(self, directory: Path) -> Record:
        result: Record = {key: read_link(directory / key) for key in ("cwd", "exe")}
        result["mountinfo"] = read_text(directory / "mountinfo")
        try:
            fds = sorted((directory / "fd").iterdir(), key=lambda path: int(path.name))
            result["fds"] = {path.name: read_link(path) for path in fds[:self.max_fds]}
            result["fds_truncated"] = len(fds) > self.max_fds
        except (OSError, ValueError) as exc:
            result["fds_error"] = type(exc).__name__
        return result

    def cgroups(self, membership: Record) -> Record:
        # cgroup v2 only, in the sampler's mount/cgroup namespace. Never escape root.
        for line in membership.get("text", "").splitlines():
            if not line.startswith("0::/"):
                continue
            path = self.cgroup_root / line[4:]
            root = self.cgroup_root.resolve()
            if not path.resolve().is_relative_to(root):
                return {"error": "cgroup path outside configured root"}
            records: Record = {}
            while True:
                records[str(path.relative_to(self.cgroup_root))] = {
                    name: read_text(path / name, 16384) for name in (
                        "memory.current", "memory.max", "memory.high", "memory.events",
                        "memory.swap.current", "memory.swap.max", "cpu.max", "cpu.stat",
                        "io.stat", "io.max", "cpu.pressure", "memory.pressure", "io.pressure")}
                if path == self.cgroup_root:
                    break
                path = path.parent
            return {"scope": "sampler namespace; includes ancestors", "values": records}
        return {"error": "cgroup v2 membership unavailable"}

    def sample(self, now: float) -> Record:
        table, unreadable = self.process_tree()
        processes, groups = [], {}
        for pid, stat in sorted(table.items()):
            directory = self.proc / str(pid)
            key = identity(stat)
            item = dict(stat)
            item.update(status=selected_status(directory / "status"),
                        io=read_text(directory / "io"), cgroup=read_text(directory / "cgroup"))
            item["environment"] = selected_environment(directory / "environ")
            membership = item["cgroup"].get("text", "")
            if membership not in groups:
                groups[membership] = self.cgroups(item["cgroup"])
            if now - self.last_details.get(key, -math.inf) >= self.details_interval:
                item["files"] = self.file_details(directory)
                self.last_details[key] = now
            threads = []
            try:
                paths = sorted((directory / "task").iterdir(), key=lambda path: int(path.name))
                item["threads_truncated"] = len(paths) > self.max_threads
            except (OSError, ValueError) as exc:
                paths = []
                item["threads_error"] = type(exc).__name__
            for path in paths[:self.max_threads]:
                task = stat_at(path)
                if task is None:
                    continue
                task["status"] = selected_status(path / "status")
                task["wchan"] = read_text(path / "wchan", 1024)
                task_key = identity(task)
                if task["state"] in ("D", "R") and now - self.last_stack.get(task_key, -math.inf) >= self.details_interval:
                    task["stack"] = read_text(path / "stack", 8192)
                    self.last_stack[task_key] = now
                after = stat_at(path)
                if after is not None and identity(after) == task_key:
                    threads.append(task)
            item["tasks"] = threads
            after = stat_at(directory)
            if after is not None and identity(after) == key:
                processes.append(item)
        active = {identity(stat) for stat in table.values()}
        self.last_details = {key: value for key, value in self.last_details.items() if key in active}
        tasks = {identity(task) for item in processes for task in item["tasks"]}
        self.last_stack = {key: value for key, value in self.last_stack.items() if key in tasks}
        return {"processes": processes, "tracked_live_processes": len(table),
                "unreadable_same_uid_stats": unreadable, "cgroups": groups,
                "system": {name: read_text(self.proc / name) for name in (
                    "stat", "meminfo", "vmstat", "diskstats", "loadavg",
                    "pressure/cpu", "pressure/memory", "pressure/io")},
                "phase": "unknown; external sampling does not identify Fortran call sites"}


class BoundedJsonl:
    def __init__(self, path: Path, max_bytes: int) -> None:
        self.stream = path.open("xb")
        self.max_bytes, self.bytes_written = max_bytes, 0

    def append(self, record: Record) -> bool:
        data = (json.dumps(record, ensure_ascii=True, separators=(",", ":")) + "\n").encode()
        if self.bytes_written + len(data) > self.max_bytes:
            return False
        self.stream.write(data)
        self.stream.flush()
        self.bytes_written += len(data)
        return True

    def close(self) -> None:
        self.stream.close()


def positive_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--pid", type=positive_int, action="append", required=True,
                        help="local launcher or rank PID; repeat for separate process trees")
    result.add_argument("--output", type=Path, required=True, help="new local output directory")
    result.add_argument("--interval", type=positive_float, default=2.0)
    result.add_argument("--duration", type=positive_float, default=86400.0, help="maximum seconds")
    result.add_argument("--details-interval", type=positive_float, default=30.0)
    result.add_argument("--max-mib", type=positive_float, default=256.0, help="samples.jsonl size cap")
    result.add_argument("--max-threads", type=positive_int, default=1024, help="per process sample")
    result.add_argument("--max-fds", type=positive_int, default=128, help="per process detail sample")
    return result


def run(args: argparse.Namespace) -> int:
    if not sys.platform.startswith("linux"):
        raise ValueError("Linux /proc is required; no job was started or modified")
    sampler = ProcSampler(args.pid, details_interval=args.details_interval,
                          max_threads=args.max_threads, max_fds=args.max_fds)
    args.output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    metadata = {"schema_version": 1, "utc_start": datetime.now(timezone.utc).isoformat(),
                "hostname": os.uname().nodename, "kernel": os.uname().release,
                "boot_id": read_text(Path("/proc/sys/kernel/random/boot_id")),
                "root_identities": sampler.roots, "clock_ticks_per_second": os.sysconf("SC_CLK_TCK"),
                "page_size_bytes": os.sysconf("SC_PAGE_SIZE"),
                "options": {key: str(value) if isinstance(value, Path) else value
                            for key, value in vars(args).items()},
                "scope": "local same-UID descendants; reparented children retained only if already seen",
                "limitations": "non-atomic snapshots; unknown phase; permissions and namespaces may restrict coverage"}
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    writer = BoundedJsonl(args.output / "samples.jsonl", int(args.max_mib * 1024**2))
    interrupted = False

    def stop(signum: int, frame: Any) -> None:
        nonlocal interrupted
        interrupted = True

    old_handlers = {sig: signal.signal(sig, stop) for sig in (signal.SIGINT, signal.SIGTERM)}
    count, states, d_samples, max_processes, cost = 0, Counter(), 0, 0, 0.0
    reason, error = "duration", None
    print(f"Sampling {args.pid} into {args.output}; only the sampler stops at its limits.", flush=True)
    try:
        while time.monotonic() - start < args.duration and not interrupted:
            tick = time.monotonic()
            utc = datetime.now(timezone.utc).isoformat()
            record = sampler.sample(tick)
            elapsed = time.monotonic() - tick
            record.update(utc=utc, elapsed_seconds=tick - start,
                          collection_seconds=elapsed)
            if not writer.append(record):
                reason = "sample_size_limit"
                break
            count += 1
            cost += elapsed
            max_processes = max(max_processes, len(record["processes"]))
            sample_states = Counter(task["state"] for proc in record["processes"] for task in proc["tasks"])
            states.update(sample_states)
            d_samples += int(sample_states["D"] > 0)
            if record["tracked_live_processes"] == 0:
                reason = "targets_exited_or_inaccessible"
                break
            remaining = min(args.interval - elapsed, args.duration - (time.monotonic() - start))
            deadline = time.monotonic() + max(0.0, remaining)
            while not interrupted and time.monotonic() < deadline:
                time.sleep(min(0.2, max(0.0, deadline - time.monotonic())))
        if interrupted:
            reason = "sampler_interrupted"
    except Exception as exc:
        reason, error = "sampler_error", f"{type(exc).__name__}: {exc}"
    finally:
        writer.close()
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        summary = {"stop_reason": reason, "error": error, "samples": count,
                   "samples_with_D": d_samples, "task_state_observations": dict(states),
                   "max_observed_processes": max_processes, "sample_bytes": writer.bytes_written,
                   "collection_seconds": cost, "elapsed_seconds": time.monotonic() - start,
                   "interpretation": "D sample counts are not durations or a root-cause diagnosis; inspect missing data and rank coverage"}
        (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if error or count == 0 else 0


def main() -> int:
    args = parser().parse_args()
    try:
        return run(args)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
