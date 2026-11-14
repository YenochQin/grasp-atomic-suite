"""Disposable /proc fixtures: sampling scope, races, missing data and output limits."""

from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sample_rmcdhf_mpi as sample


def stat_line(pid: int, parent: int, start: int, state: str = "S") -> str:
    fields = ["0"] * 50
    for index, value in {0: state, 1: parent, 7: 11, 9: 3, 11: 100, 12: 20,
                         17: 2, 19: start, 20: 8192, 21: 2, 36: 4, 39: 7}.items():
        fields[index] = str(value)
    return f"{pid} (worker ) thread) {' '.join(fields)}\n"


class SamplingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.proc, self.cg = self.base / "proc", self.base / "cgroup"
        self.proc.mkdir()
        (self.cg / "job").mkdir(parents=True)
        (self.cg / "memory.max").write_text("2000\n")
        (self.cg / "job/memory.max").write_text("max\n")

    def process(self, pid: int, parent: int, start: int, state: str = "S") -> Path:
        path = self.proc / str(pid)
        path.mkdir(exist_ok=True)
        (path / "stat").write_text(stat_line(pid, parent, start, state))
        (path / "status").write_text("Uid:\t1000\nVmRSS:\t8 kB\nCpus_allowed_list:\t0-7\n")
        (path / "environ").write_bytes(b"OMPI_COMM_WORLD_RANK=0\0SECRET=never-save\0OMP_NUM_THREADS=1\0")
        (path / "cgroup").write_text("0::/job\n")
        (path / "io").write_text("read_bytes: 4096\nwrite_bytes: 1024\n")
        (path / "mountinfo").write_text("fixture mount information\n")
        (path / "fd").mkdir(exist_ok=True)
        task = path / "task" / str(pid)
        task.mkdir(parents=True, exist_ok=True)
        (task / "stat").write_text(stat_line(pid, parent, start, state))
        (task / "wchan").write_text("folio_wait_bit_common" if state == "D" else "0")
        (task / "stack").write_text("fixture kernel frame\n")
        return path

    def sampler(self, **kwargs: object) -> sample.ProcSampler:
        return sample.ProcSampler([100], proc=self.proc, cgroup_root=self.cg, **kwargs)

    def test_stat_fields_and_parentheses(self) -> None:
        stat = sample.parse_stat(stat_line(100, 50, 900, "D"))
        self.assertEqual(stat["comm"], "worker ) thread")
        self.assertEqual((stat["ppid"], stat["starttime_ticks"], stat["rss_pages"]), (50, 900, 2))
        self.assertEqual((stat["utime_ticks"], stat["processor"], stat["blkio_delay_ticks"]), (100, 4, 7))
        self.assertEqual((stat["minflt"], stat["majflt"]), (11, 3))

    def test_tree_isolation_reparenting_and_pid_reuse(self) -> None:
        self.process(100, 1, 1000)
        child = self.process(101, 100, 1001)
        self.process(102, 101, 1002)
        self.process(200, 1, 500)
        sampler = self.sampler()
        self.assertEqual(set(sampler.process_tree()[0]), {100, 101, 102})
        (child / "stat").write_text(stat_line(101, 1, 1001))
        self.assertEqual(set(sampler.process_tree()[0]), {100, 101, 102})
        (self.proc / "100/stat").write_text(stat_line(100, 1, 2000))
        (child / "stat").write_text(stat_line(101, 200, 2001))
        (self.proc / "102/stat").unlink()
        self.assertEqual(sampler.process_tree()[0], {})

    def test_dead_and_foreign_roots_rejected(self) -> None:
        self.process(100, 1, 1000, "Z")
        with self.assertRaises(ValueError):
            self.sampler()
        self.process(100, 1, 1000)
        with patch.object(sample.os, "getuid", return_value=os.getuid() + 1):
            with self.assertRaises(ValueError):
                self.sampler()

    def test_environment_whitelist_and_conflicting_ranks(self) -> None:
        path = self.process(100, 1, 1000)
        result = sample.selected_environment(path / "environ")
        self.assertEqual(result["rank"], "0")
        self.assertNotIn("never-save", json.dumps(result))
        with (path / "environ").open("ab") as stream:
            stream.write(b"PMI_RANK=2\0")
        self.assertIsNone(sample.selected_environment(path / "environ")["rank"])

    def test_wait_and_running_samples_details_throttled(self) -> None:
        self.process(100, 1, 1000, "D")
        self.process(101, 100, 1001, "R")
        sampler = self.sampler(details_interval=30)
        first = sampler.sample(0)
        tasks = [proc["tasks"][0] for proc in first["processes"]]
        self.assertEqual([task["state"] for task in tasks], ["D", "R"])
        self.assertTrue(all("stack" in task for task in tasks))
        self.assertIn("unknown", first["phase"])
        self.assertIn("error", first["system"]["diskstats"])
        second = sampler.sample(2)
        self.assertTrue(all("files" not in proc for proc in second["processes"]))
        self.assertTrue(all("stack" not in proc["tasks"][0] for proc in second["processes"]))
        self.assertIn("files", sampler.sample(30)["processes"][0])

    def test_cgroup_ancestor_limits_and_escape_rejection(self) -> None:
        self.process(100, 1, 1000)
        sampler = self.sampler()
        result = sampler.cgroups({"text": "0::/job\n"})["values"]
        self.assertEqual(result["job"]["memory.max"]["text"], "max\n")
        self.assertEqual(result["."]["memory.max"]["text"], "2000\n")
        self.assertIn("error", sampler.cgroups({"text": "0::/../outside\n"}))
        self.assertIn("error", sampler.cgroups({"text": "2:memory:/job\n"}))

    def test_process_replaced_during_read_is_discarded(self) -> None:
        path = self.process(100, 1, 1000)
        sampler = self.sampler()
        original = sampler.file_details

        def replacing(directory: Path) -> sample.Record:
            (path / "stat").write_text(stat_line(100, 1, 2000))
            return original(directory)

        with patch.object(sampler, "file_details", side_effect=replacing):
            self.assertEqual(sampler.sample(0)["processes"], [])

    def test_thread_and_fd_limits_are_explicit(self) -> None:
        path = self.process(100, 1, 1000)
        other = path / "task/101"
        other.mkdir()
        (other / "stat").write_text(stat_line(101, 1, 1000, "D"))
        for fd in (3, 4):
            (path / "fd" / str(fd)).symlink_to(self.base / f"mcp.{fd}")
        item = self.sampler(max_threads=1, max_fds=1).sample(0)["processes"][0]
        self.assertTrue(item["threads_truncated"])
        self.assertTrue(item["files"]["fds_truncated"])
        self.assertEqual(len(item["tasks"]), 1)

    def test_missing_and_denied_data_not_zero(self) -> None:
        self.assertIn("error", sample.read_text(self.base / "missing"))
        with patch.object(Path, "open", side_effect=PermissionError(13, "denied")):
            self.assertEqual(sample.read_text(self.base / "stack")["errno"], 13)

    def test_size_cap_preserves_complete_json_records(self) -> None:
        path = self.base / "samples.jsonl"
        writer = sample.BoundedJsonl(path, 20)
        self.assertTrue(writer.append({"a": 1}))
        self.assertFalse(writer.append({"payload": "too large"}))
        writer.close()
        self.assertEqual(json.loads(path.read_text()), {"a": 1})
        self.assertLessEqual(path.stat().st_size, 20)
        with self.assertRaises(FileExistsError):
            sample.BoundedJsonl(path, 20)

    def test_run_lifecycle_and_existing_output_protection(self) -> None:
        self.process(100, 1, 1000, "D")
        sampler = self.sampler()
        args = sample.parser().parse_args(["--pid", "100", "--output", str(self.base / "output"),
                                          "--duration", "0.03", "--interval", "0.01"])
        with patch.object(sample.sys, "platform", "linux"), patch.object(sample, "ProcSampler", return_value=sampler):
            with redirect_stdout(io.StringIO()):
                self.assertEqual(sample.run(args), 0)
            with self.assertRaises(FileExistsError):
                sample.run(args)
        summary = json.loads((args.output / "summary.json").read_text())
        self.assertGreater(summary["samples_with_D"], 0)
        self.assertEqual(summary["stop_reason"], "duration")
        self.assertNotIn("never-save", (args.output / "samples.jsonl").read_text())
        self.assertTrue((self.proc / "100/stat").exists())

    def test_cli_rejects_nonfinite_intervals(self) -> None:
        for value in ("nan", "inf", "0", "-1"):
            with self.assertRaises(sample.argparse.ArgumentTypeError):
                sample.positive_float(value)

    def test_no_samples_due_to_cap_returns_failure(self) -> None:
        self.process(100, 1, 1000)
        args = sample.parser().parse_args(["--pid", "100", "--output", str(self.base / "small"),
                                          "--max-mib", "0.000001"])
        with patch.object(sample.sys, "platform", "linux"), patch.object(sample, "ProcSampler", return_value=self.sampler()):
            with redirect_stdout(io.StringIO()):
                self.assertEqual(sample.run(args), 1)
        summary = json.loads((args.output / "summary.json").read_text())
        self.assertEqual(summary["stop_reason"], "sample_size_limit")
        self.assertEqual(summary["samples"], 0)

    def test_exit_and_sampler_exception_write_summary(self) -> None:
        self.process(100, 1, 1000)
        for fail in (False, True):
            sampler = self.sampler()
            args = sample.parser().parse_args(["--pid", "100", "--output", str(self.base / str(fail))])
            empty = {"processes": [], "tracked_live_processes": 0}
            with patch.object(sample.sys, "platform", "linux"), patch.object(sample, "ProcSampler", return_value=sampler):
                with patch.object(sampler, "sample", side_effect=OSError("fixture failure") if fail else None, return_value=empty):
                    with redirect_stdout(io.StringIO()):
                        self.assertEqual(sample.run(args), int(fail))
            summary = json.loads((args.output / "summary.json").read_text())
            self.assertEqual(summary["stop_reason"], "sampler_error" if fail else "targets_exited_or_inaccessible")

    def test_sampler_signal_does_not_signal_job(self) -> None:
        self.process(100, 1, 1000)
        sampler = self.sampler()
        original = sampler.sample

        def stop_sampler(now: float) -> sample.Record:
            sample.signal.getsignal(sample.signal.SIGTERM)(sample.signal.SIGTERM, None)
            return original(now)

        args = sample.parser().parse_args(["--pid", "100", "--output", str(self.base / "signal")])
        with patch.object(sample.sys, "platform", "linux"), patch.object(sample, "ProcSampler", return_value=sampler):
            with patch.object(sampler, "sample", side_effect=stop_sampler), patch.object(os, "kill") as kill:
                with redirect_stdout(io.StringIO()):
                    self.assertEqual(sample.run(args), 0)
                kill.assert_not_called()
        summary = json.loads((args.output / "summary.json").read_text())
        self.assertEqual(summary["stop_reason"], "sampler_interrupted")


if __name__ == "__main__":
    unittest.main()
