"""Exercise the production mixing stream writer and its completion marker."""
from __future__ import annotations

import resource
import select
import signal
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


def ignore_file_size_signal() -> None:
    signal.signal(signal.SIGXFSZ, signal.SIG_IGN)


def limit_file_size() -> None:
    ignore_file_size_signal()
    resource.setrlimit(resource.RLIMIT_FSIZE, (256, 256))


def exercise(command: list[str], scenario: str, expected: str = "") -> None:
    with tempfile.TemporaryDirectory(prefix="mix-output-") as directory:
        result = subprocess.run(
            [*command, scenario],
            cwd=directory,
            capture_output=True,
            text=True,
            timeout=20,
            preexec_fn=limit_file_size if scenario == "limited" else None,
        )
        output = result.stdout + result.stderr
        payload = (Path(directory) / "fixture.mix").read_bytes()
        if expected:
            assert result.returncode != 0 and expected in output, output
            assert payload[:8] == b"GRMIXTMP", payload[:8]
            return
        assert result.returncode == 0, output
        assert len(payload) == 272 and payload[:8] == b"CSF--MIX"
        assert struct.unpack_from("<5q", payload, 8) == (2, 2, 2, 4, 2)
        assert payload[48:112] == b"a" * 64
        assert struct.unpack_from("<12q", payload, 112) == (0, 1, 2, 1, 0, 240, 2, -1, 2, 1, 1, 256)
        assert struct.unpack_from("<qdqd4d", payload, 208) == (2, -3.0, 1, -1.5, 1.0, 0.0, 0.0, 1.0)


def exercise_failed_rewrite(binary: str) -> None:
    with tempfile.TemporaryDirectory(prefix="mix-rewrite-") as directory:
        with subprocess.Popen(
            [binary, "limited-rewrite"], cwd=directory, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, preexec_fn=ignore_file_size_signal,
        ) as process:
            try:
                assert process.stdout is not None
                readable, _, _ = select.select([process.stdout], [], [], 5)
                assert readable and process.stdout.readline().strip() == "First pass complete"
                resource.prlimit(process.pid, resource.RLIMIT_FSIZE, (256, 256))
                stdout, stderr = process.communicate("\n", timeout=20)
                assert process.returncode != 0 and "Mixing output failure:" in stderr, stdout + stderr
                payload = (Path(directory) / "fixture.mix").read_bytes()
                assert payload[:8] == b"GRMIXTMP" and len(payload) < 272
            finally:
                if process.poll() is None:
                    process.kill()


if __name__ == "__main__":
    command = sys.argv[1:]
    if len(command) == 1:
        exercise(command, "valid")
        exercise(command, "incomplete", "Cannot close an incomplete mixing stream")
        exercise(command, "limited", "Mixing output failure:")
        if hasattr(resource, "prlimit"):
            exercise_failed_rewrite(command[0])
    exercise(command, "order", "Out-of-order mixing block")
