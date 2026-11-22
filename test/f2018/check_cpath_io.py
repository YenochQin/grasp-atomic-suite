"""Check MPI directory setup and run summaries with shared and split scratch roots."""
from __future__ import annotations

import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def run(work: Path, roots: list[Path], *, use_disks: bool) -> list[str]:
    command = [sys.argv[1], "-n", str(len(roots)), str(Path(sys.argv[2]).resolve())]
    if use_disks:
        (work / "disks").write_text("".join(f'"{path}"\n' for path in [work, *roots]))
        command.extend(str(root) for root in roots)

    result = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=30,
                            env=os.environ | {"MPI_TMP": str(roots[0])})
    assert result.returncode == 0, result.stdout + result.stderr
    assert all((root / f"{rank:03d}").is_dir() for rank, root in enumerate(roots))

    lines = result.stdout.splitlines()
    assert lines.count(f"Start Dir: {work}") == 1, result.stdout
    assert lines.count(f"Serial I/O Dir (rank 0): {work}") == 1, result.stdout
    node_headers = [line for line in lines if line.startswith("Participating nodes:")]
    assert len(node_headers) == 1, result.stdout
    node_header = re.fullmatch(r"Participating nodes: hosts=(\d+), MPI ranks=(\d+)", node_headers[0])
    assert node_header and int(node_header[2]) == len(roots), result.stdout
    host_counts = re.findall(r"^  .+: ranks=(\d+)$", result.stdout, re.MULTILINE)
    assert len(host_counts) == int(node_header[1]) and sum(map(int, host_counts)) == len(roots), result.stdout

    timestamp_pattern = r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d\.\d{3} [+-]\d{4}"
    for label in ("Date and Time:", "Finish Date and Time:"):
        timestamps = [line for line in lines if line.startswith(label)]
        assert len(timestamps) == 1, result.stdout
        assert re.fullmatch(f"{label} {timestamp_pattern}", timestamps[0]), result.stdout
    wall_seconds = re.findall(r"^Wall time: (\d+) seconds$", result.stdout, re.MULTILINE)
    assert len(wall_seconds) == 1 and int(wall_seconds[0]) >= (2 if use_disks else 0), result.stdout
    assert lines.count("CPATH_IO: Execution complete.") == 1, result.stdout
    assert "mpi stopped by node-" not in result.stdout, result.stdout
    return lines


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="cpath-io-") as temporary:
        work = Path(temporary).resolve()
        root_a, root_b = work / "scratch A with spaces", work / "scratch B with spaces"
        root_a.mkdir()
        root_b.mkdir()

        lines = run(work, [root_a] * 2, use_disks=False)
        assert lines.count(f"Work Dir (Parallel I/O): {root_a}") == 1, lines

        lines = run(work, [root_a, root_a, root_b, root_a], use_disks=True)
        start = lines.index("Work Dir (Parallel I/O):")
        assert lines[start + 1:start + 3] == [f"  {root_a} (ranks 0-1,3)", f"  {root_b} (ranks 2)"], lines
        print("MPI startup and completion summaries passed for shared and split scratch roots")


if __name__ == "__main__":
    main()
