"""Native sequential-record fixtures for the production RIS4 cache probe."""
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


def record(payload: bytes) -> bytes:
    marker = struct.pack("=i", len(payload))
    return marker + payload + marker


def main() -> None:
    executable = Path(sys.argv[1]).resolve()
    command = [str(executable)]
    if len(sys.argv) > 2:
        command = [sys.argv[2], "-n", "2", *command]
    end = record(struct.pack("=i", -1))
    ordinary = record(struct.pack("=id", 42, 0.125))
    cases = {"missing": None, "end": b"", "normal": ordinary, "legacy": end, "cutoff": ordinary,
             "incomplete": ordinary, "empty": b"", "truncated": b"\x04\x00",
             "short": record(b"\xff\xff"), "unreadable": None}
    with tempfile.TemporaryDirectory(prefix="ris4-cache-") as temporary:
        for whichone, name, contents in (
                (whichone, name, contents) for whichone in (1, 2) for name, contents in cases.items()):
            work = Path(temporary) / f"{whichone}-{name}"
            work.mkdir()
            filename = "cache.isotope.one-body.angular.cache" if whichone == 1 else "cache.isotope.two-body.angular.cache"
            if contents is not None:
                (work / filename).write_bytes(contents)
            if name == "unreadable":
                # A directory is deterministic even for privileged test users.
                (work / filename).mkdir()
            result = subprocess.run([*command, str(whichone), "stamp" if name in ("end", "normal") else name], cwd=work, capture_output=True, timeout=20)
            if name == "unreadable":
                assert result.returncode != 0, (name, result.stdout, result.stderr)
                assert b"ris4: angular cache" in result.stderr
                continue
            assert result.returncode == 0, (name, result.stdout, result.stderr)
            expected = b"T" if name in ("end", "normal") else b"F"
            flags = [line.split()[-1] for line in result.stdout.splitlines() if b"AVAILABLE" in line]
            assert flags == [expected, expected], (name, result.stdout)
    print("RIS4: one/two-body cache cases, each successful probe called twice, passed")


if __name__ == "__main__":
    main()
