"""Offline test-fixture conversion; never imported by production programs.

Set RCSFS_TEST_CLI to an rcsfs executable that writes descriptor format 3.
The reader under test never invokes this helper or consumes the text fixture.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from struct import pack
import shutil
import subprocess
import tempfile
from pathlib import Path


def descriptors(text: Path, stem: str | None = None) -> None:
    stem = stem or text.stem
    output = text.parent / f"{stem}_desc.parquet"
    header = text.parent / f"{stem}_header.toml"
    if output.is_file() and header.is_file():
        return
    executable = os.environ.get("RCSFS_TEST_CLI") or shutil.which("rcsfs")
    if not executable:
        raise RuntimeError("Set RCSFS_TEST_CLI to the V3 native fixture producer")
    executable = str(Path(executable).resolve(strict=True))
    with tempfile.TemporaryDirectory(prefix="csf-fixture-") as temporary:
        work = Path(temporary)
        result = subprocess.run(
            [executable, "convert", str(text.resolve()), output.name],
            cwd=work,
            check=False,
            capture_output=True,
            text=True,
            timeout=90,
        )
        if result.returncode:
            raise RuntimeError(result.stdout + result.stderr)
        shutil.copyfile(work / header.name, header)
        shutil.copyfile(work / output.name, output)

def write_radial_stream(
    path: Path, grid: Sequence[float],
    orbitals: Sequence[tuple[int, int, float, float, Sequence[float], Sequence[float]]],
) -> None:
    with path.open("wb") as stream:
        stream.write(pack("<8sii", b"RAD--WFN", len(orbitals), len(grid)))
        for n, kappa, energy, origin, p, q in orbitals:
            assert len(p) == len(q) <= len(grid)
            stream.write(pack("<4i2d", n, kappa, len(p), 0, energy, origin))
        stream.write(pack(f"<{len(grid)}d", *grid))
        for _, _, _, _, p, q in orbitals:
            for component in (p, q):
                stream.write(pack(f"<{len(component)}d", *component))
                stream.write(bytes(8 * (len(grid) - len(component))))
