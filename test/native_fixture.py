"""Offline test-fixture conversion; never imported by production programs.

Set RCSFS_TEST_CLI to a native rcsfs executable, or put it on PATH.
The reader under test never invokes this helper or consumes the text fixture.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def descriptors(text: Path, stem: str | None = None) -> None:
    stem = stem or text.stem
    output = text.parent / f"{stem}_descriptors.parquet"
    if output.exists():
        return
    executable = os.environ.get("RCSFS_TEST_CLI") or shutil.which("rcsfs")
    if not executable:
        raise RuntimeError("Set RCSFS_TEST_CLI to the native fixture producer")
    with tempfile.TemporaryDirectory(prefix="csf-fixture-") as temporary:
        work = Path(temporary)
        for args in [
            ["convert", str(text.resolve()), f"{stem}.parquet", "--num-workers", "1"],
            [
                "gen-descriptors",
                f"{stem}.parquet",
                f"{stem}_descriptors.parquet",
                "--header",
                f"{stem}_header.toml",
                "--num-workers",
                "1",
            ],
        ]:
            result = subprocess.run(
                [executable, *args],
                cwd=work,
                check=False,
                capture_output=True,
                text=True,
                timeout=90,
            )
            if result.returncode:
                raise RuntimeError(result.stdout + result.stderr)
        shutil.copyfile(work / output.name, output)
        shutil.copyfile(
            work / f"{stem}_header.toml", text.parent / f"{stem}_header.toml"
        )
