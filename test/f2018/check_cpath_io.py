"""Exercise CPATH with a defined MPI_TMP and initially unallocated helper values."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory(prefix="cpath-io-") as temporary:
    work = Path(temporary)
    rank_root = work/"ranks with spaces"
    rank_root.mkdir()
    result = subprocess.run([sys.argv[1], "-n", "2", str(Path(sys.argv[2]).resolve())],
                            cwd=work, capture_output=True, text=True, timeout=30,
                            env=os.environ | {"MPI_TMP": str(rank_root)})
    assert result.returncode == 0, result.stdout+result.stderr
    assert result.stdout.count("PASS CPATH input environment") == 2, result.stdout
    assert all((rank_root/f"{rank:03d}").is_dir() for rank in range(2))
    print("CPATH created both rank directories and retained the original input directory")
