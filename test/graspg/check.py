"""Small file-format regression against the unmodified GRASPG e157d9f outputs."""
from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main() -> None:
    binaries = Path(sys.argv[1]).resolve()
    fixtures = Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="graspg-formats-") as temporary:
        work = Path(temporary)

        def run(name: str, data: str = "") -> str:
            result = subprocess.run(
                [str(binaries / name)], input=data, cwd=work, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30,
                check=False,
            )
            assert result.returncode == 0, result.stdout
            return result.stdout

        output = run("rcsfggenerate_csfg", (fixtures / "generate.inp").read_text())
        assert "Execution complete" in output, output
        for actual, expected in [("rcsfg.out", "expected.g"), ("rlabel.out", "expected.l")]:
            assert (work / actual).read_bytes() == (fixtures / expected).read_bytes(), actual
        shutil.copyfile(work / "rcsfg.out", work / "rcsfg.inp")
        shutil.copyfile(work / "rlabel.out", work / "rlabel.inp")
        run("rcsfgexpand_csfg")
        assert (work / "rcsf.out").read_bytes() == (fixtures / "expected.c").read_bytes()
        # Every candidate is in the reference; reduction must preserve all CSFs.
        shutil.copyfile(work / "rcsfg.inp", work / "rcsfgmr.inp")
        shutil.copyfile(work / "rlabel.inp", work / "rcsfgmr.l")
        (work / "rcsfg.out").unlink()
        output = run("rcsfginteract_csfg", "1\n")
        assert (work / "rcsfg.out").read_bytes() == (fixtures / "expected.g").read_bytes(), output
        # A smaller reference exercises the changed angular-interaction path.
        shutil.copyfile(fixtures / "reference.g", work / "rcsfgmr.inp")
        (work / "rcsfg.out").unlink()
        output = run("rcsfginteract_csfg", "1\n")
        assert (work / "rcsfg.out").read_bytes() == (fixtures / "expected-interact.g").read_bytes(), output
        print("CSFG generation/expansion: 6 CSFs; two-CSF reference retains 5 CSFs")


if __name__ == "__main__":
    main()
