"""Compare serial RIS4 reports from freshly computed and cached Ni coefficients."""

import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rdensity.regression import run


def main() -> None:
    command = [str(Path(sys.argv[1]).resolve())]
    data = Path(sys.argv[2])
    for higher_order in ("n", "y"):
        with tempfile.TemporaryDirectory(prefix="ris4-cache-reuse-") as temporary:
            work = Path(temporary)
            for name in ("isodata", "test.c", "test.w", "test.m"):
                shutil.copyfile(data / name, work / name)

            answers = f"y\ntest\nn\n{higher_order}\n"
            run(command, work, answers + "y\n")
            caches = {path: path.read_bytes() for path in (work / "test.IOB", work / "test.ITB")}
            report = work / "test.i"
            expected = report.read_bytes()
            report.unlink()

            output = run(command, work, answers)
            assert "One-body angular file available" in output, output
            assert "Two-body angular file available" in output, output
            assert report.read_bytes() == expected, "Cache reuse changed the RIS4 report"
            for path, contents in caches.items():
                assert path.read_bytes() == contents, f"Cache reuse changed {path.name}"


if __name__ == "__main__":
    main()
