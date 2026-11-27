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
            for name in ("isodata", "test.mix", "test.wfn", "test_desc.parquet", "test_header.toml"):
                shutil.copyfile(data / name, work / name)

            answers = f"y\ntest\nn\n{higher_order}\n"
            run(command, work, answers + "y\n")
            caches = {path: path.read_bytes() for path in (work / "test.isotope.one-body.angular.cache", work / "test.isotope.two-body.angular.cache")}
            report = work / "test.mcdhf.isotope-shift"
            expected = report.read_bytes()
            report.unlink()

            output = run(command, work, answers)
            assert "One-body angular file available" in output, output
            assert "Two-body angular file available" in output, output
            assert report.read_bytes() == expected, "Cache reuse changed the RIS4 report"
            for path, contents in caches.items():
                assert path.read_bytes() == contents, f"Cache reuse changed {path.name}"

            legacy_caches = {
                path.rename(path.with_name(path.name.replace(".isotope.", ".ris4."))): contents
                for path, contents in caches.items()
            }
            report.unlink()
            output = run(command, work, answers)
            assert "One-body angular file available" in output, output
            assert "Two-body angular file available" in output, output
            assert report.read_bytes() == expected
            assert all(path.read_bytes() == contents for path, contents in legacy_caches.items())


if __name__ == "__main__":
    main()
