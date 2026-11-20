"""Exercise FICAL with the fixed-width rows written by RIS4 and RIS4_MPI."""

from pathlib import Path
import math
import subprocess
import sys
import tempfile
from typing import NamedTuple


class Level(NamedTuple):
    number: int
    energy: float
    nms: tuple[float, float, float]
    sms: tuple[float, float, float]
    sms_unit: str


UPPER = (
    Level(101, -9.0, (0.5, -0.1, 0.4), (50.0, 10.0, 60.0), "MHz"),
    Level(102, -8.0, (0.7, -0.2, 0.5), (0.4, -0.1, 0.3), "GHz"),
)
LOWER = (
    Level(201, -10.0, (0.2, 0.02, 0.22), (0.2, 0.02, 0.22), "GHz"),
    Level(202, -9.5, (0.3, -0.02, 0.28), (-75.0, 25.0, -50.0), "MHz"),
)
MASS_COLUMNS = (0, 2)  # Nonrelativistic and total relativistic terms.
SMS_GHZ = {101: (0.05, 0.06), 102: (0.4, 0.3), 201: (0.2, 0.22), 202: (-0.075, -0.05)}


def fields(*values: float) -> str:
    return "".join(f"{value:20.10E}".replace("E", "D") for value in values)


def ris4_report(path: Path, levels: tuple[Level, ...], *, fitted: bool) -> None:
    def prefix(level: Level) -> str:
        return f" {level.number:3d}     {'   0'}{' +  '}   "

    lines = [f"{'Number of eigenvalues: ':>24}{len(levels):3d}", "", "", " Level  J Parity  Energy"]
    lines.extend(prefix(level) + fields(level.energy) + "  (a.u.)" for level in levels)
    for name in ("Normal", "Specific"):
        lines.extend(("", "", f" Level  J Parity  {name} mass shift parameter"))
        for level in levels:
            lines.extend(("", "                             <K^1>             <K^2+K^3>          <K^1+K^2+K^3>"))
            values = level.nms if name == "Normal" else level.sms
            unit = "GHz" if name == "Normal" else level.sms_unit
            # The unconverted row and the converted row occupy separate records.
            lines.append(prefix(level) + fields(*values) + "  (a.u.)")
            label = f"  ({unit} u)" if unit else ""
            lines.append(" " * 20 + fields(*values) + label)
    if fitted:
        lines.extend(("", "", " Level  J Parity  Electron density in atomic units", "", " Dens. (a.u.)"))
    else:
        lines.extend(("", "", " Electron density in atomic units", "", " Level  J Parity        DENS (a.u.)", ""))
    lines.extend(prefix(level) + fields(0.3) for level in levels)
    if fitted:
        lines.extend(("", "", " Field shift electronic factors", "", " F0 (GHz/fm^2)"))
        for level in levels:
            f0 = level.nms[0]
            lines.append(prefix(level) + fields(f0, 2*f0, 3*f0, 4*f0) + f"{0.0:10.4f}")
        lines.extend(("", "", " Corrected field shift electronic factors", "", " F0VED0 (GHz/fm^2)"))
        lines.extend(prefix(level) + fields(level.nms[0], 2*level.nms[0]) for level in levels)
    path.write_text("\n".join(lines) + "\n")


def run_fical(binary: Path, directory: Path, unit: str, *, ci: bool = False, relativistic: bool = False,
              fitted: bool = False, upper: tuple[Level, ...] = UPPER) -> tuple[subprocess.CompletedProcess[str], Path]:
    extension = "ci" if ci else "i"
    ris4_report(directory / f"upper.{extension}", upper, fitted=fitted)
    ris4_report(directory / f"lower.{extension}", LOWER, fitted=fitted)
    (directory / "isodata").write_text("Z\n10.0\nlabel\nlabel\nlabel\n0.52\nlabel\n3.00\n")
    answers = ["n", "upper", "lower", "y" if ci else "n", "y" if fitted else "n", "para", unit,
               "y" if relativistic else "n"]
    result = subprocess.run(
        [str(binary)], cwd=directory, input="\n".join(answers) + "\n", capture_output=True, text=True, timeout=10,
    )
    return result, directory / f"upper.lower.{'cfi' if ci else 'fi'}"


def check_mixed_units(binary: Path) -> None:
    cases = (
        ("GHz", 1.0, False, False, False),
        ("MHz", 1000.0, True, False, False),
        ("meV", 1 / 241.798926, False, True, True),
    )
    for unit, scale, ci, relativistic, fitted in cases:
        with tempfile.TemporaryDirectory(prefix="fical-mixed-units-") as work:
            result, report = run_fical(binary, Path(work), unit, ci=ci, relativistic=relativistic, fitted=fitted)
            assert result.returncode == 0 and "program FICAL finished" in result.stdout, result.stderr + result.stdout
            assert result.stderr.count("Name of state:") == 2
            assert "Upper state" in result.stderr and "Lower state" in result.stderr
            block, field_block = report.read_text().split("LINE MASS SHIFT PARAMETERS", 1)[1].split(
                "LINE FIELD SHIFT PARAMETERS", 1,
            )
            assert f"SMS ({unit} u)" in block, block
            rows = [line for line in block.splitlines() if line[:4].strip().isdigit()]
            assert len(rows) == len(UPPER) * len(LOWER), rows
            column = MASS_COLUMNS[relativistic]
            for row, (upper, lower) in zip(rows, ((a, b) for a in UPPER for b in LOWER), strict=True):
                assert (int(row[0:4]), int(row[17:20])) == (upper.number, lower.number), row
                nms = (upper.nms[column] - lower.nms[column]) * scale
                sms = (SMS_GHZ[upper.number][relativistic] - SMS_GHZ[lower.number][relativistic]) * scale
                observed = [float(row[start:start + 18].replace("D", "E")) for start in (62, 80, 98)]
                for actual, expected in zip(observed, (nms, sms, nms + sms), strict=True):
                    assert math.isclose(actual, expected, rel_tol=1e-7, abs_tol=1e-12), (unit, row, actual, expected)

            first_field_row = next(line for line in field_block.splitlines() if line[:4].strip().isdigit())
            if fitted:
                f0 = float(first_field_row[44:62].replace("D", "E"))
                assert math.isclose(f0, 0.3 * scale, rel_tol=1e-7), first_field_row
                assert "NaN" not in field_block, field_block
            else:
                assert "NaN" in first_field_row, first_field_row


def check_invalid_units(binary: Path) -> None:
    for label in ("THz", ""):
        with tempfile.TemporaryDirectory(prefix="fical-invalid-unit-") as work:
            first = UPPER[0]._replace(sms_unit=label)
            result, report = run_fical(binary, Path(work), "GHz", upper=(first, UPPER[1]))
            assert result.returncode != 0, result.stdout
            assert "unsupported SMS unit in upper.i" in result.stderr, result.stderr
            assert not report.exists(), report


if __name__ == "__main__":
    tool = Path(sys.argv[1]).resolve()
    check_mixed_units(tool)
    check_invalid_units(tool)
    print("FICAL mixed MHz/GHz records, conversions, and unit-error handling passed")
