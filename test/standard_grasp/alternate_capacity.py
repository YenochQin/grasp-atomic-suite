"""Check imported wavefunction interfaces against alternate suite capacities."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> None:
    compiler = sys.argv[1]
    suite = Path(__file__).resolve().parents[2]
    sources = [
        "lib/libmod/vast_kind_param_M.f90",
        "lib/libmod/suite_parameters_M.f90",
        "lib/libmod/parameter_def_M.f90",
        "lib/libmod/grid_C.f90",
        "lib/libmod/pote_C.f90",
        "lib/libmod/tatb_C.f90",
        "appl/rwfnestimate/tail_I.f90",
        "appl/rwfnestimate/solvh_I.f90",
        "appl/rwfnestimate/frmtfp_I.f90",
        "appl/rwfnestimate/frmrwf_I.f90",
        "appl/rwfnestimate/summry_I.f90",
    ]
    with tempfile.TemporaryDirectory(prefix="standard-grasp-capacity-") as temporary:
        work = Path(temporary)
        copies: list[str] = []
        for source in sources:
            path = suite / "src" / source
            text = path.read_text()
            if path.name == "suite_parameters_M.f90":
                if "NNNP = 590" not in text or "NNNW = 127" not in text:
                    raise AssertionError("alternate-capacity test expects pristine suite defaults")
                text = text.replace("NNNP = 590", "NNNP = 790").replace("NNNW = 127", "NNNW = 149")
            copy = work / path.name
            copy.write_text(text)
            copies.append(str(copy))
        for fixture in ["shared_grid_capacity.f90", "standard_grasp/interface_capacity.f90"]:
            binary = work / Path(fixture).stem
            result = subprocess.run(
                [compiler, "-std=f2018", "-fcheck=all", *copies, str(suite / "test" / fixture), "-o", str(binary)],
                cwd=work, capture_output=True, text=True, timeout=60, check=False,
            )
            assert result.returncode == 0, result.stderr
            subprocess.run([str(binary)], cwd=work, check=True, timeout=10)
        print("Ordinary GRASP radial/orbital interfaces: 790 / 800 / 149")


if __name__ == "__main__":
    main()
