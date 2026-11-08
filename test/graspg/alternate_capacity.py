"""Compile production module/interface copies at non-default capacities."""
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
        "lib/libmod/grid_C.f90", "lib/libmod/def_C.f90", "lib/libmod/kkstart_C.f90",
        "lib/libmod_csfg/csfg_parameters_M.f90", "lib/libmod_csfg/csfg_scf_C.f90",
        "lib/libmod_csfg/csfg_tv_C.f90", "lib/libmod_csfg/symexpand_mod.f90",
        "lib/libmod_csfg/symmatrix_mod.f90", "lib/libmod_csfg/symmatrix_restart_C.f90",
        "lib/libmod_csfg/csfg_decide_C.f90",
        "appl/rwfnestimate90_csfg/tail_I.f90", "appl/rwfnestimate90_csfg/frmtfp_I.f90",
    ]
    with tempfile.TemporaryDirectory(prefix="csfg-capacity-") as temporary:
        work = Path(temporary)
        copies = []
        for source in sources:
            path = suite / "src" / source
            text = path.read_text()
            if path.name == "suite_parameters_M.f90":
                text = text.replace("NNNP = 590", "NNNP = 790").replace("NNNW = 127", "NNNW = 149")
            elif path.name == "csfg_parameters_M.f90":
                text = text.replace("CSFG_MAX_BLOCKS = 50", "CSFG_MAX_BLOCKS = 64")
            copy = work / path.name
            copy.write_text(text)
            copies.append(str(copy))
        for name in ["capacity", "interface_capacity"]:
            binary = work / name
            result = subprocess.run(
                [compiler, "-std=f2018", *copies, str(suite / "test/graspg" / f"{name}.f90"), "-o", str(binary)],
                cwd=work, capture_output=True, text=True, timeout=60, check=False,
            )
            assert result.returncode == 0, result.stderr
            result = subprocess.run([str(binary)], cwd=work, capture_output=True, text=True, check=True)
            if name == "capacity":
                assert result.stdout.split()[-4:] == ["790", "800", "149", "64"], result.stdout
        print("Alternate capacities and wavefunction interfaces: 790 / 800 / 149 / 64")


if __name__ == "__main__":
    main()
