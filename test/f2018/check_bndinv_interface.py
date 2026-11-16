"""Compare BNDINV interfaces to the complete implementations via procedure association."""
from __future__ import annotations

from pathlib import Path
import re
import subprocess
import sys
import tempfile


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    for variant in ("rbiotransform", "rbiotransform_mpi"):
        app = root / "src/appl" / variant
        implementation = (app / "bndinv.f90").read_text()
        interface = (app / "bndinv_I.f90").read_text()
        # GNU's procedure association check does not compare explicit bounds.
        # Check these separately instead of treating compilation as proof.
        for name in ("a", "el"):
            pattern = r"::\s*" + name + r"\(([^)]+)\)"
            actual = re.search(pattern, implementation, re.I)
            declared = re.search(pattern, interface, re.I)
            assert actual is not None and declared is not None
            assert actual.group(1).replace(" ", "") == declared.group(1).replace(" ", ""), (variant, name)
        source = ("module implementation_contract\nuse bndinv_i, only: contract => bndinv\ncontains\n"
                  + implementation + "\nsubroutine verify\nprocedure(contract), pointer :: p\n"
                  "p => bndinv\nend subroutine\nend module\n")
        with tempfile.TemporaryDirectory(prefix="bndinv-interface-") as temporary:
            work = Path(temporary)
            (work / "implementation.f90").write_text(source)
            subprocess.run([sys.argv[1], "-std=f2018", "-pedantic-errors", "-fimplicit-none",
                            "-fsyntax-only", str(app / "bndinv_I.f90"), "implementation.f90"],
                           cwd=work, check=True, capture_output=True, text=True)
        print(f"{variant}: actual implementation matches public BNDINV procedure contract")


if __name__ == "__main__":
    main()
