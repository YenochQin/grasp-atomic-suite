"""Exercise MPI startup with only its output filename replaced by a test seam."""
from __future__ import annotations

import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> None:
    compiler, mpifort, mpiexec = sys.argv[1:4]
    root = Path(__file__).resolve().parents[2]
    text = (root / "src/lib/mpi/mpi_utilities.f90").read_text()
    declaration = re.search(r"(?im)^[ \t]*subroutine[ \t]+start_mpi_application[ \t]*\(", text)
    assert declaration is not None, "MPI startup routine declaration missing"
    start = declaration.start()
    ending = re.search(r"(?im)^\s*end(?:\s+subroutine(?:\s+start_mpi_application)?)?\s*$", text[start:])
    assert ending is not None, "MPI startup routine ending missing"
    body = text[start:start + ending.end()] + "\n"
    body, count = re.subn(r"file\s*=\s*(['\"])/dev/null\1", "file=output_path", body)
    assert count == 1, "Expected exactly one worker output filename"
    source = ("module startup_probe\nuse mpi_runtime, only: initialize_mpi_runtime\n"
              "implicit none\ncharacter(1024) :: output_path\ncontains\n" + body + "end module\n")
    driver = """program output_probe
use startup_probe
use mpi_f08
use, intrinsic :: iso_fortran_env, only: int64, output_unit
implicit none
integer::rank,ranks,length,status
integer(int64)::count
character(MPI_MAX_PROCESSOR_NAME)::host
character(16)::mode
call get_command_argument(1,mode)
call get_command_argument(2,output_path)
call start_mpi_application(rank,ranks,host,length,count,'probe','probe','none','none',mode/='visible')
write(output_unit,'(a,i0)') 'visible rank ',rank
call mpi_barrier(MPI_COMM_WORLD,status)
call mpi_finalize(status)
end program
"""
    env = dict(os.environ, OMPI_FC=compiler)
    with tempfile.TemporaryDirectory(prefix="mpi-output-") as temporary:
        work = Path(temporary)
        (work / "probe.f90").write_text(source + driver)
        subprocess.run([mpifort, "-std=f2018", "-fcheck=all",
                        str(root / "src/lib/core/run_timer.f90"),
                        str(root / "src/lib/mpi/mpi_runtime.f90"),
                        str(root / "src/lib/mpi/mpi_run_logging.f90"), "probe.f90", "-o", "probe"],
                       cwd=work, env=env, check=True, capture_output=True, text=True)
        for mode, path in (("quiet", "/dev/null"), ("failure", str(work)), ("visible", str(work))):
            result = subprocess.run([mpiexec, "-n", "2", str(work / "probe"), mode, path],
                                    cwd=work, env=env, capture_output=True, text=True, timeout=20)
            assert result.returncode == 0, (mode, result.stdout, result.stderr)
            assert "visible rank 0" in result.stdout
            assert result.stdout.count("Date and Time:") == 1, result.stdout
            assert result.stdout.count("Participating nodes:") == 1, result.stdout
            if mode == "failure":
                assert "MPI worker output warning: rank=1; iostat=" in result.stderr, result.stderr
            else:
                assert "MPI worker output warning" not in result.stderr, result.stderr
            if mode == "quiet":
                assert "visible rank 1" not in result.stdout, result.stdout
            if mode == "visible":
                assert "visible rank 1" in result.stdout, result.stdout
    print("MPI worker output: quiet, visible and nonfatal open failure passed on two ranks")


if __name__ == "__main__":
    main()
