"""Run the Suite CMake ABI probes against same-symbol LP64 and ILP64 fixtures."""

from pathlib import Path
import subprocess
import sys
import tempfile


def make_archive(work: Path, compiler: str, archiver: str, width: int) -> Path:
    fixture = Path(__file__).with_name("blas_lapack_width_fixture.f90")
    obj = work / f"width{width}.o"
    archive = work / f"libwidth{width}.a"
    flags = ["-fdefault-integer-8"] if width == 8 else []
    subprocess.run([compiler, "-std=f2018", *flags, "-c", str(fixture), "-o", str(obj)], cwd=work, check=True)
    subprocess.run([archiver, "rcs", str(archive), str(obj)], cwd=work, check=True)
    return archive


def main() -> None:
    cmake, compiler, archiver = sys.argv[1:]
    module = Path(__file__).resolve().parents[2] / "cmake/BlasLapackLP64.cmake"
    with tempfile.TemporaryDirectory(prefix="suite-lp64-") as temporary:
        work = Path(temporary)
        lp64 = make_archive(work, compiler, archiver, 4)
        ilp64 = make_archive(work, compiler, archiver, 8)
        project = work / "source"
        project.mkdir()
        (project / "CMakeLists.txt").write_text(
            'cmake_minimum_required(VERSION 3.31)\n'
            'project(blas_lapack_width LANGUAGES Fortran)\n'
            'list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")\n'
            'include("${SUITE_ABI_MODULE}")\n'
        )
        (project / "FindBLAS.cmake").write_text(
            'if(NOT BLA_SIZEOF_INTEGER STREQUAL "4")\n'
            '    message(FATAL_ERROR "Expected the Suite to request 32-bit BLAS integers")\n'
            'endif()\n'
            'add_library(BLAS::BLAS INTERFACE IMPORTED)\n'
            'target_link_libraries(BLAS::BLAS INTERFACE "${SUITE_TEST_BLAS}")\n'
        )
        (project / "FindLAPACK.cmake").write_text(
            'if(NOT BLA_SIZEOF_INTEGER STREQUAL "4")\n'
            '    message(FATAL_ERROR "Expected the Suite to request 32-bit LAPACK integers")\n'
            'endif()\n'
            'add_library(LAPACK::LAPACK INTERFACE IMPORTED)\n'
            'target_link_libraries(LAPACK::LAPACK INTERFACE "${SUITE_TEST_LAPACK}" BLAS::BLAS)\n'
        )

        def configure(blas: Path, lapack: Path, *, size: str = "4", flags: str = "") -> subprocess.CompletedProcess[str]:
            # Reuse the successful configure directory to check that an old
            # affirmative result cannot mask a replacement library.
            return subprocess.run([
                cmake, "-S", str(project), "-B", str(work / "build"),
                f"-DSUITE_ABI_MODULE={module}", f"-DCMAKE_Fortran_COMPILER={compiler}",
                "-DCMAKE_BUILD_TYPE=Release", f"-DBLA_SIZEOF_INTEGER={size}",
                f"-DCMAKE_Fortran_FLAGS={flags}", f"-DSUITE_TEST_BLAS={blas}", f"-DSUITE_TEST_LAPACK={lapack}",
            ], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=45)

        normal = configure(lp64, lp64, size="ANY")
        if normal.returncode or any(f"Suite {kind} LP64 ABI verified" not in normal.stdout for kind in ("BLAS", "LAPACK")):
            raise RuntimeError(normal.stdout)
        print("LP64 BLAS and LAPACK accepted")
        for name, recognition, result in (
            ("requested ILP64", "Requested BLA_SIZEOF_INTEGER=8", configure(lp64, lp64, size="8")),
            ("64-bit Suite integers", "32-bit default integers", configure(lp64, lp64, flags="-fdefault-integer-8")),
            ("same-symbol ILP64 BLAS", "DGEMV did not accept the LP64 argument layout", configure(ilp64, lp64)),
            ("same-symbol ILP64 LAPACK", "DSPEVX wrote outside the 32-bit integer outputs", configure(lp64, ilp64)),
        ):
            if result.returncode == 0 or recognition not in result.stdout:
                raise RuntimeError(f"{name}: expected configuration rejection with {recognition!r}\n{result.stdout}")
            print(f"{name} rejected")


if __name__ == "__main__":
    main()
