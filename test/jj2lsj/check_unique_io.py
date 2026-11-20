"""Exercise the complete production uniquelsj workspace, without resizing it."""
from pathlib import Path
import subprocess
import sys
import tempfile


def level(number: int, j: str = "0") -> str:
    return f"{number:3d} {j:>4}     +        {-number:16.9f}     {100:7.3f}%\n"


def component(mix: float, weight: float, label: str) -> str:
    return f"       {mix:12.8f}   {weight:11.8f}   {label}\n"


def main() -> None:
    if sys.argv[1] == "--compile":
        compiler, destination = sys.argv[2:4]
        root = Path(__file__).resolve().parents[2]
        work = Path(destination).resolve()
        work.mkdir(parents=True, exist_ok=True)
        source = (root / "src/appl/jj2lsj/jj2lsj_code.f90").read_text()
        start = source.index("      subroutine uniquelsj\n")
        end = source.index("      end subroutine uniquelsj", start) + len("      end subroutine uniquelsj")
        # Only the private routine's host is reduced. The entire procedure,
        # including its production capacities and algorithm, is copied verbatim.
        (work / "unique.f90").write_text(
            "module jj2lsj_code\nuse iso_fortran_env, only: real64\n"
            "use jj2lsj_profile\nimplicit none\ncontains\n"
            + source[start:end] + "\nend module\n")
        subprocess.run([compiler, "-std=f2018", "-fcheck=all", "-fbacktrace",
                        str(root / "src/lib/grasp_shared/input_error_handling.f90"),
                        str(root / "src/appl/jj2lsj/jj2lsj_profile.f90"),
                        "unique.f90", str(root / "test/jj2lsj/unique_io.f90"),
                        "-o", "unique"], cwd=work, check=True)
        return
    executable = Path(sys.argv[1]).resolve()
    header = " Pos   J   Parity      Energy Total      Comp. of ASF\n"
    one = level(1) + component(1.0, 1.0, "label A")
    # Collision requires reassignment, preserving the original label algorithm.
    collision = (level(1) + component(0.8, 0.64, "label A")
                 + component(0.6, 0.36, "label B")
                 + level(2) + component(0.9, 0.81, "label A")
                 + component(0.4, 0.16, "label C"))
    valid = {"one": one, "collision": collision,
             "groups": one + " \n" + level(2, "1") + component(1.0, 1.0, "label B")}
    invalid = {"empty": "", "header_only": header,
               "bad_level": header + level(1).replace("-1.000", "xx.000"),
               "no_components": header + level(1),
               "empty_group": header + level(1) + " \n",
               "bad_component": header + level(1) + "       invalid number\n",
               "short_component": header + level(1) + "       1.0\n",
               "zero_weight": header + level(1) + component(1.0, 0.0, "label A"),
               "nan_weight": header + level(1) + component(1.0, float('nan'), "label A")}
    invalid["exhausted_labels"] = header + one + level(2) + component(1.0, 1.0, "label A")
    invalid["zero_after_collision"] = (header + one + level(2)
                                       + component(0.5, 0.25, "label A")
                                       + component(0.0, 0.0, "label B"))
    with tempfile.TemporaryDirectory(prefix="unique-io-") as temporary:
        root = Path(temporary)
        for name, body in valid.items():
            work = root / name
            work.mkdir()
            (work / "input.lbl").write_text(header + body)
            result = subprocess.run([executable], cwd=work, capture_output=True, timeout=30)
            assert result.returncode == 0, (name, result.stdout, result.stderr)
            output = (work / "unique.lbl").read_text()
            assert "label A" in output
            if name == "collision":
                assert "label B" in (work / "unique.summary").read_text()
            reference = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).with_name("unique_expected")
            for filename in ("unique.lbl", "unique.summary"):
                target = reference / f"{name}-{filename}"
                actual = (work / filename).read_bytes()
                assert actual == target.read_bytes(), (name, filename)
        for name, body in invalid.items():
            work = root / name
            work.mkdir()
            (work / "input.lbl").write_text(body)
            result = subprocess.run([executable], cwd=work, capture_output=True, timeout=30)
            assert result.returncode != 0, (name, "invalid input accepted")
            assert b"uniquelsj:" in result.stderr, (name, result.stderr)
    print(f"uniquelsj: {len(valid)} normal and {len(invalid)} malformed complete-workspace cases passed")


if __name__ == "__main__":
    main()
