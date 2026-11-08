"""Negative checks for the acceptance gate; no Fortran calculations."""
from pathlib import Path
import tempfile
import unittest

from check import check_commands, check_junit, check_roster, requirements


class AcceptanceTests(unittest.TestCase):
    def test_roster_requires_each_configured_group(self) -> None:
        for features in (set(), {"serial"}, {"mpi"}, {"serial", "mpi", "csfg"}):
            targets, tests = requirements(features)
            check_roster(features, targets, tests)
            with self.assertRaises(ValueError):
                check_roster(features, targets, tests - {"lib9290_quad"})
            with self.assertRaises(ValueError):
                check_roster(features, targets - {"mod"}, tests)

    def test_disabled_targets_and_empty_tests_are_rejected(self) -> None:
        targets, tests = requirements(set())
        with self.assertRaises(ValueError):
            check_roster(set(), targets | {"rmcdhf_mpi"}, tests)
        with self.assertRaises(ValueError):
            check_roster(set(), targets, set())

    def test_every_compile_command_is_checked(self) -> None:
        compiler = Path("/compiler/gfortran-15")
        flags = ["-std=f2018", "-pedantic-errors", "-Wall", "-Wextra", "-Werror=implicit-interface",
                 "-Werror=implicit-procedure", "-fno-fast-math", "-ffp-contract=off", "-O0", "-g",
                 "-fcheck=all", "-fbacktrace", "-ffpe-trap=invalid,zero,overflow"]
        good = {"file": "source.f90", "arguments": [str(compiler), *flags]}
        check_commands([good], "f2018-check", compiler)
        for bad_flag in ("-fallow-argument-mismatch", "-fno-automatic", "-std=legacy", "-O3",
                         "-ffast-math", "-fdefault-integer-8", "-ffp-contract=fast", "-Wno-implicit-interface",
                         "-fcheck=none", "-ffpe-trap=none"):
            with self.subTest(flag=bad_flag), self.assertRaises(ValueError):
                check_commands([good, {**good, "arguments": [*good["arguments"], bad_flag]}], "f2018-check", compiler)
        with self.assertRaises(ValueError):
            check_commands([{**good, "arguments": [str(compiler), "-std=f2018"]}], "f2018-check", compiler)
        with self.assertRaises(ValueError):
            check_commands([{**good, "arguments": ["/compiler/gfortran-16", *flags]}], "f2018-check", compiler)

    def test_skips_failures_missing_and_duplicate_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tests.xml"
            path.write_text('<testsuite><testcase name="a"/></testsuite>')
            self.assertEqual(check_junit(path, {"a"}), {"a": "passed"})
            for content in ('<testcase name="a"><skipped/></testcase>',
                            '<testcase name="a"><failure/></testcase>', '',
                            '<testcase name="a"/><testcase name="a"/>'):
                path.write_text(f"<testsuite>{content}</testsuite>")
                with self.assertRaises(ValueError):
                    check_junit(path, {"a"})


if __name__ == "__main__":
    unittest.main()
