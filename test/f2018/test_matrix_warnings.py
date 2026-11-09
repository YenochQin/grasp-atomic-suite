"""Negative tests for matrix coverage and lossless warning budgets; no SCF."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from matrix import cases, configure_command
from warning_inventory import (category, command_parts, compare, digest, inventory,
                               make_baseline, parse_diagnostics, record_path, write_json)


class MatrixWarningTests(unittest.TestCase):
    def test_cartesian_matrix_and_explicit_switches(self) -> None:
        configurations = cases()
        self.assertEqual(len(configurations), 32)
        self.assertEqual(len({case["id"] for case in configurations}), 32)
        for profile in ("f2018-check", "f2018-release"):
            self.assertEqual(len({tuple(case["features"].values()) for case in configurations
                                  if case["profile"] == profile}), 16)
        for case in configurations:
            command = configure_command(case, Path("/build space"), Path("/compiler/gfortran-15"),
                                        Path("/mpi/fort"), Path("/mpi/exec"))
            for axis, flag in (("serial", "GRASP_BUILD_SERIAL_APPS"), ("mpi", "GRASP_ENABLE_MPI"),
                               ("csfg", "GRASP_BUILD_CSFG"), ("openmp", "GRASP_JJ2LSJ_OPENMP")):
                self.assertIn(f"-D{flag}={'ON' if case['features'][axis] else 'OFF'}", command)

    def test_parser_never_silently_discards_output(self) -> None:
        self.assertEqual(parse_diagnostics("[]\n[]\n"), [])
        self.assertEqual(len(parse_diagnostics('[{"kind":"warning","message":"one"}]\n[]')), 1)
        for text in ('warning text', '[{"message":"no kind"}]', '{}', '[]\ntrailing output'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_diagnostics(text)

    def observation(self) -> dict:
        return {"case": "check-s1-m1-c1-o1", "profile": "f2018-check",
                "features": {key: True for key in ("serial", "mpi", "csfg", "openmp")},
                "compiler_version": "15.3.0", "compilation_count": 2,
                "target_compilation_counts": {"demo": 2},
                "files": {"src/a.f90": "a", "src/quiet.f90": "b"},
                "diagnostics": {"known": {"file": "src/a.f90", "target": "demo", "line": 3, "column": 1,
                    "option": "-Wmaybe-uninitialized", "message": "x may be used uninitialized",
                    "category": "correctness-review"}},
                "counts": {"known": 1}, "category_counts": {"correctness-review": 1}}

    def test_known_risk_is_not_certified_as_resolved(self) -> None:
        observation = self.observation()
        result = compare(observation, make_baseline([observation]))
        self.assertEqual(result["status"], "passed")
        self.assertFalse(result["correctness_warnings_clear"])

    def test_new_warning_cannot_be_offset_by_a_removed_one(self) -> None:
        observation = self.observation()
        baseline = make_baseline([observation])
        changed = deepcopy(observation)
        changed["diagnostics"]["new"] = changed["diagnostics"].pop("known")
        changed["counts"] = {"new": 1}
        with self.assertRaisesRegex(ValueError, "New/increased"):
            compare(changed, baseline)
        changed = deepcopy(observation)
        changed["counts"]["known"] = 2
        with self.assertRaisesRegex(ValueError, "New/increased"):
            compare(changed, baseline)

    def test_missing_coverage_and_version_change_are_rejected(self) -> None:
        observation = self.observation()
        baseline = make_baseline([observation])
        for field, value in (("files", {"src/a.f90": "a"}), ("compilation_count", 1),
                             ("target_compilation_counts", {"another-target": 2}),
                             ("compiler_version", "15.4.0"), ("case", "absent")):
            changed = deepcopy(observation)
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                compare(changed, baseline)
        with self.assertRaisesRegex(ValueError, "schema"):
            compare(observation, {**baseline, "schema_version": 2})
        changed = deepcopy(observation)
        changed["case"] = "another"
        changed["files"]["src/a.f90"] = "modified"
        with self.assertRaisesRegex(ValueError, "Source changed"):
            make_baseline([observation, changed])

    def test_warning_categories_keep_correctness_separate(self) -> None:
        for option, text in (("-Warray-bounds", "bad index"), ("", "INTENT mismatch in argument"),
                             ("-Wsurprising", "Array larger than limit")):
            self.assertEqual(category(option, text), "correctness-review")
        self.assertEqual(category("", "Fortran 2018 obsolescent feature: Labeled DO"), "obsolescent-standard")
        self.assertEqual(category("-Wunused-variable", "Unused x"), "general")

    @patch("warning_inventory.subprocess.check_output", return_value="15.3.0\n")
    def test_records_must_cover_current_commands_sources_and_objects(self, _version) -> None:
        with tempfile.TemporaryDirectory(prefix="warning test ") as temp:
            source = Path(temp) / "source"
            build = source / "build"
            build.mkdir(parents=True)
            path = source / "unit.f90"
            path.write_text("end\n")
            obj = build / "CMakeFiles/demo.dir/unit.f90.o"
            obj.parent.mkdir(parents=True)
            obj.write_bytes(b"object")
            flags = ["-std=f2018", "-pedantic-errors", "-Wall", "-Wextra", "-Werror=implicit-interface",
                     "-Werror=implicit-procedure", "-fno-fast-math", "-ffp-contract=off", "-O2"]
            command = ["/compiler/gfortran-15", *flags, "-o", str(obj), "-c", str(path)]
            filename, object_name, normalized = command_parts(command, build, source, build)
            self.assertIn("<build>/CMakeFiles/demo.dir/unit.f90.o", normalized)
            write_json(build / "f2018-features.json", {"profile": "f2018-release", "serial": 0,
                       "mpi": 0, "csfg": 0, "openmp": 0})
            (build / "CMakeCache.txt").write_text(f"CMAKE_HOME_DIRECTORY:INTERNAL={source}\n"
                "CMAKE_Fortran_COMPILER:FILEPATH=/compiler/gfortran-15\n")
            commands = [{"directory": str(build), "file": str(path), "arguments": command}]
            (build / "compile_commands.json").write_text(json.dumps(commands))
            record = record_path(build, object_name)
            with self.assertRaisesRegex(ValueError, "Missing warning record"):
                inventory(build)
            write_json(record, {"source": filename, "object": object_name, "command": normalized,
                "source_sha256": digest(path), "object_sha256": digest(obj), "diagnostics": []})
            self.assertEqual(inventory(build)["counts"], {})
            path.write_text("! changed\nend\n")
            with self.assertRaisesRegex(ValueError, "Stale warning record"):
                inventory(build)
            path.write_text("end\n")
            obj.write_bytes(b"foreign object")
            with self.assertRaisesRegex(ValueError, "Stale warning record"):
                inventory(build)
            obj.write_bytes(b"object")
            commands[0]["arguments"].append("-g")
            (build / "compile_commands.json").write_text(json.dumps(commands))
            with self.assertRaisesRegex(ValueError, "Stale warning record"):
                inventory(build)


if __name__ == "__main__":
    unittest.main()
