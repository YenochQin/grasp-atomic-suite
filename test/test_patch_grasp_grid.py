"""Offline patch checks, with a Fortran compile check when gfortran is available."""

from __future__ import annotations

import importlib.util
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "patch_grasp_grid.py"
SPEC = importlib.util.spec_from_file_location("patch_grasp_grid", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
grid = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = grid
SPEC.loader.exec_module(grid)

POINT_AND_FINITE = """      IF (NPARM == 0) THEN
         RNT = EXP((-65.D0/16.D0))/Z
         H = 0.0625D0
         N = MIN(220,NNNP)
      ELSE
         RNT = 2.D-6/Z
         H = 0.05D0
         N = NNNP
      ENDIF
      HP = 0.D0
"""
FINITE = """      RNT = 2.D-6/Z
      H = 0.05D0
      N = NNNP
      HP = 0.D0
"""


def fixture(root: Path, layout: str = "grasp2018") -> None:
    suite = layout == "atomic-suite"
    if suite:
        repository = SCRIPT.parents[1]
        for relative in (
            grid.SUITE_CONFIG_FILE,
            grid.SUITE_DEFAULTS_FILE,
            grid.CAPACITY_FILES[0],
        ):
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            text = (repository / relative).read_text()
            if relative == grid.SUITE_CONFIG_FILE:
                # Fixture values stay independent of the user's selected grid.
                for name, expression in {
                    "NNNP": "590",
                    "FINITE_N": "NNNP",
                    "POINT_N": "MIN(220,NNNP)",
                    "FINITE_H": "0.05D0",
                    "POINT_H": "0.0625D0",
                    "FINITE_RNT_SCALE": "2.0D-6",
                    "POINT_RNT_SCALE": "EXP(-65.0D0/16.0D0)",
                    "DEFAULT_HP": "0.0D0",
                    "DEFAULT_ACCY": "0.0D0",
                }.items():
                    text = re.sub(
                        rf"(\b{name}\s*=\s*)[^!\n]+", lambda m: m[1] + expression, text
                    )
            path.write_text(text)
        for relative in grid.SUITE_GRID_FILES:
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("      CALL SET_RADIAL_DEFAULTS(NPARM, Z)\n")
        return
    for relative in grid.CAPACITY_FILES:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "      integer, parameter :: NNNP = 590\n      integer, parameter :: NNN1 = 600\n"
        )
    for relative in grid.GRID_FILES:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        block = FINITE if "/rci90" in relative else POINT_AND_FINITE
        if relative.endswith("rwfnrelabel.f90"):
            block = "\n".join("!" + line for line in block.splitlines()) + "\n"
        text = block + "      ACCY = H**6\n"
        if "/rmcdhf90" in relative:
            text += "      IF (NDEF /= 0) THEN\n         WRITE (*,*) 'Revise the default ACCY = ', ACCY\n         READ *, ACCY\n      ENDIF\n"
        text += "      CALL SETQIC\n! NNNP = 590 is a historical comment\n"
        if relative.endswith("rwfnrelabel.f90"):
            # Match upstream's host IMPLICIT NONE and contained GETHFD scope.
            text = (
                "      PROGRAM RWFNRELABEL\n"
                "      USE parameter_def\n"
                "      IMPLICIT NONE\n"
                "      CONTAINS\n"
                "      SUBROUTINE GETHFD(NAME)\n"
                "      IMPLICIT DOUBLEPRECISION (A-H,O-Z)\n"
                "      INTEGER N\n"
                "      CHARACTER*24 NAME\n"
                + text
                + "      END SUBROUTINE\n      END PROGRAM\n"
            )
        path.write_text(text)
    for relative in grid.RESTART_FILES:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("      READ (IMCDF) RNT, H, HP\n      ACCY = H**6\n")


class PatchGridTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "grasp"
        fixture(self.root)
        self.script_directory = Path(self.temporary.name) / "scripts"
        self.script_directory.mkdir()
        script_location = patch.object(
            grid, "__file__", str(self.script_directory / SCRIPT.name)
        )
        script_location.start()
        self.addCleanup(script_location.stop)

    def write_config(self, directory: Path | None = None, **overrides: object) -> Path:
        config = {
            "GRASP_SOURCE": str(self.root),
            "SOURCE_LAYOUT": "grasp2018",
            "RUN_MODE": "preview",
            "PRINT_DIFF": False,
            "GRID_PARAMETERS": {
                "nnnp": 3990,
                "n": 2946,
                "h": 0.010,
                "rnt_scale": 2e-6,
                "hp": 0.0,
                "accy": 1e-10,
            },
            **overrides,
        }
        lines = []
        for name, value in config.items():
            if value is not None and not isinstance(value, dict):
                lines.append(f"{name} = {json.dumps(value)}\n")
        for name, value in config.items():
            if isinstance(value, dict):
                lines.append(f"\n[{name}]\n")
                lines.extend(
                    f"{key} = {json.dumps(item)}\n"
                    for key, item in value.items()
                    if item is not None
                )
        path = (directory or self.script_directory) / "config.toml"
        path.write_text("".join(lines))
        return path

    def test_toml_preview_apply_check_and_restore(self) -> None:
        originals = {p: p.read_bytes() for p in self.root.rglob("*.f90")}
        backup = Path(self.temporary.name) / "config-backup"
        parameters = {
            "nnnp": 1990,
            "n": 1179,
            "h": 0.025,
            "rnt_scale": 2e-6,
            "hp": 0.0,
            "accy": None,
        }
        config = {"GRID_PARAMETERS": parameters, "BACKUP_DIRECTORY": str(backup)}
        self.write_config(**config)
        with redirect_stdout(io.StringIO()):
            self.assertEqual(grid.main([]), 0)
            self.assertFalse(backup.exists())
            for path, data in originals.items():
                self.assertEqual(path.read_bytes(), data)
            self.write_config(**config, RUN_MODE="check")
            self.assertEqual(grid.main([]), 1)
            self.write_config(**config, RUN_MODE="apply")
            self.assertEqual(grid.main([]), 0)
            # The already-existing backup must not break repeated application.
            self.assertEqual(grid.main([]), 0)
            self.assertTrue((backup / "manifest.json").is_file())
            self.assertEqual(grid.plan(self.root, grid.Settings(**parameters)), [])
            self.write_config(**config, RUN_MODE="check")
            self.assertEqual(grid.main([]), 0)
            self.write_config(**config, RESTORE_DIRECTORY=str(backup))
            self.assertEqual(grid.main([]), 0)  # Preview restoration.
            self.assertEqual(grid.plan(self.root, grid.Settings(**parameters)), [])
            self.write_config(
                **config,
                GRASP_SOURCE=None,
                RESTORE_DIRECTORY=str(backup),
                RUN_MODE="apply",
            )
            self.assertEqual(grid.main([]), 0)
            for path, data in originals.items():
                self.assertEqual(path.read_bytes(), data)

    def test_copied_script_relative_root_without_automatic_backup(self) -> None:
        root = Path(self.temporary.name) / "copied-suite"
        fixture(root, "atomic-suite")
        self.write_config(
            root,
            GRASP_SOURCE=".",
            SOURCE_LAYOUT="atomic-suite",
            RUN_MODE="apply",
            GRID_PARAMETERS={"nnnp": 1990, "h": 0.025, "accy": 0},
        )
        with (
            patch.object(grid, "__file__", str(root / SCRIPT.name)),
            patch.object(grid.tempfile, "tempdir", str(self.temporary.name)),
            redirect_stdout(io.StringIO()),
        ):
            files_before = set(Path(self.temporary.name).rglob("*"))
            self.assertEqual(grid.main([]), 0)
            self.assertEqual(set(Path(self.temporary.name).rglob("*")), files_before)
            self.assertEqual(
                grid.plan(
                    root, grid.Settings(nnnp=1990, h=0.025, accy=0.0), "atomic-suite"
                ),
                [],
            )

    def test_cli_apply_without_backup_can_switch_grid_parameters(self) -> None:
        files_before = set(Path(self.temporary.name).rglob("*"))
        with (
            patch.object(grid.tempfile, "tempdir", str(self.temporary.name)),
            redirect_stdout(io.StringIO()),
        ):
            for capacity, count, step in ((1990, 1179, 0.025), (590, 590, 0.05)):
                self.assertEqual(
                    grid.main(
                        [
                            "--grasp",
                            str(self.root),
                            "--nnnp",
                            str(capacity),
                            "--n",
                            str(count),
                            "--h",
                            str(step),
                            "--apply",
                        ]
                    ),
                    0,
                )
                self.assertEqual(
                    grid.plan(self.root, grid.Settings(nnnp=capacity, n=count, h=step)),
                    [],
                )
                self.assertEqual(
                    set(Path(self.temporary.name).rglob("*")), files_before
                )

    def test_invalid_toml_configuration_does_not_write_source(self) -> None:
        original = (self.root / grid.CAPACITY_FILES[0]).read_bytes()
        with (
            redirect_stderr(io.StringIO()),
            redirect_stdout(io.StringIO()),
        ):
            invalid = (
                {"GRASP_SOURCE": None},
                {"RUN_MODE": "unknown"},
                {"SOURCE_LAYOUT": "unknown"},
                {"GRID_PARAMETERS": {"nnnp": 1990, "n": 2000}},
                {"GRID_PARAMETERS": {"nnnp": 1990.5}},
                {"GRID_PARAMETERS": {"typo": 1990}},
                {"GRID_PARAMETERS": {"nnnp": None}},
                {"PRINT_DIFF": "False"},
                {"GRID_PARAMETERS": {"nnnp": True}},
                {"GRID_PARAMETERS": {"h": "0.025"}},
                {"GRID_PARAMETERS": "not a table"},
                {"GRID_PARAMETERS": {"point_n": 220.0}},
                {"GRASP_SOURCE": ""},
                {"BACKUP_DIRECTORY": False},
                {"RESTORE_DIRECTORY": ""},
                {"TYPO": 1},
            )
            for values in invalid:
                with self.subTest(values=values):
                    self.write_config(
                        **{
                            "RUN_MODE": "apply",
                            "GRID_PARAMETERS": {"nnnp": 1990},
                            **values,
                        }
                    )
                    with self.assertRaises(SystemExit) as error:
                        grid.main([])
                    self.assertEqual(error.exception.code, 2)
                    self.assertEqual(
                        (self.root / grid.CAPACITY_FILES[0]).read_bytes(), original
                    )
            self.assertFalse((self.root / "grid-backups").exists())

    def test_missing_or_malformed_toml_fails_without_writes(self) -> None:
        originals = {p: p.read_bytes() for p in self.root.rglob("*.f90")}
        config = self.script_directory / "config.toml"
        with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
            for content in (
                None,
                'RUN_MODE = "apply"\n[GRID_PARAMETERS\n',
                "PRINT_DIFF = True\n",
            ):
                with self.subTest(content=content):
                    if content is not None:
                        config.write_text(content)
                    with self.assertRaises(SystemExit) as error:
                        grid.main([])
                    self.assertEqual(error.exception.code, 2)
                    for path, original in originals.items():
                        self.assertEqual(path.read_bytes(), original)

    def test_copied_script_loads_config_from_another_working_directory(self) -> None:
        copied = self.root / SCRIPT.name
        shutil.copyfile(SCRIPT, copied)
        shutil.copyfile(
            SCRIPT.with_name("config.example.toml"), self.root / "config.toml"
        )
        original = (self.root / grid.CAPACITY_FILES[0]).read_bytes()
        result = subprocess.run(
            [sys.executable, str(copied)],
            cwd=self.script_directory,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"Source: {self.root.resolve()}", result.stdout)
        self.assertIn(
            f"Configuration: {self.root.resolve() / 'config.toml'}", result.stdout
        )
        self.assertEqual((self.root / grid.CAPACITY_FILES[0]).read_bytes(), original)
        self.assertIn("NNNP = 3990", result.stdout)
        self.assertIn("N = 2946", result.stdout)
        self.assertIn("H = 1D-2", result.stdout)
        self.assertIn("ACCY = 1D-10", result.stdout)

    def test_toml_defaults_and_omitted_parameters_preserve_values(self) -> None:
        (self.script_directory / "config.toml").write_text(
            'GRASP_SOURCE = "../grasp"\n[GRID_PARAMETERS]\nh = 0.025\n'
        )
        original = (self.root / grid.CAPACITY_FILES[0]).read_bytes()
        with (
            patch.object(grid, "plan", wraps=grid.plan) as planner,
            redirect_stdout(io.StringIO()),
        ):
            self.assertEqual(grid.main([]), 0)
            self.assertEqual(
                planner.call_args.args,
                (self.root.resolve(), grid.Settings(h=0.025), "grasp2018"),
            )
        self.assertEqual((self.root / grid.CAPACITY_FILES[0]).read_bytes(), original)
        self.assertFalse((self.root / "grid-backups").exists())

    def test_explicit_cli_keeps_its_values_independent_of_toml(self) -> None:
        config = self.script_directory / "config.toml"
        with (
            patch.object(grid, "plan", wraps=grid.plan) as planner,
            redirect_stdout(io.StringIO()),
        ):
            # Explicit CLI works both without a config and with an invalid one.
            for content in (None, "invalid TOML"):
                if content is not None:
                    config.write_text(content)
                self.assertEqual(
                    grid.main(["--grasp", str(self.root), "--h", "0.025"]), 0
                )
                self.assertEqual(planner.call_args.args[1], grid.Settings(h=0.025))

    def test_suite_layout_patches_all_present_programs_and_restores(self) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        originals = {p: p.read_bytes() for p in root.rglob("*.f90")}
        capacity_only = grid.plan(root, grid.Settings(nnnp=2990), "atomic-suite")
        self.assertEqual([c.relative for c in capacity_only], [grid.SUITE_CONFIG_FILE])
        self.assertIn(b"FINITE_N = NNNP", capacity_only[0].after)
        self.assertIn(b"NNN1 = NNNP + 10", capacity_only[0].after)
        settings = grid.Settings(nnnp=1990, n=1179, h=0.025)
        changes = grid.plan(root, settings, "atomic-suite")
        self.assertEqual(len(changes), 1)
        self.assertIn(b"FINITE_N = 1179", changes[0].after)
        self.assertIn(b"FINITE_H = 2.5D-2", changes[0].after)
        backup = grid.apply(root, changes, Path(self.temporary.name) / "suite-backup")
        self.assertEqual(grid.plan(root, settings, "atomic-suite"), [])
        restored_root, restoration = grid.restore_plan(backup, root)
        grid.write_changes(restored_root, restoration)
        for path, data in originals.items():
            self.assertEqual(path.read_bytes(), data)

    def test_partial_layout_still_requires_every_supported_program(self) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        with self.assertRaisesRegex(ValueError, "missing required source"):
            grid.plan(root, grid.Settings(nnnp=1990))
        (root / grid.SUITE_GRID_FILES[-1]).unlink()
        with self.assertRaisesRegex(ValueError, "missing required source"):
            grid.plan(root, grid.Settings(nnnp=1990), "atomic-suite")

    def test_suite_unknown_initializer_blocks_changes(self) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        unknown = root / "src/appl/new_program.f90"
        unknown.write_text("      RNT = 2.D-6/Z\n")
        with self.assertRaisesRegex(ValueError, "unrecognized grid initializer"):
            grid.plan(root, grid.Settings(nnnp=1990), "atomic-suite")

    def test_suite_duplicate_capacity_and_literal_derived_limit_are_rejected(
        self,
    ) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        duplicate = root / "src/duplicate.f90"
        duplicate.write_text("      integer, parameter :: NNNP = 590\n")
        with self.assertRaisesRegex(ValueError, "duplicate capacity"):
            grid.plan(root, grid.Settings(h=0.025), "atomic-suite")
        duplicate.unlink()
        config = root / grid.SUITE_CONFIG_FILE
        config.write_text(config.read_text().replace("NNN1 = NNNP + 10", "NNN1 = 600"))
        with self.assertRaisesRegex(ValueError, "must remain derived"):
            grid.plan(root, grid.Settings(h=0.025), "atomic-suite")

    def test_suite_local_math_step_is_not_a_radial_grid_default(self) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        local = root / "src/local_math.f90"
        local.write_text(
            "subroutine fit\nreal*8 :: h\nh = 0.01d0\nend subroutine fit\n"
        )
        changes = grid.plan(root, grid.Settings(h=0.025), "atomic-suite")
        self.assertEqual([c.relative for c in changes], [grid.SUITE_CONFIG_FILE])
        local.write_text("use grid_C, only: h\nh = 0.01d0\n")
        with self.assertRaisesRegex(ValueError, "unrecognized grid initializer"):
            grid.plan(root, grid.Settings(h=0.025), "atomic-suite")

    def test_suite_fixed_accuracy_can_return_to_adaptive(self) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        fixed = grid.plan(root, grid.Settings(accy=1e-10), "atomic-suite")
        grid.write_changes(root, fixed)
        self.assertIn(b"DEFAULT_ACCY = 1D-10", fixed[0].after)
        reset = grid.plan(root, grid.Settings(accy=0), "atomic-suite")
        grid.write_changes(root, reset)
        self.assertEqual(grid.plan(root, grid.Settings(accy=0), "atomic-suite"), [])

    def test_suite_capacity_reduction_rejects_existing_runtime_point_count(
        self,
    ) -> None:
        root = Path(self.temporary.name) / "suite"
        fixture(root, "atomic-suite")
        grid.write_changes(
            root, grid.plan(root, grid.Settings(nnnp=1990, n=1179), "atomic-suite")
        )
        with self.assertRaisesRegex(ValueError, "n must be between"):
            grid.plan(root, grid.Settings(nnnp=590), "atomic-suite")

    def test_finite_patch_preserves_point_model_and_comments(self) -> None:
        changes = grid.plan(
            self.root, grid.Settings(nnnp=1990, n=1179, h=0.025, rnt_scale=2e-6)
        )
        output = {change.relative: change.after.decode() for change in changes}
        for relative in grid.GRID_FILES:
            self.assertIn("N = 1179", output[relative])
            self.assertIn("H = 2.5D-2", output[relative])
            self.assertIn("! NNNP = 590 is a historical comment", output[relative])
            if "/rci90" not in relative:
                self.assertIn("N = MIN(220,NNNP)", output[relative])
                self.assertIn("H = 0.0625D0", output[relative])
        self.assertTrue(
            all("NNNP = 1990" in output[path] for path in grid.CAPACITY_FILES)
        )
        # Planning is read-only, including on the previously broken helper.
        self.assertIn("!      IF", (self.root / "src/tool/rwfnrelabel.f90").read_text())

    def test_apply_is_idempotent_and_restore_is_byte_exact(self) -> None:
        before = {
            path: (self.root / path).read_bytes()
            for path in (*grid.GRID_FILES, *grid.CAPACITY_FILES, *grid.RESTART_FILES)
        }
        settings = grid.Settings(nnnp=1990, n=1179, h=0.025, accy=1e-10)
        changes = grid.plan(self.root, settings)
        backup = grid.apply(self.root, changes, Path(self.temporary.name) / "backup")
        self.assertEqual(grid.plan(self.root, settings), [])
        # Changing only the grid must preserve a previously selected fixed ACCY.
        reconfigured = grid.plan(self.root, grid.Settings(h=0.02))
        rmcdhf_text = next(
            change.after.decode()
            for change in reconfigured
            if change.relative == "src/appl/rmcdhf90/getscd.f90"
        )
        self.assertEqual(rmcdhf_text.count("ACCY = 1D-10"), 2)
        root, restoration = grid.restore_plan(backup, self.root)
        grid.write_changes(root, restoration)
        for relative, data in before.items():
            self.assertEqual((root / relative).read_bytes(), data)

    def test_relabel_imports_shared_nuclear_model_and_repairs_old_patch(self) -> None:
        relative = "src/tool/rwfnrelabel.f90"
        path = self.root / relative
        settings = grid.Settings(
            nnnp=2990, n=1179, h=0.025, rnt_scale=2e-6, hp=0.0, accy=1e-10
        )
        changes = grid.plan(self.root, settings)
        output = next(c.after.decode() for c in changes if c.relative == relative)
        preamble = output.split("SUBROUTINE GETHFD(NAME)", 1)[1].split("IMPLICIT", 1)[0]
        self.assertIn("USE npar_C, ONLY: NPARM", preamble)
        grid.write_changes(self.root, changes)
        self.assertEqual(grid.plan(self.root, settings), [])
        # The previous patcher already activated the block but omitted the USE.
        path.write_text(output.replace("      USE npar_C, ONLY: NPARM\n", ""))
        repair = grid.plan(self.root, settings)
        self.assertEqual([c.relative for c in repair], [relative])
        self.assertEqual(repair[0].after.decode(), output)

    def test_relabel_keeps_existing_nuclear_model_import(self) -> None:
        path = self.root / "src/tool/rwfnrelabel.f90"
        path.write_text(
            path.read_text().replace(
                "      IMPLICIT DOUBLEPRECISION",
                "      USE npar_C\n      IMPLICIT DOUBLEPRECISION",
            )
        )
        changes = grid.plan(self.root, grid.Settings(h=0.025))
        output = next(
            c.after.decode() for c in changes if c.relative.endswith("rwfnrelabel.f90")
        )
        self.assertEqual(output.lower().count("use npar_c"), 1)

    @unittest.skipUnless(shutil.which("gfortran"), "gfortran is not available")
    def test_d1_relabel_compiles_with_implicit_none_host(self) -> None:
        changes = grid.plan(
            self.root,
            grid.Settings(
                nnnp=2990, n=1179, h=0.025, rnt_scale=2e-6, hp=0.0, accy=1e-10
            ),
        )
        output = next(
            c.after for c in changes if c.relative.endswith("rwfnrelabel.f90")
        )
        build = Path(self.temporary.name) / "compile"
        build.mkdir()
        # Real module ownership: NPARM is shared state, not a local declaration.
        (build / "modules.f90").write_text(
            "module parameter_def\ninteger, parameter :: NNNP = 2990\nend module\n"
            "module npar_C\ninteger :: NPARM\nend module\n"
        )
        (build / "rwfnrelabel.f90").write_bytes(output)
        result = subprocess.run(
            ["gfortran", "-c", "modules.f90", "rwfnrelabel.f90"],
            cwd=build,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_unknown_capacity_blocks_all_writes(self) -> None:
        path = self.root / "src/lib/new_grid.f90"
        path.write_text("      integer, parameter :: NNNP = 590\n")
        original = (self.root / grid.CAPACITY_FILES[0]).read_bytes()
        with self.assertRaisesRegex(ValueError, "unmatched capacity"):
            grid.plan(self.root, grid.Settings(nnnp=1990))
        self.assertEqual((self.root / grid.CAPACITY_FILES[0]).read_bytes(), original)

    def test_unknown_initializer_and_malformed_branch_are_rejected(self) -> None:
        path = self.root / "src/lib/new_grid.f90"
        path.write_text("      RNT = 2e-6/Z\n")
        with self.assertRaisesRegex(ValueError, "unrecognized grid initializer"):
            grid.plan(self.root, grid.Settings(h=0.025))
        path.write_text("! no active assignments\n")
        current = self.root / grid.GRID_FILES[0]
        current.write_text(current.read_text().replace("H = 0.05D0", "HH = 0.05D0"))
        with self.assertRaisesRegex(ValueError, "expected one H assignment"):
            grid.plan(self.root, grid.Settings(h=0.025))

    def test_invalid_numbers_and_runtime_capacity_conflict(self) -> None:
        for settings in (
            grid.Settings(nnnp=10),
            grid.Settings(nnnp=590, n=591),
            grid.Settings(h=float("nan")),
            grid.Settings(h=-0.01),
            grid.Settings(hp=-1),
            grid.Settings(accy=0),
        ):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                grid.plan(self.root, settings)
        path = self.root / grid.GRID_FILES[0]
        path.write_text(path.read_text().replace("N = NNNP", "N = 1990"))
        with self.assertRaisesRegex(ValueError, "incompatible"):
            grid.plan(self.root, grid.Settings(nnnp=590))

    def test_point_model_and_common_hp_can_be_selected_explicitly(self) -> None:
        changes = grid.plan(
            self.root,
            grid.Settings(
                nnnp=1990, point_n=400, point_h=0.03, point_rnt_scale=0.01, hp=0.1
            ),
        )
        text = next(
            change.after.decode()
            for change in changes
            if change.relative == grid.GRID_FILES[0]
        )
        self.assertIn("N = 400", text)
        self.assertIn("H = 0.05D0", text)  # finite model unchanged
        self.assertIn("HP = 1D-1", text)

    def test_restore_refuses_subsequent_source_edits(self) -> None:
        changes = grid.plan(self.root, grid.Settings(nnnp=1990))
        backup = grid.apply(self.root, changes, Path(self.temporary.name) / "backup")
        path = self.root / changes[0].relative
        path.write_bytes(path.read_bytes() + b"! user's later edit\n")
        with self.assertRaisesRegex(ValueError, "subsequent edits"):
            grid.restore_plan(backup, self.root)
        self.assertTrue(path.read_bytes().endswith(b"! user's later edit\n"))

    def test_restore_refuses_corrupted_backup(self) -> None:
        changes = grid.plan(self.root, grid.Settings(nnnp=1990))
        backup = grid.apply(self.root, changes, Path(self.temporary.name) / "backup")
        saved = backup / changes[0].relative
        saved.write_bytes(saved.read_bytes() + b"! corrupted backup\n")
        with self.assertRaisesRegex(ValueError, "backup checksum mismatch"):
            grid.restore_plan(backup, self.root)

    def test_write_failure_rolls_back_completed_files(self) -> None:
        changes = grid.plan(self.root, grid.Settings(nnnp=1990))
        original_write = grid.atomic_write
        calls = 0

        def failing_write(path: Path, data: bytes) -> None:
            nonlocal calls
            calls += 1
            if calls == 2:
                raise OSError("simulated disk failure")
            original_write(path, data)

        with patch.object(grid, "atomic_write", side_effect=failing_write):
            with self.assertRaisesRegex(OSError, "simulated disk failure"):
                grid.write_changes(self.root, changes)
        for change in changes:
            self.assertEqual((self.root / change.relative).read_bytes(), change.before)


if __name__ == "__main__":
    unittest.main()
