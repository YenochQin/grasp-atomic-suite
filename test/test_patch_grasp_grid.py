"""Offline patch checks, with a Fortran compile check when gfortran is available."""

from __future__ import annotations

import importlib.util
from contextlib import redirect_stderr, redirect_stdout
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "patch_grasp_grid.py"
FORTRAN_COMPILER = os.environ.get("FC", "gfortran")
SPEC = importlib.util.spec_from_file_location("patch_grasp_grid", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
grid = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = grid
SPEC.loader.exec_module(grid)

POINT_AND_FINITE = """      IF (nuclear_model_parameter_count == 0) THEN
         radial_grid_scale = EXP((-65.D0/16.D0))/nuclear_charge
         transformed_grid_step = 0.0625D0
         radial_point_count = MIN(220,max_radial_points)
      ELSE
         radial_grid_scale = 2.D-6/nuclear_charge
         transformed_grid_step = 0.05D0
         radial_point_count = max_radial_points
      ENDIF
      asymptotic_grid_step = 0.D0
"""
FINITE = """      radial_grid_scale = 2.D-6/nuclear_charge
      transformed_grid_step = 0.05D0
      radial_point_count = max_radial_points
      asymptotic_grid_step = 0.D0
"""


def fixture(root: Path) -> None:
    for relative in (grid.SUITE_CONFIG_FILE, grid.SUITE_EXPORT_FILE):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(SCRIPT.parents[1] / relative, path)
    for relative in grid.SUITE_GRID_FILES:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        block = FINITE if "/rci" in relative else POINT_AND_FINITE
        text = block + "      working_tolerance = transformed_grid_step**6\n"
        if "/rmcdhf" in relative:
            text += (
                "      IF (NDEF /= 0) THEN\n"
                "         WRITE (*,*) 'Revise the default ACCY = ', working_tolerance\n"
                "         READ *, working_tolerance\n"
                "      ENDIF\n"
            )
        path.write_text(text + "! max_radial_points = 590 is a historical comment\n")
    for relative in grid.SUITE_RESTART_FILES:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("      working_tolerance = transformed_grid_step**6\n")


class PatchGridTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name) / "suite"
        fixture(self.root)
        self.script_directory = self.root / "scripts"
        self.script_directory.mkdir()
        script_location = patch.object(
            grid, "__file__", str(self.script_directory / SCRIPT.name)
        )
        script_location.start()
        self.addCleanup(script_location.stop)

    def write_config(self, directory: Path | None = None, **overrides: object) -> Path:
        config = {
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
                RESTORE_DIRECTORY=str(backup),
                RUN_MODE="apply",
            )
            self.assertEqual(grid.main([]), 0)
            for path, data in originals.items():
                self.assertEqual(path.read_bytes(), data)

    def test_copied_script_patches_its_own_source_without_automatic_backup(self) -> None:
        root = Path(self.temporary.name) / "copied-suite"
        fixture(root)
        scripts = root / "scripts"
        scripts.mkdir()
        self.write_config(
            scripts,
            RUN_MODE="apply",
            GRID_PARAMETERS={"nnnp": 1990, "h": 0.025, "accy": 0},
        )
        with (
            patch.object(grid, "__file__", str(scripts / SCRIPT.name)),
            patch.object(grid.tempfile, "tempdir", str(self.temporary.name)),
            redirect_stdout(io.StringIO()),
        ):
            files_before = set(Path(self.temporary.name).rglob("*"))
            self.assertEqual(grid.main([]), 0)
            self.assertEqual(set(Path(self.temporary.name).rglob("*")), files_before)
            self.assertEqual(
                grid.plan(root, grid.Settings(nnnp=1990, h=0.025, accy=0.0)),
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
        original = (self.root / grid.SUITE_CONFIG_FILE).read_bytes()
        with (
            redirect_stderr(io.StringIO()),
            redirect_stdout(io.StringIO()),
        ):
            invalid = (
                {"RUN_MODE": "unknown"},
                {"SOURCE_LAYOUT": "atomic-suite"},
                {"GRID_PARAMETERS": {"nnnp": 1990, "n": 2000}},
                {"GRID_PARAMETERS": {"nnnp": 1990.5}},
                {"GRID_PARAMETERS": {"typo": 1990}},
                {"GRID_PARAMETERS": {"nnnp": None}},
                {"PRINT_DIFF": "False"},
                {"CAPACITY_ONLY": "true"},
                {"GRID_PARAMETERS": {"nnnp": True}},
                {"GRID_PARAMETERS": {"h": "0.025"}},
                {"GRID_PARAMETERS": "not a table"},
                {"GRID_PARAMETERS": {"point_n": 220.0}},
                {"GRASP_SOURCE": ".."},
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
                        (self.root / grid.SUITE_CONFIG_FILE).read_bytes(), original
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
        copied = self.script_directory / SCRIPT.name
        shutil.copyfile(SCRIPT, copied)
        shutil.copyfile(
            SCRIPT.with_name("config.example.toml"), copied.with_name("config.toml")
        )
        original = (self.root / grid.SUITE_CONFIG_FILE).read_bytes()
        other_root = Path(self.temporary.name) / "other-suite"
        fixture(other_root)
        other_capacity = other_root / grid.SUITE_CONFIG_FILE
        script_link = other_root / SCRIPT.name
        script_link.symlink_to(copied)
        result = subprocess.run(
            [sys.executable, str(script_link)],
            cwd=other_root,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f"Source: {self.root.resolve()}", result.stdout)
        self.assertIn(
            f"Configuration: {copied.resolve().with_name('config.toml')}", result.stdout
        )
        self.assertEqual((self.root / grid.SUITE_CONFIG_FILE).read_bytes(), original)
        self.assertIn("max_radial_points = 3990", result.stdout)
        self.assertIn("radial_point_count = 2946", result.stdout)
        self.assertIn("transformed_grid_step = 1d-2", result.stdout)
        self.assertIn("working_tolerance = 1d-10", result.stdout)
        result = subprocess.run(
            [sys.executable, str(script_link), "--nnnp", "1990", "--apply"],
            cwd=other_root,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(
            b"max_radial_points = 1990", (self.root / grid.SUITE_CONFIG_FILE).read_bytes()
        )
        self.assertEqual(other_capacity.read_bytes(), original)

    def test_toml_defaults_and_omitted_parameters_preserve_values(self) -> None:
        (self.script_directory / "config.toml").write_text(
            "[GRID_PARAMETERS]\nh = 0.025\n"
        )
        original = (self.root / grid.SUITE_CONFIG_FILE).read_bytes()
        with (
            patch.object(grid, "plan", wraps=grid.plan) as planner,
            redirect_stdout(io.StringIO()),
        ):
            self.assertEqual(grid.main([]), 0)
            self.assertEqual(
                planner.call_args.args,
                (self.root.resolve(), grid.Settings(h=0.025)),
            )
        self.assertEqual((self.root / grid.SUITE_CONFIG_FILE).read_bytes(), original)
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
                    grid.main(["--h", "0.025"]), 0
                )
                self.assertEqual(planner.call_args.args[1], grid.Settings(h=0.025))

    def test_cli_and_toml_capacity_options_change_only_max_radial_points(self) -> None:
        originals = {p: p.read_bytes() for p in self.root.rglob("*.f90")}
        config = self.root / grid.SUITE_CONFIG_FILE
        self.write_config(RUN_MODE="apply", GRID_PARAMETERS={"nnnp": 2990})
        with redirect_stdout(io.StringIO()):
            self.assertEqual(grid.main(["--nnnp", "1990", "--apply"]), 0)
            self.assertEqual(
                {p: p.read_bytes() for p in originals},
                originals | {config: originals[config].replace(b"= 590", b"= 1990")},
            )
            self.assertEqual(grid.main([]), 0)
        self.assertEqual(
            {p: p.read_bytes() for p in originals},
            originals | {config: originals[config].replace(b"= 590", b"= 2990")},
        )

    def test_cli_rejects_removed_source_options(self) -> None:
        for option, value in (("--grasp", str(self.root)), ("--layout", "atomic-suite")):
            with self.subTest(option=option), redirect_stderr(io.StringIO()) as stderr:
                with self.assertRaises(SystemExit) as error:
                    grid.main([option, value, "--nnnp", "1990", "--apply"])
                self.assertEqual(error.exception.code, 2)
                self.assertIn(f"unrecognized arguments: {option}", stderr.getvalue())

    def test_legacy_layout_is_rejected_for_patching_and_restoration(self) -> None:
        root = Path(self.temporary.name) / "legacy"
        relative = "src/lib/libmod/parameter_def_M.f90"
        capacity_file = root / relative
        capacity_file.parent.mkdir(parents=True)
        before = b"      integer, parameter :: NNNP = 590\n"
        after = before.replace(b"590", b"1990")
        capacity_file.write_bytes(after)
        backup = Path(self.temporary.name) / "legacy-backup"
        saved = backup / relative
        saved.parent.mkdir(parents=True)
        saved.write_bytes(before)
        (backup / "manifest.json").write_text(json.dumps({
            "root": str(root),
            "files": [{
                "path": relative, "before": grid.digest(before), "after": grid.digest(after)
            }],
        }))
        for arguments in (
            ["--nnnp", "590", "--apply"],
            ["--restore", str(backup), "--apply"],
        ):
            with (
                self.subTest(arguments=arguments),
                patch.object(grid, "__file__", str(root / "scripts" / SCRIPT.name)),
                redirect_stderr(io.StringIO()),
            ):
                with self.assertRaises(SystemExit) as error:
                    grid.main(arguments)
                self.assertEqual(error.exception.code, 2)
                self.assertEqual(capacity_file.read_bytes(), after)

    def test_restore_rejects_backup_from_another_suite_copy(self) -> None:
        other_root = Path(self.temporary.name) / "other-suite"
        fixture(other_root)
        settings = grid.Settings(nnnp=1990)
        changes = grid.plan(other_root, settings)
        backup = Path(self.temporary.name) / "other-backup"
        with redirect_stdout(io.StringIO()):
            grid.apply(other_root, changes, backup)
        grid.write_changes(self.root, grid.plan(self.root, settings))
        with redirect_stderr(io.StringIO()) as stderr:
            with self.assertRaises(SystemExit) as error:
                grid.main(["--restore", str(backup), "--apply"])
        self.assertEqual(error.exception.code, 2)
        self.assertIn("backup belongs to a different source tree", stderr.getvalue())
        for root in (self.root, other_root):
            self.assertEqual((root / grid.SUITE_CONFIG_FILE).read_bytes(), changes[0].after)

    def test_capacity_only_preserves_every_other_source_line(self) -> None:
        originals = {
            p.relative_to(self.root).as_posix(): p.read_bytes()
            for p in self.root.rglob("*.f90")
        }
        expected = originals.copy()
        expected[grid.SUITE_CONFIG_FILE] = expected[grid.SUITE_CONFIG_FILE].replace(
            b"max_radial_points = 590", b"max_radial_points = 1990"
        )
        for relative in grid.SUITE_GRID_FILES:
            expected[relative] = expected[relative].replace(
                b"radial_point_count = max_radial_points", b"radial_point_count = 590"
            )
        settings = grid.Settings(nnnp=1990, n=590)
        changes = grid.plan(self.root, settings, capacity_only=True)
        self.assertEqual(
            {c.relative: c.after for c in changes},
            {p: data for p, data in expected.items() if data != originals[p]},
        )
        backup = grid.apply(
            self.root, changes, Path(self.temporary.name) / "capacity-only-backup"
        )
        self.assertEqual(grid.plan(self.root, settings, capacity_only=True), [])
        self.assertEqual(
            {p: (self.root / p).read_bytes() for p in originals}, expected
        )
        restoration = grid.restore_plan(backup, self.root)
        grid.write_changes(self.root, restoration)
        self.assertEqual(
            {p: (self.root / p).read_bytes() for p in originals}, originals
        )

    def test_capacity_only_toml_and_cli_apply_and_check(self) -> None:
        parameters = {"nnnp": 1990, "n": 590}
        originals = {p: p.read_bytes() for p in self.root.rglob("*.f90")}
        backup = Path(self.temporary.name) / "capacity-only-backup"
        config = {
            "CAPACITY_ONLY": True,
            "GRID_PARAMETERS": parameters,
            "BACKUP_DIRECTORY": str(backup),
        }
        with redirect_stdout(io.StringIO()):
            self.write_config(**config)
            self.assertEqual(grid.main([]), 0)
            self.assertEqual({p: p.read_bytes() for p in originals}, originals)
            self.write_config(**config, RUN_MODE="check")
            self.assertEqual(grid.main([]), 1)
            self.write_config(**config, RUN_MODE="apply")
            self.assertEqual(grid.main([]), 0)
            self.assertEqual(grid.main([]), 0)
            self.write_config(**config, RUN_MODE="check")
            self.assertEqual(grid.main([]), 0)
            self.assertEqual(
                grid.main(
                    [
                        "--capacity-only",
                        "--nnnp",
                        "1990",
                        "--n",
                        "590",
                        "--apply",
                    ]
                ),
                0,
            )
            # Restoring from TOML ignores patch settings, including this mode.
            self.write_config(**config, RUN_MODE="apply", RESTORE_DIRECTORY=str(backup))
            self.assertEqual(grid.main([]), 0)
            self.assertEqual({p: p.read_bytes() for p in originals}, originals)

    def test_capacity_only_rejects_missing_or_extra_parameters_without_writes(
        self,
    ) -> None:
        originals = {p: p.read_bytes() for p in self.root.rglob("*.f90")}
        invalid = [{"nnnp": 1990}, {"n": 590}]
        invalid.extend(
            {"nnnp": 1990, "n": 590, name: 1}
            for name in grid.Settings.__dataclass_fields__
            if name not in ("nnnp", "n")
        )
        with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
            for parameters in invalid:
                with self.subTest(parameters=parameters):
                    self.write_config(
                        CAPACITY_ONLY=True, GRID_PARAMETERS=parameters, RUN_MODE="apply"
                    )
                    with self.assertRaises(SystemExit) as error:
                        grid.main([])
                    self.assertEqual(error.exception.code, 2)
                    self.assertEqual({p: p.read_bytes() for p in originals}, originals)

    def test_capacity_only_preserves_preexisting_accy_patch(self) -> None:
        settings = grid.Settings(nnnp=1990, n=590, accy=1e-10)
        grid.write_changes(self.root, grid.plan(self.root, settings))
        self.assertEqual(
            grid.plan(self.root, grid.Settings(nnnp=1990, n=590), capacity_only=True),
            [],
        )

    def test_suite_layout_patches_all_present_programs_and_restores(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        originals = {p: p.read_bytes() for p in root.rglob("*.f90")}
        capacity_only = grid.plan(root, grid.Settings(nnnp=2990))
        self.assertEqual([c.relative for c in capacity_only], [grid.SUITE_CONFIG_FILE])
        self.assertIn(b"max_radial_points = 2990", capacity_only[0].after)
        self.assertIn(b"radial_workspace_capacity = max_radial_points + 10", capacity_only[0].after)
        settings = grid.Settings(nnnp=1990, n=1179, h=0.025)
        changes = grid.plan(root, settings)
        self.assertEqual(len(changes), 1 + len(grid.SUITE_GRID_FILES))
        self.assertNotIn(grid.SUITE_EXPORT_FILE, {c.relative for c in changes})
        output = {c.relative: c.after for c in changes}
        for entry in grid.SUITE_GRID_FILES:
            self.assertIn(b"radial_point_count = 1179", output[entry])
            self.assertIn(b"transformed_grid_step = 2.5d-2", output[entry])
            if "/rci" not in entry:
                self.assertIn(b"radial_point_count = MIN(220,max_radial_points)", output[entry])
        backup = grid.apply(root, changes, Path(self.temporary.name) / "suite-backup")
        self.assertEqual(grid.plan(root, settings), [])
        restoration = grid.restore_plan(backup, root)
        grid.write_changes(root, restoration)
        for path, data in originals.items():
            self.assertEqual(path.read_bytes(), data)

    def test_partial_layout_still_requires_every_supported_program(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        (root / grid.SUITE_GRID_FILES[-1]).unlink()
        with self.assertRaisesRegex(ValueError, "missing required source"):
            grid.plan(root, grid.Settings(nnnp=1990))

    def test_suite_real_entries_patch_and_restore_preserving_wrapper(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        repository = SCRIPT.parents[1]
        for relative in (
            grid.SUITE_CONFIG_FILE,
            grid.SUITE_EXPORT_FILE,
            *grid.SUITE_GRID_FILES,
            *grid.SUITE_RESTART_FILES,
        ):
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(repository / relative, path)
        originals = {p: p.read_bytes() for p in root.rglob("*.f90")}
        self.assertEqual(grid.plan(root, grid.Settings()), [])
        settings = grid.Settings(nnnp=1990, n=1179, h=0.025, accy=1e-10)
        changes = grid.plan(root, settings)
        self.assertEqual(
            len(changes), 1 + len(grid.SUITE_GRID_FILES) + len(grid.SUITE_RESTART_FILES)
        )
        self.assertNotIn(grid.SUITE_EXPORT_FILE, {c.relative for c in changes})
        for change in changes:
            if "/rmcdhf" in change.relative:
                self.assertEqual(change.after.lower().count(b"working_tolerance = 1d-10"), 2)
        backup = grid.apply(root, changes, Path(self.temporary.name) / "real-backup")
        self.assertEqual(grid.plan(root, settings), [])
        restoration = grid.restore_plan(backup, root)
        grid.write_changes(root, restoration)
        self.assertEqual({p: p.read_bytes() for p in originals}, originals)

    def test_suite_unknown_initializer_blocks_changes(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        unknown = root / "src/apps/new_program.f90"
        unknown.write_text("      RNT = 2.D-6/Z\n")
        with self.assertRaisesRegex(ValueError, "unrecognized grid initializer"):
            grid.plan(root, grid.Settings(nnnp=1990))

    def test_suite_requires_complete_reexports_and_single_definition(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        wrapper = root / grid.SUITE_EXPORT_FILE
        original = wrapper.read_text()
        incomplete, count = re.subn(
            r", max_orbitals_minus_two\b", "", original,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        self.assertEqual(count, 1, "fixture must remove exactly one required re-export")
        wrapper.write_text(incomplete)
        with self.assertRaisesRegex(ValueError, "must re-export all six"):
            grid.plan(root, grid.Settings(nnnp=1990))
        wrapper.write_text(original.replace("max_radial_points,", "max_radial_points => max_orbitals,"))
        with self.assertRaisesRegex(ValueError, "must use the semantic capacity names directly"):
            grid.plan(root, grid.Settings(nnnp=1990))
        wrapper.write_text(original)
        config = root / grid.SUITE_CONFIG_FILE
        config.write_text(config.read_text() + "integer, parameter :: radial_workspace_capacity = 600\n")
        with self.assertRaisesRegex(ValueError, "one central definition of radial_workspace_capacity"):
            grid.plan(root, grid.Settings(nnnp=1990))
        config.write_text(config.read_text().replace(
            "    integer, parameter :: radial_workspace_capacity = max_radial_points + 10\n", ""
        ))
        with self.assertRaisesRegex(ValueError, "must be derived as max_radial_points"):
            grid.plan(root, grid.Settings(nnnp=1990))

    def test_suite_rejects_old_capacity_exports(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        wrapper = root / grid.SUITE_EXPORT_FILE
        original = wrapper.read_text()
        for obsolete_export, message in (
            ("use compile_time_capacities, only: nnnp => max_radial_points", "semantic capacity names directly"),
            ("integer, parameter :: nnnp = max_radial_points", "only direct capacity re-exports"),
        ):
            with self.subTest(obsolete_export=obsolete_export):
                wrapper.write_text(original.replace("end module", f"    {obsolete_export}\nend module"))
                with self.assertRaisesRegex(ValueError, message):
                    grid.plan(root, grid.Settings(nnnp=1990))

    def test_suite_capacity_exports_honor_explicit_visibility(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        wrapper = root / grid.SUITE_EXPORT_FILE
        names = ", &\n        ".join(grid.SUITE_PARAMETER_NAMES)
        source = wrapper.read_text().replace("implicit none (type, external)", "implicit none(type,external)\n    private")
        source = source.replace("end module grasp_parameters", f"    public :: {names}\nend module")
        wrapper.write_text(source)
        changes = grid.plan(root, grid.Settings(nnnp=1990))
        self.assertEqual([change.relative for change in changes], [grid.SUITE_CONFIG_FILE])
        wrapper.write_text(source.replace(f"public :: {names}", "public :: max_radial_points"))
        with self.assertRaisesRegex(ValueError, "must re-export all six"):
            grid.plan(root, grid.Settings(nnnp=1990))

    @unittest.skipUnless(shutil.which(FORTRAN_COMPILER), "Fortran compiler is not available")
    def test_suite_capacity_edit_reaches_reexports_and_shared_arrays(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        source = SCRIPT.parents[1] / "src/lib/grasp_shared"
        module_names = [
            "compile_time_capacities.f90",
            "grasp_parameters.f90",
            "radial_grid_state.f90",
            "scf_radial_potentials.f90",
            "radial_quadrature_workspace.f90",
        ]
        for name in module_names:
            destination = root / "src/lib/grasp_shared" / name
            if not destination.exists():
                shutil.copyfile(source / name, destination)
        for capacity in (590, 1990):
            with self.subTest(capacity=capacity):
                changes = grid.plan(root, grid.Settings(nnnp=capacity))
                self.assertTrue(
                    all(c.relative == grid.SUITE_CONFIG_FILE for c in changes)
                )
                grid.write_changes(root, changes)
                build = Path(self.temporary.name) / f"compile-{capacity}"
                build.mkdir()
                result = subprocess.run(
                    [
                        FORTRAN_COMPILER,
                        "-o",
                        str(build / "capacity"),
                        *(str(root / "src/lib/grasp_shared" / name) for name in module_names),
                        str(SCRIPT.parents[1] / "test/shared_grid_capacity.f90"),
                    ],
                    cwd=build,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                result = subprocess.run(
                    [str(build / "capacity")],
                    cwd=build,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(
                    result.stdout.split()[-2:], [str(capacity), str(capacity + 10)]
                )

    def test_suite_duplicate_capacity_is_rejected(
        self,
    ) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        duplicate = root / "src/duplicate.f90"
        duplicate.write_text("      integer, parameter :: NNNP = 590\n")
        with self.assertRaisesRegex(ValueError, "duplicate capacity"):
            grid.plan(root, grid.Settings(h=0.025))

    def test_suite_local_math_step_is_not_a_radial_grid_default(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        local = root / "src/local_math.f90"
        local.write_text(
            "subroutine fit\nreal*8 :: h\nh = 0.01d0\nend subroutine fit\n"
        )
        changes = grid.plan(root, grid.Settings(h=0.025))
        self.assertEqual({c.relative for c in changes}, set(grid.SUITE_GRID_FILES))
        local.write_text(
            "use radial_grid_state, only: transformed_grid_step\n"
            "transformed_grid_step = 0.01d0\n"
        )
        with self.assertRaisesRegex(ValueError, "unrecognized grid initializer"):
            grid.plan(root, grid.Settings(h=0.025))

    def test_suite_fixed_accuracy_can_return_to_adaptive(self) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        fixed = grid.plan(root, grid.Settings(accy=1e-10))
        grid.write_changes(root, fixed)
        self.assertEqual(
            {c.relative for c in fixed},
            set(grid.SUITE_GRID_FILES + grid.SUITE_RESTART_FILES),
        )
        self.assertTrue(all(b"working_tolerance = 1d-10" in c.after for c in fixed))
        reset = grid.plan(root, grid.Settings(accy=0))
        self.assertTrue(all(b"working_tolerance = transformed_grid_step**6" in c.after for c in reset))
        grid.write_changes(root, reset)
        self.assertEqual(grid.plan(root, grid.Settings(accy=0)), [])

    def test_suite_capacity_reduction_rejects_existing_runtime_point_count(
        self,
    ) -> None:
        root = Path(self.temporary.name) / "suite-fixture"
        fixture(root)
        grid.write_changes(
            root, grid.plan(root, grid.Settings(nnnp=1990, n=1179))
        )
        with self.assertRaisesRegex(ValueError, "incompatible"):
            grid.plan(root, grid.Settings(nnnp=590))

    def test_finite_patch_preserves_point_model_and_comments(self) -> None:
        changes = grid.plan(
            self.root, grid.Settings(nnnp=1990, n=1179, h=0.025, rnt_scale=2e-6)
        )
        output = {change.relative: change.after.decode() for change in changes}
        for relative in grid.SUITE_GRID_FILES:
            self.assertIn("radial_point_count = 1179", output[relative])
            self.assertIn("transformed_grid_step = 2.5d-2", output[relative])
            self.assertIn("! max_radial_points = 590 is a historical comment", output[relative])
            if "/rci" not in relative:
                self.assertIn("radial_point_count = MIN(220,max_radial_points)", output[relative])
                self.assertIn("transformed_grid_step = 0.0625D0", output[relative])
        self.assertIn("max_radial_points = 1990", output[grid.SUITE_CONFIG_FILE])

    def test_apply_is_idempotent_and_restore_is_byte_exact(self) -> None:
        before = {
            path: (self.root / path).read_bytes()
            for path in grid.PATCHABLE_FILES
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
            if change.relative == "src/apps/rmcdhf/read_scf_options.f90"
        )
        self.assertEqual(rmcdhf_text.count("working_tolerance = 1d-10"), 2)
        restoration = grid.restore_plan(backup, self.root)
        grid.write_changes(self.root, restoration)
        for relative, data in before.items():
            self.assertEqual((self.root / relative).read_bytes(), data)

    def test_unknown_initializer_and_malformed_branch_are_rejected(self) -> None:
        path = self.root / "src/lib/new_grid.f90"
        path.write_text("      radial_grid_scale = 2e-6/nuclear_charge\n")
        with self.assertRaisesRegex(ValueError, "unrecognized grid initializer"):
            grid.plan(self.root, grid.Settings(h=0.025))
        path.write_text("! no active assignments\n")
        current = self.root / grid.SUITE_GRID_FILES[0]
        current.write_text(current.read_text().replace("transformed_grid_step = 0.05D0", "unknown_step = 0.05D0"))
        with self.assertRaisesRegex(ValueError, "expected one transformed_grid_step assignment"):
            grid.plan(self.root, grid.Settings(h=0.025))

    def test_invalid_numbers_and_runtime_capacity_conflict(self) -> None:
        for settings in (
            grid.Settings(nnnp=10),
            grid.Settings(nnnp=590, n=591),
            grid.Settings(h=float("nan")),
            grid.Settings(h=-0.01),
            grid.Settings(hp=-1),
            grid.Settings(accy=-1),
        ):
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                grid.plan(self.root, settings)
        path = self.root / grid.SUITE_GRID_FILES[0]
        path.write_text(path.read_text().replace("radial_point_count = max_radial_points", "radial_point_count = 1990"))
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
            if change.relative == grid.SUITE_GRID_FILES[0]
        )
        self.assertIn("radial_point_count = 400", text)
        self.assertIn("transformed_grid_step = 0.05D0", text)  # finite model unchanged
        self.assertIn("asymptotic_grid_step = 1d-1", text)

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
        changes = grid.plan(self.root, grid.Settings(nnnp=1990, n=1179))
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
