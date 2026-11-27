#!/usr/bin/env python3
"""Preview, apply, check, or restore radial-grid edits to a grasp-atomic-suite tree.

This patches known source locations, not calculation files or installed binaries.
The script in scripts/ edits the src/ directory alongside scripts/.
Unrecognized layouts fail before any source file is written.
Run without arguments to read config.toml beside this script (Python >= 3.11).
Explicit command-line requests remain independent of that configuration.
Backups are created only when a backup directory is explicitly requested.
Capacity-only mode edits max_radial_points and the finite-nucleus point count.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
import difflib
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import tomllib


SUITE_GRID_FILES = (
    "src/apps/rwfnestimate/read_calculation_options.f90",
    "src/apps/rci/read_ci_options.f90",
    "src/apps/rci_mpi/read_ci_options.f90",
    "src/apps/biorthtransform/configure_radial_grid.f90",
    "src/apps/biorthtransform_mpi/configure_radial_grid.f90",
    "src/apps/rtransition/read_transition_options.f90",
    "src/apps/rtransition_phase/read_transition_options.f90",
    "src/apps/rtransition_mpi/read_transition_options.f90",
    "src/apps/rmcdhf/read_scf_options.f90",
    "src/apps/rmcdhf_mpi/read_scf_options_mpi.f90",
    "src/apps/rmcdhf_mem/read_scf_options.f90",
    "src/apps/rmcdhf_mem_mpi/read_scf_options_mpi.f90",
    "src/apps/hyperfine/read_hyperfine_options.f90",  # also compiled into gj and gj_mpi
    "src/apps/isotope/prepare_isotope_calculation.f90",
    "src/apps/rdensity/read_radial_density_input.f90",
)
SUITE_RESTART_FILES = ("src/lib/core/rci_restart_records.f90",)  # shared by rci and rci_mpi
SUITE_CONFIG_FILE = "src/lib/grasp_shared/cal_constants.f90"
SUITE_EXPORT_FILE = "src/lib/grasp_shared/grasp_parameters.f90"
PATCHABLE_FILES = (SUITE_CONFIG_FILE, *SUITE_GRID_FILES, *SUITE_RESTART_FILES)
SUITE_PARAMETER_NAMES = (
    "orbital_index_radix",
    "max_radial_points",
    "radial_workspace_capacity",
    "max_orbitals",
    "max_orbitals_minus_one",
    "max_orbitals_minus_two",
)
LEGACY_CAPACITY_RE = re.compile(r"\b(?:keyorb|nnnp|nnn1|nnnw|nnnwm1|nnnwm2)\b")
CAPACITY_RE = re.compile(
    r"\b(?:nnnp|nnn1|max_radial_points|radial_workspace_capacity)\s*=(?!>)",
    re.IGNORECASE,
)
MAX_RADIAL_POINTS_RE = re.compile(
    r"^[^!\n]*\bmax_radial_points[ \t]*=[ \t]*(?P<value>\d+)(?=[ \t]*(?:!|$))",
    re.IGNORECASE | re.MULTILINE,
)


@dataclass(frozen=True)
class Settings:
    nnnp: int | None = None
    n: int | None = None
    h: float | None = None
    rnt_scale: float | None = None
    hp: float | None = None
    accy: float | None = None
    point_n: int | None = None
    point_h: float | None = None
    point_rnt_scale: float | None = None


@dataclass(frozen=True)
class Change:
    relative: str
    before: bytes
    after: bytes


def active(line: str) -> str:
    """Drop Fortran comments for the known assignment/declaration locations."""
    return line.partition("!")[0].strip()


def assignment(line: str, name: str) -> bool:
    return re.match(rf"^{name}\s*=(?!>)", active(line), re.IGNORECASE) is not None


def set_assignment(line: str, name: str, expression: str) -> str:
    code, marker, comment = line.partition("!")
    if ";" in code or "&" in code:
        raise ValueError(f"unsupported continued/multiple assignment: {line.strip()}")
    indent = code[: len(code) - len(code.lstrip())]
    suffix = f" !{comment.rstrip()}" if marker else ""
    return f"{indent}{name} = {expression.lower()}{suffix}\n"


def real_literal(value: float) -> str:
    return format(Decimal(str(value)), "E").replace("E", "D")


def validate(settings: Settings, capacity: int) -> None:
    if capacity < 13:
        raise ValueError("max_radial_points must be at least 13")
    for name in ("n", "point_n"):
        value = getattr(settings, name)
        if value is not None and not 13 <= value <= capacity:
            raise ValueError(f"{name} must be between 13 and max_radial_points={capacity}")
    for name in ("h", "rnt_scale", "accy", "point_h", "point_rnt_scale", "hp"):
        value = getattr(settings, name)
        nonnegative = name in ("hp", "accy")
        if value is not None and (
            not math.isfinite(value) or value < 0 or (not nonnegative and value == 0)
        ):
            raise ValueError(
                f"{name} must be finite and {'nonnegative' if nonnegative else 'positive'}"
            )


def patch_grid(
    text: str,
    relative: str,
    settings: Settings,
    capacity: int,
    *,
    sync_accy: bool,
) -> str:
    lines = text.splitlines(keepends=True)
    starts = [i for i, line in enumerate(lines) if assignment(line, "radial_grid_scale")]
    expected = 1 if "/rci" in relative else 2
    if len(starts) != expected:
        raise ValueError(
            f"expected {expected} default radial_grid_scale assignments, found {len(starts)}"
        )
    if expected == 2:
        prefix = [active(line) for line in lines[: starts[0]] if active(line)]
        between = [
            active(line).upper()
            for line in lines[starts[0] : starts[1]]
            if active(line)
        ]
        if (
            not prefix
            or not re.fullmatch(
                r"if\s*\(nuclear_model_parameter_count\s*(?:==|\.eq\.)\s*0\)\s*then",
                prefix[-1],
                re.IGNORECASE,
            )
            or between[-1] != "ELSE"
        ):
            raise ValueError("unrecognized point/finite nucleus branch layout")
    for position, start in enumerate(starts):
        point = expected == 2 and position == 0
        values = {
            "radial_grid_scale": settings.point_rnt_scale if point else settings.rnt_scale,
            "transformed_grid_step": settings.point_h if point else settings.h,
            "radial_point_count": settings.point_n if point else settings.n,
        }
        end = (
            starts[position + 1]
            if position + 1 < len(starts)
            else min(start + 12, len(lines))
        )
        for name, value in values.items():
            matches = [i for i in range(start, end) if assignment(lines[i], name)]
            if len(matches) != 1:
                raise ValueError(
                    f"expected one {name} assignment in {'point' if point else 'finite'} branch"
                )
            i = matches[0]
            if value is not None:
                expression = (
                    str(int(value))
                    if name == "radial_point_count"
                    else real_literal(float(value))
                )
                if name == "radial_grid_scale":
                    expression += "/nuclear_charge"
                lines[i] = set_assignment(lines[i], name, expression)
            if name == "radial_point_count":
                expression = (
                    re.sub(r"\s+", "", active(lines[i]).split("=", 1)[1]).lower()
                )
                limited = re.fullmatch(r"min\((\d+),max_radial_points\)", expression)
                if expression == "max_radial_points":
                    actual_n = capacity
                elif limited:
                    actual_n = min(int(limited[1]), capacity)
                elif expression.isdigit():
                    actual_n = int(expression)
                else:
                    raise ValueError(f"unsupported radial_point_count expression: {expression}")
                if not 13 <= actual_n <= capacity:
                    raise ValueError(
                        f"existing radial_point_count={actual_n} is incompatible with max_radial_points={capacity}"
                    )
    hp_lines = [i for i, line in enumerate(lines) if assignment(line, "asymptotic_grid_step")]
    if len(hp_lines) != 1:
        raise ValueError("expected exactly one common HP assignment")
    if settings.hp is not None:
        i = hp_lines[0]
        lines[i] = set_assignment(lines[i], "asymptotic_grid_step", real_literal(settings.hp))

    accy_lines = [i for i, line in enumerate(lines) if assignment(line, "working_tolerance")]
    expected_accy = (1, 2) if "/rmcdhf" in relative else (1,)
    if len(accy_lines) not in expected_accy:
        raise ValueError("unexpected default ACCY assignments")
    if len(accy_lines) == 2:
        prompts = [
            i for i, line in enumerate(lines) if "Revise the default ACCY" in line
        ]
        if len(prompts) != 1 or accy_lines[1] != prompts[0] - 1:
            raise ValueError("unrecognized additional RMCDHF ACCY assignment")
    i = accy_lines[0]
    if settings.accy is not None:
        expression = (
            "transformed_grid_step**6"
            if settings.accy == 0
            else real_literal(settings.accy)
        )
        lines[i] = set_assignment(lines[i], "working_tolerance", expression)
    # RMCDHF computes ACCY before optional grid input. Recompute before its
    # ACCY prompt, so the user can still override ACCY afterwards.
    if "/rmcdhf" in relative and sync_accy:
        prompts = [
            i for i, line in enumerate(lines) if "Revise the default ACCY" in line
        ]
        if len(prompts) != 1:
            raise ValueError("expected exactly one RMCDHF ACCY prompt")
        i = prompts[0]
        expression = active(lines[accy_lines[0]]).split("=", 1)[1].strip()
        replacement = set_assignment(lines[accy_lines[0]], "working_tolerance", expression)
        if assignment(lines[i - 1], "working_tolerance"):
            lines[i - 1] = replacement
        else:
            lines.insert(i, replacement)
    return "".join(lines)


def source_path(root: Path, relative: str) -> Path:
    root = root.resolve()
    path = root / relative
    if path.is_symlink() or not path.resolve().is_relative_to(root):
        raise ValueError(f"source path escapes tree or is a symlink: {relative}")
    if not path.is_file():
        raise ValueError(f"missing required source: {relative}")
    return path


def validate_suite_exports(config: str, wrapper: str) -> None:
    """Require the suite's central definitions and six direct public imports."""
    names = set(SUITE_PARAMETER_NAMES)
    config_code = "\n".join(line.partition("!")[0].lower() for line in config.splitlines())
    lines = (line.partition("!")[0].strip().lower() for line in wrapper.splitlines())
    wrapper_code = re.sub(r"&\n&?", " ", "\n".join(line for line in lines if line))
    imported: set[str] = set()
    explicit_public: set[str] = set()
    explicit_private: set[str] = set()
    default_public = True
    for line in wrapper_code.splitlines():
        if re.fullmatch(r"(?:module\s+grasp_parameters|end(?:\s*module(?:\s+grasp_parameters)?)?"
                        r"|implicit\s+none(?:\s*\([\w\s,]+\))?)", line):
            continue
        if line in {"public", "private"}:
            default_public = line == "public"
            continue
        access = re.fullmatch(r"(public|private)\s*(?:::\s*|\s+)(.+)", line)
        if access:
            visible = explicit_public if access[1] == "public" else explicit_private
            visible.update(item.strip() for item in access[2].split(","))
            continue
        match = re.fullmatch(r"use\s*(?:::\s*|\s+)cal_constants\s*,\s*only\s*:\s*(.+)", line)
        if not match:
            raise ValueError("grasp_parameters may contain only direct capacity re-exports")
        imported.update(item.strip() for item in match[1].split(","))
    descriptor_constants = {
        "csf_total_j_slot", "csf_parity_slot", "csf_angular_code_count",
        "csf_orbital_letters", "maximum_orbital_two_j",
    }
    if not (imported | explicit_public | explicit_private) <= names | descriptor_constants:
        raise ValueError("grasp_parameters must use the semantic capacity names directly")
    exported = (imported if default_public else explicit_public) - explicit_private
    if not names <= imported & exported:
        raise ValueError("grasp_parameters must re-export all six cal_constants limits")

    legacy_name = LEGACY_CAPACITY_RE.search(config_code)
    if legacy_name:
        raise ValueError(f"legacy capacity name {legacy_name[0]} in cal_constants")
    for name in SUITE_PARAMETER_NAMES:
        if len(re.findall(rf"\b{name}\s*=(?!>)", config_code)) != 1:
            raise ValueError(f"expected one central definition of {name}")
    if not re.search(
        r"\bradial_workspace_capacity\s*=\s*max_radial_points\s*\+\s*10[ 	]*$",
        config_code,
        re.MULTILINE,
    ):
        raise ValueError("radial_workspace_capacity must be derived as max_radial_points + 10")


def plan(
    root: Path,
    settings: Settings,
    *,
    capacity_only: bool = False,
) -> list[Change]:
    if capacity_only:
        if settings.nnnp is None or settings.n is None:
            raise ValueError("--capacity-only requires both --nnnp and --n")
        if any(
            value is not None
            for name, value in vars(settings).items()
            if name not in ("nnnp", "n")
        ):
            raise ValueError("--capacity-only permits only nnnp and n parameters")

    required = (*PATCHABLE_FILES, SUITE_EXPORT_FILE)
    originals = {
        relative: source_path(root, relative).read_bytes() for relative in required
    }
    texts = {relative: data.decode("utf-8") for relative, data in originals.items()}
    if any("\r" in text for text in texts.values()):
        raise ValueError(
            "this patcher expects LF source files; convert line endings first"
        )
    config = texts[SUITE_CONFIG_FILE]
    validate_suite_exports(config, texts[SUITE_EXPORT_FILE])
    capacity_match = MAX_RADIAL_POINTS_RE.search(config)
    if not capacity_match:
        raise ValueError("cannot read integer max_radial_points")
    current_capacity = int(capacity_match["value"])
    capacity = settings.nnnp if settings.nnnp is not None else current_capacity
    validate(settings, capacity)
    if capacity != current_capacity:
        start, end = capacity_match.span("value")
        texts[SUITE_CONFIG_FILE] = config[:start] + str(capacity) + config[end:]

    sync_accy = any(
        value is not None for value in (settings.h, settings.point_h, settings.accy)
    )
    for relative in SUITE_GRID_FILES:
        try:
            texts[relative] = patch_grid(
                texts[relative], relative, settings, capacity, sync_accy=sync_accy
            )
        except ValueError as error:
            raise ValueError(f"{relative}: {error}") from error
    if settings.accy is not None:
        expression = (
            "transformed_grid_step**6"
            if settings.accy == 0
            else real_literal(settings.accy)
        )
        for relative in SUITE_RESTART_FILES:
            lines = texts[relative].splitlines(keepends=True)
            indices = [
                i for i, line in enumerate(lines) if assignment(line, "working_tolerance")
            ]
            if len(indices) != 1:
                raise ValueError(f"{relative}: expected one restart ACCY assignment")
            i = indices[0]
            lines[i] = set_assignment(lines[i], "working_tolerance", expression)
            texts[relative] = "".join(lines)

    # Reject unregistered grid initializers and duplicate capacity definitions.
    for path in sorted((root / "src").rglob("*.f90")):
        relative = path.relative_to(root).as_posix()
        text = texts.get(relative)
        if text is None:
            text = source_path(root, relative).read_text(encoding="utf-8")
        code_lines = [active(line) for line in text.splitlines()]
        code = "\n".join(code_lines)
        uses_grid = re.search(r"\buse\s+radial_grid_state\b", code, re.IGNORECASE)
        uses_accuracy = re.search(r"\buse\s+atomic_calculation_state\b", code, re.IGNORECASE)
        for number, line in enumerate(code_lines, 1):
            if relative != SUITE_CONFIG_FILE and CAPACITY_RE.search(line):
                raise ValueError(f"duplicate capacity: {relative}:{number}")
            if relative in SUITE_GRID_FILES:
                continue
            if (
                assignment(line, "radial_grid_scale")
                or assignment(line, "rnt")
                or (
                    uses_grid
                    and any(
                        assignment(line, name)
                        for name in ("transformed_grid_step", "asymptotic_grid_step")
                    )
                )
                or (
                    uses_accuracy
                    and relative not in SUITE_RESTART_FILES
                    and assignment(line, "working_tolerance")
                )
            ):
                raise ValueError(f"unrecognized grid initializer: {relative}:{number}")
    return [
        Change(relative, originals[relative], text.encode("utf-8"))
        for relative, text in texts.items()
        if text.encode("utf-8") != originals[relative]
    ]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def atomic_write(path: Path, data: bytes) -> None:
    mode = path.stat().st_mode & 0o777
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as output:
            temporary = Path(output.name)
            output.write(data)
        temporary.chmod(mode)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def write_changes(root: Path, changes: list[Change]) -> None:
    for change in changes:
        if source_path(root, change.relative).read_bytes() != change.before:
            raise ValueError(f"source changed after preview: {change.relative}")
    written: list[Change] = []
    try:
        for change in changes:
            atomic_write(root / change.relative, change.after)
            written.append(change)
    except Exception:
        for change in reversed(written):
            atomic_write(root / change.relative, change.before)
        raise


def apply(root: Path, changes: list[Change], backup: Path | None) -> Path | None:
    root = root.resolve()
    if backup is None:
        write_changes(root, changes)
        return None
    backup = backup.resolve()
    if backup.is_relative_to(root / "src"):
        raise ValueError("backup directory must be outside the source directory")
    backup.mkdir(parents=True, exist_ok=False)
    entries = []
    for change in changes:
        saved = backup / change.relative
        saved.parent.mkdir(parents=True, exist_ok=True)
        saved.write_bytes(change.before)
        entries.append(
            {
                "path": change.relative,
                "before": digest(change.before),
                "after": digest(change.after),
            }
        )
    manifest = {
        "root": str(root),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "files": entries,
    }
    (backup / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Backup: {backup}")
    write_changes(root, changes)
    return backup


def restore_plan(backup: Path, root: Path) -> list[Change]:
    manifest = json.loads((backup / "manifest.json").read_text(encoding="utf-8"))
    if Path(manifest["root"]).resolve() != root.resolve():
        raise ValueError("backup belongs to a different source tree")
    for relative in (SUITE_CONFIG_FILE, SUITE_EXPORT_FILE):
        source_path(root, relative)
    changes = []
    for entry in manifest["files"]:
        relative = entry["path"]
        if relative not in PATCHABLE_FILES:
            raise ValueError(f"backup contains a file outside the atomic-suite grid layout: {relative}")
        path = source_path(root, relative)
        saved_path = source_path(backup, relative)
        current, saved = path.read_bytes(), saved_path.read_bytes()
        if digest(saved) != entry["before"]:
            raise ValueError(f"backup checksum mismatch: {relative}")
        if digest(current) != entry["after"]:
            raise ValueError(
                f"source has subsequent edits; refusing to overwrite: {relative}"
            )
        changes.append(Change(relative, current, saved))
    return changes


def configured_arguments(configuration: Path) -> list[str]:
    """Translate the sibling TOML file into the same validated CLI request."""
    try:
        with configuration.open("rb") as source:
            config = tomllib.load(source)
    except FileNotFoundError as error:
        raise ValueError(
            f"missing configuration: {configuration}; copy config.example.toml "
            "to config.toml beside this script"
        ) from error
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"invalid TOML in {configuration}: {error}") from error
    unknown = config.keys() - {
        "RUN_MODE",
        "PRINT_DIFF",
        "CAPACITY_ONLY",
        "GRID_PARAMETERS",
        "BACKUP_DIRECTORY",
        "RESTORE_DIRECTORY",
    }
    if unknown:
        raise ValueError(f"unknown config.toml keys: {', '.join(sorted(unknown))}")
    mode = config.get("RUN_MODE", "preview")
    print_diff = config.get("PRINT_DIFF", True)
    capacity_only = config.get("CAPACITY_ONLY", False)
    if mode not in ("preview", "apply", "check"):
        raise ValueError("RUN_MODE must be preview, apply, or check")
    if not isinstance(print_diff, bool):
        raise ValueError("PRINT_DIFF must be true or false")
    if not isinstance(capacity_only, bool):
        raise TypeError("CAPACITY_ONLY must be true or false")

    def config_path(name: str) -> Path | None:
        value = config.get(name)
        if value is None:
            return None
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be a nonempty path string")
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = configuration.parent / path
        return path.resolve()

    arguments = []
    restore = config_path("RESTORE_DIRECTORY")
    if restore is not None:
        arguments.extend(("--restore", str(restore)))
    else:
        if capacity_only:
            arguments.append("--capacity-only")
        parameters = config.get("GRID_PARAMETERS", {})
        if not isinstance(parameters, dict):
            raise ValueError("GRID_PARAMETERS must be a TOML table: [GRID_PARAMETERS]")
        unknown = parameters.keys() - Settings.__dataclass_fields__.keys()
        if unknown:
            raise ValueError(
                f"unknown GRID_PARAMETERS keys: {', '.join(sorted(unknown))}"
            )
        for name, value in parameters.items():
            types = (int,) if name in ("nnnp", "n", "point_n") else (int, float)
            if type(value) not in types:
                raise ValueError(
                    f"GRID_PARAMETERS.{name} must be {'an integer' if types == (int,) else 'a number'}; omit it to preserve the current value"
                )
            arguments.extend((f"--{name.replace('_', '-')}", str(value)))
        configured_backup = config_path("BACKUP_DIRECTORY")
        if mode == "apply" and configured_backup is not None:
            arguments.extend(("--backup-dir", str(configured_backup)))
    if mode != "preview":
        arguments.append(f"--{mode}")
    if print_diff:
        arguments.append("--diff")
    return arguments


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nnnp",
        type=int,
        help="compile-time radial grid capacity (max_radial_points)",
    )
    parser.add_argument("--n", type=int, help="finite-nucleus runtime point count")
    parser.add_argument(
        "--capacity-only",
        action="store_true",
        help="edit only max_radial_points and the finite-nucleus point count; requires --nnnp and --n, "
        "rejects other grid parameters, and preserves existing initialization",
    )
    parser.add_argument(
        "--h", type=float, help="finite-nucleus exponential-coordinate step"
    )
    parser.add_argument(
        "--rnt-scale", type=float, help="finite-nucleus RNT coefficient: RNT=value/Z"
    )
    parser.add_argument(
        "--hp",
        type=float,
        help="common HP for both nuclear models; zero selects exponential grid",
    )
    parser.add_argument(
        "--accy",
        type=float,
        help="explicit grid-related tolerance; use 0 for adaptive H**6",
    )
    parser.add_argument(
        "--point-n",
        type=int,
        help="optional point-nucleus point count; the initial default is min(220,max_radial_points)",
    )
    parser.add_argument("--point-h", type=float, help="optional point-nucleus step")
    parser.add_argument(
        "--point-rnt-scale",
        type=float,
        help="optional point-nucleus RNT coefficient: value/Z",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--apply", action="store_true", help="write edits (default is preview only)"
    )
    action.add_argument(
        "--check",
        action="store_true",
        help="exit 1 if the requested patch is not already applied",
    )
    parser.add_argument("--diff", action="store_true", help="print full unified diffs")
    parser.add_argument(
        "--backup-dir",
        type=Path,
        help="optional new backup directory for --apply; no backup by default",
    )
    parser.add_argument(
        "--restore",
        type=Path,
        help="preview backup restoration; add --apply to restore",
    )
    script_directory = Path(__file__).resolve().parent
    configuration = script_directory / "config.toml"
    root = script_directory.parent
    arguments = list(sys.argv[1:] if argv is None else argv)
    file_config = not arguments
    if file_config:
        try:
            arguments = configured_arguments(configuration)
        except (ValueError, TypeError, OSError) as error:
            parser.error(str(error))
    args = parser.parse_args(arguments)
    settings = Settings(
        **{name: getattr(args, name) for name in Settings.__dataclass_fields__}
    )
    requested = any(value is not None for value in vars(settings).values())
    if args.restore and (
        requested or args.backup_dir or args.check or args.capacity_only
    ):
        parser.error(
            "--restore cannot be combined with parameter, capacity-only, backup-dir, or check options"
        )
    if not args.restore and not requested:
        parser.error(
            "provide at least one grid parameter, or --restore (omitted parameters stay unchanged)"
        )
    if args.backup_dir and not args.apply:
        parser.error("--backup-dir requires --apply")
    try:
        if args.restore:
            changes = restore_plan(args.restore.resolve(), root)
        else:
            changes = plan(root, settings, capacity_only=args.capacity_only)
        if file_config:
            print(f"Configuration: {configuration}")
        print(f"Source: {root}\nFiles requiring changes: {len(changes)}")
        for change in changes:
            print(f"  {change.relative}")
            if args.diff:
                print(
                    "".join(
                        difflib.unified_diff(
                            change.before.decode().splitlines(keepends=True),
                            change.after.decode().splitlines(keepends=True),
                            fromfile=change.relative,
                            tofile=change.relative,
                        )
                    ),
                    end="",
                )
        if args.apply and changes:
            if args.restore:
                write_changes(root, changes)
                print("Restored original sources.")
            else:
                apply(root, changes, args.backup_dir)
                print(
                    "Applied. Rebuild ALL libraries and programs, including MPI, and install them "
                    "using the grasp-atomic-suite README CMake workflow."
                )
        elif not args.apply:
            print(
                "Check only; no files written."
                if args.check
                else (
                    'Preview only; set RUN_MODE="apply" in config.toml to write.'
                    if file_config
                    else "Preview only; add --apply to write."
                )
            )
        if not args.restore:
            if args.capacity_only:
                print(
                    "Capacity-only: update max_radial_points and the finite-nucleus grid size; "
                    "preserve existing grid initialization and ACCY timing. "
                    "Earlier patches are not reverted."
                )
            print(
                "Scope: atomic-suite radial applications. Unspecified numerical values stay unchanged."
            )
            print(
                "Excluded: independent HF grids; calculation files; old RCI .res grid records; installed binaries."
            )
            print("Runtime input can override defaults.")
        return 1 if args.check and changes else 0
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
