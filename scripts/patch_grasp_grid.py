#!/usr/bin/env python3
"""Preview, apply, check, or restore radial-grid edits to a GRASP2018 tree.

This patches known source locations, not calculation files or installed binaries.
Unrecognized layouts fail before any source file is written.
Run without arguments to read config.toml beside this script (Python >= 3.11).
Explicit command-line requests remain independent of that configuration.
Backups are created only when a backup directory is explicitly requested.
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


GRID_FILES = (
    "src/appl/rwfnestimate90/getinfo.f90",
    "src/appl/rmcdhf90/getscd.f90",
    "src/appl/rmcdhf90_mpi/getscdmpi.f90",
    "src/appl/rmcdhf90_mem/getscd.f90",
    "src/appl/rmcdhf90_mem_mpi/getscdmpi.f90",
    "src/appl/rci90/getcid.f90",
    "src/appl/rci90_mpi/getcid.f90",
    "src/appl/rbiotransform90/radpar.f90",
    "src/appl/rbiotransform90_mpi/radparmpi.f90",
    "src/appl/rtransition90/getosd.f90",
    "src/appl/rtransition90_phase/getosd.f90",
    "src/appl/rtransition90_mpi/getosdmpi.f90",
    "src/appl/rhfs90/gethfd.f90",
    "src/appl/rhfszeeman95/gethfd.f90",
    "src/appl/sms90/getsmd.f90",
    "src/appl/ris4/getsmd.f90",
    "src/appl/rdensity/getsmd.f90",
    "src/tool/rwfnrotate.f90",
    "src/tool/rwfnrelabel.f90",
)
CAPACITY_FILES = (
    "src/lib/libmod/parameter_def_M.f90",
    "src/appl/rwfnestimate90/frmtfp.f90",
    "src/appl/rwfnestimate90/solvh_I.f90",
    "src/appl/rwfnestimate90/tail_I.f90",
    "src/appl/rwfnestimate90/frmrwf_I.f90",
    "src/appl/rwfnestimate90/summry_I.f90",
)
RESTART_FILES = ("src/appl/rci90/lodres.f90", "src/appl/rci90_mpi/lodres.f90")
SUITE_GRID_FILES = (
    "src/appl/rmcdhf90/getscd.f90",
    "src/appl/rmcdhf90_mpi/getscdmpi.f90",
    "src/appl/rmcdhf90_mem/getscd.f90",
    "src/appl/rmcdhf90_mem_mpi/getscdmpi.f90",
    "src/appl/rhfs90/gethfd.f90",  # also compiled into gj90 and gj90_mpi
    "src/appl/ris4/getsmd.f90",
    "src/appl/rdensity/getsmd.f90",
)
SUITE_CONFIG_FILE = "src/lib/libmod/suite_parameters_M.f90"
SUITE_DEFAULTS_FILE = "src/lib/libmod/radial_grid_defaults_M.f90"
CAPACITY_RE = re.compile(r"\b(NNNP|NNN1)\s*=\s*([^)!\n]+)", re.IGNORECASE)


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
    return line.split("!", 1)[0].strip()


def assignment(line: str, name: str) -> bool:
    return re.match(rf"^{name}\s*=", active(line), re.IGNORECASE) is not None


def set_assignment(line: str, name: str, expression: str) -> str:
    code, marker, comment = line.partition("!")
    if ";" in code or "&" in code:
        raise ValueError(f"unsupported continued/multiple assignment: {line.strip()}")
    indent = code[: len(code) - len(code.lstrip())]
    suffix = f" !{comment.rstrip()}" if marker else ""
    return f"{indent}{name} = {expression}{suffix}\n"


def real_literal(value: float) -> str:
    return format(Decimal(str(value)), "E").replace("E", "D")


def validate(
    settings: Settings, capacity: int, *, allow_adaptive_accy: bool = False
) -> None:
    if capacity < 13:
        raise ValueError("NNNP must be at least 13")
    for name in ("n", "point_n"):
        value = getattr(settings, name)
        if value is not None and not 13 <= value <= capacity:
            raise ValueError(f"{name} must be between 13 and NNNP={capacity}")
    for name in ("h", "rnt_scale", "accy", "point_h", "point_rnt_scale", "hp"):
        value = getattr(settings, name)
        nonnegative = name == "hp" or (name == "accy" and allow_adaptive_accy)
        if value is not None and (
            not math.isfinite(value) or value < 0 or (not nonnegative and value == 0)
        ):
            raise ValueError(
                f"{name} must be finite and {'nonnegative' if nonnegative else 'positive'}"
            )


def patch_capacity(text: str, capacity: int) -> str:
    lines = text.splitlines(keepends=True)
    hits = 0
    for i, line in enumerate(lines):
        if not active(line):
            continue

        def replace(match: re.Match[str]) -> str:
            nonlocal hits
            hits += 1
            name = match[1].upper()
            return f"{name} = {capacity if name == 'NNNP' else capacity + 10}"

        lines[i] = CAPACITY_RE.sub(replace, line)
    if not hits:
        raise ValueError("expected radial capacity declarations are missing")
    return "".join(lines)


def patch_grid(text: str, relative: str, settings: Settings, capacity: int) -> str:
    lines = text.splitlines(keepends=True)
    starts = [i for i, line in enumerate(lines) if assignment(line, "RNT")]
    expected = 1 if "/rci90" in relative else 2
    if len(starts) != expected:
        raise ValueError(
            f"expected {expected} default RNT assignments, found {len(starts)}"
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
                r"IF\s*\(NPARM\s*(==|\.EQ\.)\s*0\)\s*THEN", prefix[-1], re.IGNORECASE
            )
            or between[-1] != "ELSE"
        ):
            raise ValueError("unrecognized point/finite nucleus branch layout")
    for position, start in enumerate(starts):
        point = expected == 2 and position == 0
        values = {
            "RNT": settings.point_rnt_scale if point else settings.rnt_scale,
            "H": settings.point_h if point else settings.h,
            "N": settings.point_n if point else settings.n,
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
                    str(int(value)) if name == "N" else real_literal(float(value))
                )
                if name == "RNT":
                    expression += "/Z"
                lines[i] = set_assignment(lines[i], name, expression)
            if name == "N":
                expression = (
                    active(lines[i]).split("=", 1)[1].strip().upper().replace(" ", "")
                )
                if expression == "NNNP":
                    actual_n = capacity
                elif re.fullmatch(r"MIN\(\d+,NNNP\)", expression):
                    actual_n = min(int(expression[4:].split(",")[0]), capacity)
                elif expression.isdigit():
                    actual_n = int(expression)
                else:
                    raise ValueError(f"unsupported N expression: {expression}")
                if not 13 <= actual_n <= capacity:
                    raise ValueError(
                        f"existing N={actual_n} is incompatible with NNNP={capacity}"
                    )
    hp_lines = [i for i, line in enumerate(lines) if assignment(line, "HP")]
    if len(hp_lines) != 1:
        raise ValueError("expected exactly one common HP assignment")
    if settings.hp is not None:
        i = hp_lines[0]
        lines[i] = set_assignment(lines[i], "HP", real_literal(settings.hp))

    accy_lines = [i for i, line in enumerate(lines) if assignment(line, "ACCY")]
    expected_accy = (1, 2) if "/rmcdhf90" in relative else (1,)
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
        lines[i] = set_assignment(lines[i], "ACCY", real_literal(settings.accy))
    # RMCDHF computes ACCY before optional grid input. Recompute before its
    # ACCY prompt, so the user can still override ACCY afterwards.
    if "/rmcdhf90" in relative:
        prompts = [
            i for i, line in enumerate(lines) if "Revise the default ACCY" in line
        ]
        if len(prompts) != 1:
            raise ValueError("expected exactly one RMCDHF ACCY prompt")
        i = prompts[0]
        expression = active(lines[accy_lines[0]]).split("=", 1)[1].strip()
        replacement = set_assignment(lines[accy_lines[0]], "ACCY", expression)
        if assignment(lines[i - 1], "ACCY"):
            lines[i - 1] = replacement
        else:
            lines.insert(i, replacement)
    return "".join(lines)


def initialize_relabel(text: str, reference: str) -> str:
    """Restore grid defaults and import the shared nuclear-model selector."""
    if not any(assignment(line, "RNT") for line in text.splitlines()):
        start = re.search(r"^!\s*IF\s*\(NPARM\b.*$", text, re.MULTILINE | re.IGNORECASE)
        end = re.search(r"^!\s*HP\s*=.*$", text, re.MULTILINE | re.IGNORECASE)
        block = re.search(
            r"^[ \t]*IF \(NPARM == 0\) THEN\n.*?^[ \t]*HP\s*=[^\n]*",
            reference,
            re.MULTILINE | re.DOTALL,
        )
        if not start or not end or not block or end.start() < start.start():
            raise ValueError("unrecognized rwfnrelabel initialization block")
        text = text[: start.start()] + block[0].lstrip("\n") + text[end.end() :]

    # SETISO loads NPARM in npar_C. A local INTEGER would compile but would
    # leave the nucleus branch controlled by an uninitialized, unrelated value.
    # Check even an already-active block to repair trees made by older patchers.
    routine = re.search(
        r"^[ \t]*SUBROUTINE\s+GETHFD\s*\([^\n)]*\)[^\n]*\n"
        r"(?P<preamble>.*?)^(?P<indent>[ \t]*)IMPLICIT\b",
        text,
        re.MULTILINE | re.DOTALL | re.IGNORECASE,
    )
    if not routine:
        raise ValueError("unrecognized rwfnrelabel GETHFD declarations")
    for line in routine["preamble"].splitlines():
        imported = re.fullmatch(
            r"USE\s+npar_C\s*(?:,\s*ONLY\s*:\s*(.*))?",
            active(line),
            re.IGNORECASE,
        )
        if imported and (
            imported[1] is None
            or "NPARM" in [name.strip().upper() for name in imported[1].split(",")]
        ):
            return text
    position = routine.start("indent")
    return (
        text[:position]
        + routine["indent"]
        + "USE npar_C, ONLY: NPARM\n"
        + text[position:]
    )


def source_path(root: Path, relative: str) -> Path:
    root = root.resolve()
    path = root / relative
    if path.is_symlink() or not path.resolve().is_relative_to(root):
        raise ValueError(f"source path escapes tree or is a symlink: {relative}")
    if not path.is_file():
        raise ValueError(f"missing required source: {relative}")
    return path


def plan_suite(root: Path, settings: Settings) -> list[Change]:
    """Edit one configuration file and audit every consumer before writing."""
    path = source_path(root, SUITE_CONFIG_FILE)
    before = path.read_bytes()
    text = before.decode("utf-8")
    if "\r" in text:
        raise ValueError("this patcher expects LF source files")
    names = (
        "NNNP",
        "NNN1",
        "NNNWM1",
        "NNNWM2",
        "FINITE_N",
        "POINT_N",
        "FINITE_H",
        "POINT_H",
        "FINITE_RNT_SCALE",
        "POINT_RNT_SCALE",
        "DEFAULT_HP",
        "DEFAULT_ACCY",
    )
    expressions: dict[str, str] = {}
    for name in names:
        matches = [
            match
            for line in text.splitlines()
            if (
                match := re.search(
                    rf"\b{name}\s*=\s*(.+)$", active(line), re.IGNORECASE
                )
            )
        ]
        if len(matches) != 1:
            raise ValueError(f"expected one central definition of {name}")
        expressions[name] = matches[0][1].strip()
    for name, expression in {
        "NNN1": "NNNP+10",
        "NNNWM1": "NNNW-1",
        "NNNWM2": "NNNW-2",
    }.items():
        if re.sub(r"\s+", "", expressions[name]).upper() != expression:
            raise ValueError(f"{name} must remain derived as {expression}")
    capacity = settings.nnnp if settings.nnnp is not None else int(expressions["NNNP"])
    validate(settings, capacity, allow_adaptive_accy=True)
    requested = {
        "NNNP": settings.nnnp,
        "FINITE_N": settings.n,
        "POINT_N": settings.point_n,
        "FINITE_H": settings.h,
        "POINT_H": settings.point_h,
        "FINITE_RNT_SCALE": settings.rnt_scale,
        "POINT_RNT_SCALE": settings.point_rnt_scale,
        "DEFAULT_HP": settings.hp,
        "DEFAULT_ACCY": settings.accy,
    }
    for name, value in requested.items():
        if value is None:
            continue
        expression = (
            str(int(value))
            if name in ("NNNP", "FINITE_N", "POINT_N")
            else real_literal(value)
        )
        pattern = re.compile(rf"(\b{name}\s*=\s*)([^!\n]+)", re.IGNORECASE)
        lines = text.splitlines(keepends=True)
        for i, line in enumerate(lines):
            if re.search(rf"\b{name}\s*=", active(line), re.IGNORECASE):
                lines[i] = pattern.sub(
                    lambda m: m[1] + expression + (" " if "!" in line else ""), line
                )
        text = "".join(lines)
        expressions[name] = expression

    def point_count(expression: str) -> int:
        compact = re.sub(r"\s+", "", expression).upper()
        if compact == "NNNP":
            return capacity
        if compact.isdigit():
            return int(compact)
        if match := re.fullmatch(r"MIN\((\d+),NNNP\)", compact):
            return min(int(match[1]), capacity)
        raise ValueError(f"unsupported central point count: {expression}")

    def real_value(expression: str) -> float:
        compact = re.sub(r"\s+", "", expression).upper()
        if compact == "EXP(-65.0D0/16.0D0)":
            return math.exp(-65 / 16)
        return float(compact.replace("D", "E"))

    validate(
        Settings(
            n=point_count(expressions["FINITE_N"]),
            point_n=point_count(expressions["POINT_N"]),
            h=real_value(expressions["FINITE_H"]),
            point_h=real_value(expressions["POINT_H"]),
            rnt_scale=real_value(expressions["FINITE_RNT_SCALE"]),
            point_rnt_scale=real_value(expressions["POINT_RNT_SCALE"]),
            hp=real_value(expressions["DEFAULT_HP"]),
            accy=real_value(expressions["DEFAULT_ACCY"]),
        ),
        capacity,
        allow_adaptive_accy=True,
    )
    source_path(root, SUITE_DEFAULTS_FILE)
    for relative in SUITE_GRID_FILES:
        consumer = source_path(root, relative).read_text(encoding="utf-8")
        calls = [
            line
            for line in consumer.splitlines()
            if re.fullmatch(
                r"CALL\s+SET_RADIAL_DEFAULTS\s*\(\s*NPARM\s*,\s*Z\s*\)",
                active(line),
                re.IGNORECASE,
            )
        ]
        if len(calls) != 1:
            raise ValueError(f"{relative}: expected one shared radial initializer")
    for source in sorted((root / "src").rglob("*.f90")):
        relative = source.relative_to(root).as_posix()
        content = (
            text
            if relative == SUITE_CONFIG_FILE
            else source_path(root, relative).read_text(encoding="utf-8")
        )
        code = "\n".join(active(line) for line in content.splitlines())
        uses_grid = re.search(r"\bUSE\s+GRID_C\b", code, re.IGNORECASE) is not None
        uses_accuracy = re.search(r"\bUSE\s+DEF_C\b", code, re.IGNORECASE) is not None
        for number, line in enumerate(content.splitlines(), 1):
            if relative != SUITE_CONFIG_FILE and CAPACITY_RE.search(active(line)):
                raise ValueError(f"duplicate capacity: {relative}:{number}")
            if relative != SUITE_DEFAULTS_FILE and (
                assignment(line, "RNT")
                or (uses_grid and any(assignment(line, name) for name in ("H", "HP")))
                or (uses_accuracy and assignment(line, "ACCY"))
            ):
                raise ValueError(f"unrecognized grid initializer: {relative}:{number}")
    after = text.encode("utf-8")
    return [Change(SUITE_CONFIG_FILE, before, after)] if before != after else []


def plan(root: Path, settings: Settings, layout: str = "grasp2018") -> list[Change]:
    if layout == "grasp2018":
        capacity_files, grid_files, restart_files = (
            CAPACITY_FILES,
            GRID_FILES,
            RESTART_FILES,
        )
    elif layout == "atomic-suite":
        return plan_suite(root, settings)
    else:
        raise ValueError(f"unsupported layout: {layout}")
    required = (*capacity_files, *grid_files, *restart_files)
    originals = {
        relative: source_path(root, relative).read_bytes() for relative in required
    }
    texts = {relative: data.decode("utf-8") for relative, data in originals.items()}
    if any("\r" in text for text in texts.values()):
        raise ValueError(
            "this patcher expects LF source files; convert line endings first"
        )
    capacity_match = re.search(
        r"\bNNNP\s*=\s*(\d+)", texts[CAPACITY_FILES[0]], re.IGNORECASE
    )
    if not capacity_match:
        raise ValueError("cannot read global NNNP")
    capacity = settings.nnnp if settings.nnnp is not None else int(capacity_match[1])
    validate(settings, capacity)
    for relative in capacity_files:
        texts[relative] = patch_capacity(texts[relative], capacity)
    relabel = "src/tool/rwfnrelabel.f90"
    if relabel in texts:
        texts[relabel] = initialize_relabel(texts[relabel], texts[GRID_FILES[0]])
    for relative in grid_files:
        try:
            texts[relative] = patch_grid(texts[relative], relative, settings, capacity)
        except ValueError as error:
            raise ValueError(f"{relative}: {error}") from error
    if settings.accy is not None:
        for relative in restart_files:
            lines = texts[relative].splitlines(keepends=True)
            indices = [i for i, line in enumerate(lines) if assignment(line, "ACCY")]
            if len(indices) != 1:
                raise ValueError(f"{relative}: expected one restart ACCY assignment")
            i = indices[0]
            lines[i] = set_assignment(lines[i], "ACCY", real_literal(settings.accy))
            texts[relative] = "".join(lines)
    # Inspect the entire tree, not just files we intend to edit. Unknown radial
    # capacity/default locations may be an upstream version we do not support.
    for path in sorted((root / "src").rglob("*.f90")):
        relative = path.relative_to(root).as_posix()
        text = texts.get(relative)
        if text is None:
            text = source_path(root, relative).read_text(encoding="utf-8")
        for number, line in enumerate(text.splitlines(), 1):
            code = active(line)
            for match in CAPACITY_RE.finditer(code):
                expression = match[2].replace(" ", "").upper()
                wanted = capacity if match[1].upper() == "NNNP" else capacity + 10
                if expression != str(wanted) and not (
                    match[1].upper() == "NNN1" and expression == "NNNP+10"
                ):
                    raise ValueError(f"unmatched capacity: {relative}:{number}: {code}")
            if (
                assignment(line, "RNT")
                and relative not in grid_files
                and relative != "src/tool/rwfnmchfmcdf.f90"
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


def restore_plan(backup: Path, root: Path | None) -> tuple[Path, list[Change]]:
    manifest = json.loads((backup / "manifest.json").read_text(encoding="utf-8"))
    recorded_root = Path(manifest["root"]).resolve()
    if root is not None and root.resolve() != recorded_root:
        raise ValueError("--grasp does not match backup's source tree")
    changes = []
    for entry in manifest["files"]:
        relative = entry["path"]
        path = source_path(recorded_root, relative)
        saved_path = source_path(backup, relative)
        current, saved = path.read_bytes(), saved_path.read_bytes()
        if digest(saved) != entry["before"]:
            raise ValueError(f"backup checksum mismatch: {relative}")
        if digest(current) != entry["after"]:
            raise ValueError(
                f"source has subsequent edits; refusing to overwrite: {relative}"
            )
        changes.append(Change(relative, current, saved))
    return recorded_root, changes


def configured_arguments() -> list[str]:
    """Translate the sibling TOML file into the same validated CLI request."""
    configuration = Path(__file__).resolve().with_name("config.toml")
    try:
        with configuration.open("rb") as source:
            config = tomllib.load(source)
    except FileNotFoundError as error:
        raise ValueError(
            f"missing configuration: {configuration}; copy config.example.toml "
            "to config.toml beside this script and set GRASP_SOURCE"
        ) from error
    except tomllib.TOMLDecodeError as error:
        raise ValueError(f"invalid TOML in {configuration}: {error}") from error
    unknown = config.keys() - {
        "GRASP_SOURCE",
        "SOURCE_LAYOUT",
        "RUN_MODE",
        "PRINT_DIFF",
        "GRID_PARAMETERS",
        "BACKUP_DIRECTORY",
        "RESTORE_DIRECTORY",
    }
    if unknown:
        raise ValueError(f"unknown config.toml keys: {', '.join(sorted(unknown))}")
    mode = config.get("RUN_MODE", "preview")
    layout = config.get("SOURCE_LAYOUT", "grasp2018")
    print_diff = config.get("PRINT_DIFF", True)
    if mode not in ("preview", "apply", "check"):
        raise ValueError("RUN_MODE must be preview, apply, or check")
    if layout not in ("grasp2018", "atomic-suite"):
        raise ValueError("SOURCE_LAYOUT must be grasp2018 or atomic-suite")
    if not isinstance(print_diff, bool):
        raise ValueError("PRINT_DIFF must be true or false")

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

    arguments = ["--layout", layout]
    root = config_path("GRASP_SOURCE")
    restore = config_path("RESTORE_DIRECTORY")
    if root is not None:
        arguments.extend(("--grasp", str(root)))
    if restore is not None:
        arguments.extend(("--restore", str(restore)))
    else:
        if root is None:
            raise ValueError(
                "set GRASP_SOURCE in config.toml to your GRASP source directory"
            )
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
        "--grasp",
        type=Path,
        help="GRASP2018 source root (default: workspace sibling grasp)",
    )
    parser.add_argument(
        "--layout",
        choices=("grasp2018", "atomic-suite"),
        default="grasp2018",
        help="required source layout; atomic-suite defaults to this repository",
    )
    parser.add_argument(
        "--nnnp",
        type=int,
        help="compiled capacity; also sets every duplicate NNN1 to NNNP+10",
    )
    parser.add_argument("--n", type=int, help="finite-nucleus runtime point count")
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
        help="explicit grid-related tolerance; atomic-suite also accepts 0 for adaptive H**6",
    )
    parser.add_argument(
        "--point-n",
        type=int,
        help="optional point-nucleus point count; otherwise preserve MIN(220,NNNP)",
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
    arguments = list(sys.argv[1:] if argv is None else argv)
    file_config = not arguments
    if file_config:
        try:
            arguments = configured_arguments()
        except (ValueError, TypeError, OSError) as error:
            parser.error(str(error))
    args = parser.parse_args(arguments)
    settings = Settings(
        **{name: getattr(args, name) for name in Settings.__dataclass_fields__}
    )
    requested = any(value is not None for value in vars(settings).values())
    if args.restore and (requested or args.backup_dir or args.check):
        parser.error(
            "--restore cannot be combined with parameter, backup-dir, or check options"
        )
    if not args.restore and not requested:
        parser.error(
            "provide at least one grid parameter, or --restore (omitted parameters stay unchanged)"
        )
    if args.backup_dir and not args.apply:
        parser.error("--backup-dir requires --apply")
    try:
        root = args.grasp.resolve() if args.grasp else None
        if args.restore:
            root, changes = restore_plan(args.restore.resolve(), root)
        else:
            root = root or (
                Path(__file__).resolve().parents[1]
                if args.layout == "atomic-suite"
                else Path(__file__).resolve().parents[2] / "grasp"
            )
            root = root.resolve()
            changes = plan(root, settings, args.layout)
        if file_config:
            print(f"Configuration: {Path(__file__).resolve().with_name('config.toml')}")
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
                    "using the GRASP README CMake workflow."
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
            print(
                f"Scope: {args.layout} radial applications. Unspecified numerical values stay unchanged."
            )
            print(
                "Excluded: independent HF/wfnplot/rwfnmchfmcdf grids; calculation files; old RCI .res grid records; installed binaries."
            )
            print("Runtime input can override defaults.")
            if args.layout == "grasp2018":
                print(
                    "Without --rnt-scale, rwfnrotate retains its historical RNT=2e-6 difference."
                )
        return 1 if args.check and changes else 0
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
