"""Capture GNU Fortran diagnostics per compilation and enforce reviewed budgets.

Compiler-launcher records live in the build tree. The portable baseline records
diagnostics, not source snippets or machine paths. Known correctness warnings
remain unresolved migration work even when their regression budget is satisfied.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

BASELINE = Path(__file__).resolve().parents[2] / "docs/reference/f2018/gfortran15-warning-baseline.json"
AXES = ("serial", "mpi", "csfg", "openmp")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def case_id(profile: str, features: dict) -> str:
    result = profile.removeprefix("f2018-") + "-" + "-".join(
        f"{axis[0]}{int(bool(features[axis]))}" for axis in AXES)
    default_angular = bool(features["serial"] and features["openmp"])
    actual_angular = bool(features.get("rangular_openmp", default_angular))
    if actual_angular != default_angular:
        result += f"-a{int(actual_angular)}"
    return result


def command_parts(args: list[str], cwd: Path, source: Path, build: Path) -> tuple[str, str, list[str]]:
    if "-c" not in args or "-o" not in args:
        raise ValueError("Expected a Fortran compilation with -c and -o")
    input_file = (cwd / args[args.index("-c") + 1]).resolve()
    output_file = (cwd / args[args.index("-o") + 1]).resolve()
    relative_source = input_file.relative_to(source).as_posix()
    relative_object = output_file.relative_to(build).as_posix()
    # Paths in -I/-J can be embedded in a flag. Longest prefix first.
    replacements = sorted(((str(source), "<source>"), (str(build), "<build>")),
                          key=lambda pair: len(pair[0]), reverse=True)
    normalized = []
    for arg in args:
        for old, new in replacements:
            arg = arg.replace(old, new)
        normalized.append(arg)
    return relative_source, relative_object, normalized


def record_path(build: Path, object_name: str) -> Path:
    return build / "audit/warning-records" / (hashlib.sha256(object_name.encode()).hexdigest() + ".json")


def parse_diagnostics(stderr: str) -> list[dict]:
    # gfortran emits JSON arrays, possibly more than one from its subprocesses.
    decoder = json.JSONDecoder()
    remaining = stderr.strip()
    result = []
    while remaining:
        value, end = decoder.raw_decode(remaining)
        if not isinstance(value, list) or any(not isinstance(item, dict) or "kind" not in item for item in value):
            raise ValueError("Unexpected compiler diagnostic JSON")
        result.extend(value)
        remaining = remaining[end:].strip()
    return result


def capture(source: Path, build: Path, compiler: Path, command: list[str]) -> int:
    """CMake compiler launcher; compile normally, retaining every diagnostic."""
    source, build = source.resolve(), build.resolve()
    cwd = Path.cwd()
    # Preserve an existing launcher (for example ccache), but compare the
    # compiler invocation itself with compile_commands.json.
    compiler_command = command[command.index(str(compiler)):]
    filename, obj, normalized = command_parts(compiler_command, cwd, source, build)
    record = record_path(build, obj)
    record.unlink(missing_ok=True)
    result = subprocess.run([*command, "-fdiagnostics-format=json", "-fdiagnostics-color=never"],
                            capture_output=True, text=True, env=os.environ | {"LC_ALL": "C"})
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    try:
        diagnostics = parse_diagnostics(result.stderr)
    except (ValueError, json.JSONDecodeError) as error:
        (build / obj).unlink(missing_ok=True)
        print(f"Warning inventory: unparsed compiler output: {error}", file=sys.stderr)
        return result.returncode or 1
    if result.returncode:
        (build / obj).unlink(missing_ok=True)
        return result.returncode
    write_json(record, {"schema_version": 1, "source": filename, "object": obj,
                       "command": normalized, "source_sha256": digest(source / filename),
                       "object_sha256": digest(build / obj), "diagnostics": diagnostics})
    return 0


def category(option: str, message: str) -> str:
    if option in {"-Wuninitialized", "-Wmaybe-uninitialized", "-Warray-bounds", "-Wdo-subscript",
                  "-Wcharacter-truncation", "-Wconversion", "-Wconversion-extra", "-Waliasing",
                  "-Wsurprising", "-Wargument-mismatch", "-Wreturn-type"} or re.search(
            r"uninitialized|out of bounds|(?:INTENT|Type|Rank) mismatch|truncat|overflow|larger than limit",
            message, re.I):
        return "correctness-review"
    if "Obsolescent feature" in message or "obsolescent feature" in message:
        return "obsolescent-standard"
    return "general"


def inventory(build: Path) -> dict:
    build = build.resolve()
    features = json.loads((build / "f2018-features.json").read_text())
    profile = features.pop("profile")
    if profile == "f2018-release":
        raise ValueError("Production builds do not establish a warning baseline; use f2018-audit")
    cache = dict(re.findall(r"^([^/#][^:\n]*):[^=\n]+=(.*)$", (build / "CMakeCache.txt").read_text(), re.M))
    source = Path(cache["CMAKE_HOME_DIRECTORY"]).resolve()
    compiler = cache["CMAKE_Fortran_COMPILER"]
    version = subprocess.check_output([compiler, "-dumpfullversion"], text=True).strip()
    from check import check_commands, check_feature_commands
    commands = json.loads((build / "compile_commands.json").read_text())
    check_commands(commands, profile, Path(compiler))
    check_feature_commands(commands, features)
    records, files, diagnostics, counts, target_counts = [], {}, {}, Counter(), Counter()
    for entry in commands:
        if Path(entry["file"]).suffix.lower() not in {".f90", ".f", ".for", ".f95"}:
            continue
        args = entry.get("arguments") or shlex.split(entry["command"])
        filename, obj, normalized = command_parts(args, Path(entry["directory"]), source, build)
        path = record_path(build, obj)
        if not path.is_file():
            raise ValueError(f"Missing warning record: {filename} ({obj}); rebuild with the strict launcher")
        record = json.loads(path.read_text())
        if (record["command"] != normalized or record["source_sha256"] != digest(source / filename)
                or record["object_sha256"] != digest(build / obj)):
            raise ValueError(f"Stale warning record: {filename} ({obj}); rebuild")
        records.append(obj)
        files[filename] = record["source_sha256"]
        target_match = re.search(r"CMakeFiles/([^/]+)\.dir/", obj)
        if target_match is None:
            raise ValueError(f"Cannot identify compilation target: {obj}")
        target_counts[target_match.group(1)] += 1

        def visit(items: list[dict]) -> None:
            for diagnostic in items:
                kind = diagnostic["kind"]
                if kind not in {"warning", "note"}:
                    raise ValueError(f"Unexpected diagnostic kind in a successful compilation: {kind}")
                if kind == "warning":
                    locations = diagnostic.get("locations", [])
                    caret = locations[0].get("caret", {}) if locations else {}
                    location = caret.get("file", filename)
                    if location.startswith("<"):
                        location = filename
                    location = (source / location).resolve().relative_to(source).as_posix()
                    message = diagnostic["message"].replace(str(build), "<build>").replace(str(source), "<source>")
                    item = {"file": location, "target": target_match.group(1), "line": caret.get("line", 0),
                            "column": caret.get("column", 0), "option": diagnostic.get("option", ""),
                            "message": message}
                    key = hashlib.sha256(json.dumps(item, sort_keys=True, ensure_ascii=False).encode()).hexdigest()[:20]
                    item["category"] = category(item["option"], message)
                    if key in diagnostics and diagnostics[key] != item:
                        raise ValueError("Diagnostic identity collision")
                    diagnostics[key] = item
                    counts[key] += 1
                visit(diagnostic.get("children", []))
        visit(record["diagnostics"])
    if not records or len(records) != len(set(records)):
        raise ValueError("Missing or duplicate Fortran compilation inventory")
    categories = Counter()
    for key, count in counts.items():
        categories[diagnostics[key]["category"]] += count
    return {"schema_version": 1, "case": case_id(profile, features), "profile": profile,
            "features": features, "compiler_version": version, "compilation_count": len(records),
            "target_compilation_counts": dict(sorted(target_counts.items())),
            "files": dict(sorted(files.items())), "diagnostics": dict(sorted(diagnostics.items())),
            "counts": dict(sorted(counts.items())), "category_counts": dict(sorted(categories.items()))}


def compare(observation: dict, baseline: dict) -> dict:
    case = observation["case"]
    if baseline.get("schema_version") != 1:
        raise ValueError("Unsupported warning baseline schema")
    if observation["compiler_version"] != baseline["compiler_version"]:
        raise ValueError("Warning baseline compiler version differs; establish a separate reviewed baseline")
    if case not in baseline["cases"]:
        raise ValueError(f"No reviewed warning baseline for {case}")
    expected = baseline["cases"][case]
    if set(observation["files"]) != set(expected["files"]):
        raise ValueError(f"Warning source coverage changed for {case}; review the baseline explicitly")
    if observation["compilation_count"] != expected["compilation_count"]:
        raise ValueError(f"Warning compilation coverage changed for {case}")
    if observation["target_compilation_counts"] != expected["target_compilation_counts"]:
        raise ValueError(f"Warning target coverage changed for {case}")
    increases = {key: count - expected["counts"].get(key, 0) for key, count in observation["counts"].items()
                 if count > expected["counts"].get(key, 0)}
    if increases:
        details = [{**observation["diagnostics"][key], "increase": count} for key, count in increases.items()]
        raise ValueError(f"New/increased warnings for {case}: {json.dumps(details, ensure_ascii=False)}")
    return {"status": "passed", "warning_count": sum(observation["counts"].values()),
            "category_counts": observation["category_counts"],
            "correctness_warnings_clear": observation["category_counts"].get("correctness-review", 0) == 0}


def make_baseline(observations: list[dict]) -> dict:
    versions = {item["compiler_version"] for item in observations}
    if len(versions) != 1:
        raise ValueError("Baseline requires exactly one compiler version")
    cases, diagnostics, source_files = {}, {}, {}
    for item in observations:
        if item["case"] in cases:
            raise ValueError("Duplicate baseline configuration")
        for name, sha in item["files"].items():
            if name in source_files and source_files[name] != sha:
                raise ValueError(f"Source changed between matrix builds: {name}")
            source_files[name] = sha
        cases[item["case"]] = {"profile": item["profile"], "features": item["features"],
                              "compilation_count": item["compilation_count"],
                              "target_compilation_counts": item["target_compilation_counts"],
                              "files": sorted(item["files"]), "counts": item["counts"]}
        diagnostics.update(item["diagnostics"])
    return {"schema_version": 1, "compiler_version": versions.pop(),
            "scope": "GNU Fortran CMake compilation diagnostics; known correctness warnings remain unresolved",
            "source_files": dict(sorted(source_files.items())),
            "diagnostics": dict(sorted(diagnostics.items())), "cases": dict(sorted(cases.items()))}


def render_report(baseline: dict, output: Path) -> None:
    """Human-readable per-file index, including sources with zero warnings."""
    by_file = {name: {} for name in baseline["source_files"]}
    for key, diagnostic in baseline["diagnostics"].items():
        by_file.setdefault(diagnostic["file"], {})[key] = diagnostic
    lines = ["# GNU Fortran 15 逐文件警告基线", "",
             f"编译器：{baseline['compiler_version']}；配置数：{len(baseline['cases'])}。", "",
             "机器可读清单见 [gfortran15-warning-baseline.json](gfortran15-warning-baseline.json)。",
             "该表由 `test/f2018/warning_inventory.py report` 生成，包含无警告源文件。",
             "唯一诊断按目标/文件/行/列/选项/消息区分；O0/O2 次数为已登记 Check/Audit 配置中的预算最大次数，",
             "包含同一源文件编入多个目标的重复。正确性复核项尚未逐条解决，登记不是豁免。", "",
             "“未初始化”指变量或数组元素被读取时尚未获得确定值；Fortran 声明本身不保证零值。",
             "`-Wmaybe-uninitialized` 也可能表示编译器无法证明赋值路径或分支不可达，须结合代码和有效输入复核。",
             "诊断数不等于变量数或已确认缺陷数；修复须符合实际数据契约，不能统一补零消除告警。", "",
             "| 文件 | 唯一诊断 | 正确性复核 | 过时标准语法 | 一般 | O0 最大次数 | O2 最大次数 |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    source = Path(__file__).resolve().parents[2]
    for name, items in sorted(by_file.items()):
        categories = Counter(item["category"] for item in items.values())
        maxima = [max((sum(case["counts"].get(key, 0) for key in items) for case in baseline["cases"].values()
                       if case["profile"] == profile), default=0) for profile in ("f2018-check", "f2018-audit")]
        link = Path(os.path.relpath(source / name, output.resolve().parent)).as_posix()
        lines.append(f"| [{name}]({link}) | {len(items)} | {categories['correctness-review']} | "
                     f"{categories['obsolescent-standard']} | {categories['general']} | {maxima[0]} | {maxima[1]} |")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    launcher = sub.add_parser("capture")
    launcher.add_argument("--source", type=Path, required=True)
    launcher.add_argument("--build", type=Path, required=True)
    launcher.add_argument("--compiler", type=Path, required=True)
    launcher.add_argument("command", nargs=argparse.REMAINDER)
    audit = sub.add_parser("audit")
    audit.add_argument("build", type=Path)
    audit.add_argument("--baseline", type=Path, default=BASELINE)
    audit.add_argument("--record-only", action="store_true")
    record = sub.add_parser("record", help="Explicitly create a candidate baseline; never auto-update it")
    record.add_argument("matrix_root", type=Path)
    record.add_argument("--output", type=Path, required=True)
    report = sub.add_parser("report")
    report.add_argument("--baseline", type=Path, default=BASELINE)
    report.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.action == "capture":
            command = args.command[1:] if args.command[:1] == ["--"] else args.command
            return capture(args.source, args.build, args.compiler, command)
        if args.action == "record":
            from matrix import cases
            observations = [inventory(args.matrix_root / case["id"]) for case in cases()
                            if case["profile"] != "f2018-release"]
            if args.output.exists():
                raise ValueError("Refusing to overwrite a baseline; write a separate candidate and review its diff")
            baseline = make_baseline(observations)
            source = Path(__file__).resolve().parents[2]
            baseline["source_revision"] = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
            baseline["source_state"] = "Working tree; exact compiled source identities are in source_files"
            write_json(args.output, baseline)
            print(f"Recorded {len(observations)} configurations in {args.output}")
        elif args.action == "report":
            render_report(json.loads(args.baseline.read_text()), args.output)
        else:
            observation = inventory(args.build)
            write_json(args.build / "audit/warnings.json", observation)
            if not args.record_only:
                print(json.dumps(compare(observation, json.loads(args.baseline.read_text()))))
        return 0
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
