"""Keep the compiled OpenMP closure in sync with ordinary angular callees."""

from __future__ import annotations

import argparse
import json
import re
import shlex
from pathlib import Path

READ_ONLY_MODULES = {
    "cons_c",
    "facts_c",
    "orb_c",
    "stat_c",
    "debug_c",
    "iounit_c",
    "ribojj_c",
    "ribojj9_c",
    "ribojj11_c",
    "mtjj_c",
    "mtjj2_c",
    "vast_kind_param",
    "parameter_def",
}
PRIVATE_MODULES = {"m_c", "dumx_c", "trk_c", "buffer_c"}


def audit(source_root: Path, build: Path) -> dict[str, int]:
    files = {}
    for library in ("librang90", "lib9290", "libmcp90"):
        for path in (source_root / "src/lib" / library).glob("*.f90"):
            if not path.name.lower().endswith("_i.f90"):
                files[path.stem.lower()] = path
    source_text = {
        name: "\n".join(
            line.split("!")[0] for line in path.read_text().splitlines()
        ).lower()
        for name, path in files.items()
    }
    closure = set()
    pending = [
        "onescalar",
        "rkco_gg",
        "cord",
    ]  # Include the actual procedure callback explicitly.
    while pending:
        name = pending.pop()
        if name in closure:
            continue
        closure.add(name)
        direct_calls = set(re.findall(r"\bcall\s+(\w+)", source_text[name]))
        pending.extend(
            other
            for other in files
            if other not in closure
            and (
                other in direct_calls
                or re.search(r"\b" + re.escape(other) + r"\s*\(", source_text[name])
            )
        )
    manifest = (
        source_root / "src/appl/rangular/openmp/concurrent_sources.cmake"
    ).read_text()
    listed = {
        source_root / relative
        for relative in re.findall(r'\$\{PROJECT_SOURCE_DIR\}/([^"\n]+)', manifest)
    }
    expected = {files[name] for name in closure}
    if listed != expected:
        raise ValueError(
            f"Concurrent closure changed: missing={expected - listed}, extra={listed - expected}"
        )
    allowed = READ_ONLY_MODULES | PRIVATE_MODULES | {"memory_man"}
    for name in closure:
        for module in re.findall(
            r"^\s*use\s+([a-z0-9_]+)", source_text[name], re.MULTILINE
        ):
            if module not in allowed and not module.endswith("_i"):
                raise ValueError(f"New unaudited module {module} in {files[name]}")
        persistent_state_source = source_text[name]
        if name == "sixj":
            persistent_state_source = re.sub(
                r"^\s*save\s*::\s*raca\s*$",
                "",
                persistent_state_source,
                flags=re.MULTILINE,
            )
        if re.search(
            r"^\s*save\s+(?!\s*=)|,\s*save\b|^\s*common\b|^\s*(?:#\s*)?include\b",
            persistent_state_source,
            re.MULTILINE,
        ):
            raise ValueError(f"New persistent state requires review: {files[name]}")
    commands = json.loads((build / "compile_commands.json").read_text())
    compiled = set()
    for entry in commands:
        arguments = entry.get("arguments") or shlex.split(entry["command"])
        targets = set(re.findall(r"CMakeFiles/([^/]+)\.dir/", " ".join(arguments)))
        if targets & {
            "rangular_thread_state",
            "rangular_thread_kernels",
            "rangular_omp",
        }:
            if "-fopenmp" not in arguments or "-fno-automatic" in arguments:
                raise ValueError(
                    f"Incorrect concurrent compilation flags: {entry['file']}"
                )
            if "rangular_thread_state" not in targets:
                private_include = "-I" + str(build / "modules/rangular_thread_state")
                shared_include = "-I" + str(build / "modules/mod")
                if arguments.index(private_include) > arguments.index(shared_include):
                    raise ValueError(
                        "Thread-private module interfaces must precede the legacy module directory"
                    )
        if "rangular_thread_kernels" in targets:
            compiled.add(Path(entry["file"]).resolve())
    if compiled != {path.resolve() for path in expected}:
        raise ValueError(
            "Not every audited kernel is actually compiled in the OpenMP target"
        )
    for module in PRIVATE_MODULES:
        content = (
            source_root / "src/lib/libmod" / (module[:-2] + "_C.f90")
        ).read_text()
        if "!$OMP THREADPRIVATE" not in content:
            raise ValueError(f"Missing thread-private declaration: {module}")
    return {
        "concurrent_procedures": len(closure),
        "private_modules": len(PRIVATE_MODULES),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("build", type=Path)
    arguments = parser.parse_args()
    print(
        json.dumps(
            audit(Path(__file__).resolve().parents[2], arguments.build.resolve())
        )
    )
