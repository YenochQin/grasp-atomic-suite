"""Negative checks for the reviewed provenance chain, on disposable copies."""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import shutil
import tempfile

from check_upstream import main as verify


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    manifest = json.loads((root / "docs/reference/rmcdhf_upstream_snapshot.json").read_text())
    with tempfile.TemporaryDirectory(prefix="upstream-tamper-") as temporary:
        work = Path(temporary)
        paths = [root / name for name in manifest["files"]]
        paths += list((root / "docs/reference").glob("rmcdhf_*.json"))
        for source in paths:
            target = work / source.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        with contextlib.redirect_stdout(io.StringIO()):
            verify(work)
        for relative in ("src/appl/rmcdhf/scf.f90", "src/appl/rmcdhf/solve.f90",
                         "docs/reference/rmcdhf_declarations_snapshot.json",
                         "docs/reference/rmcdhf_numeric_kind_snapshot.json",
                         "docs/reference/rmcdhf_spacing_snapshot.json",
                         "docs/reference/rmcdhf_quote_end_snapshot.json",
                         "docs/reference/rmcdhf_decimal_snapshot.json",
                         "docs/reference/rmcdhf_parser_style_snapshot.json",
                         "docs/reference/rmcdhf_operator_style_snapshot.json",
                         "docs/reference/rmcdhf_use_only_snapshot.json",
                         "docs/reference/rmcdhf_use_only_state_snapshot.json",
                         "docs/reference/rmcdhf_use_only_state2_snapshot.json",
                         "docs/reference/rmcdhf_use_only_state3_snapshot.json",
                         "docs/reference/rmcdhf_use_only_debug_snapshot.json",
                         "docs/reference/rmcdhf_use_only_state4_snapshot.json",
                         "docs/reference/rmcdhf_use_only_pack_snapshot.json",
                         "docs/reference/rmcdhf_use_only_pos_snapshot.json",
                         "docs/reference/rmcdhf_use_only_eigv_snapshot.json",
                         "docs/reference/rmcdhf_use_only_radial_snapshot.json",
                         "docs/reference/rmcdhf_use_only_scf_snapshot.json",
                         "docs/reference/rmcdhf_use_only_memory_snapshot.json",
                         "docs/reference/rmcdhf_use_only_def_orb_snapshot.json",
                         "docs/reference/rmcdhf_use_only_orb_label_snapshot.json"):
            target = work / relative
            original = target.read_bytes()
            target.write_bytes(original + b"\n")
            try:
                verify(work)
            except AssertionError:
                print(f"Rejected tampering: {relative}")
            else:
                raise AssertionError(f"Unreviewed change accepted: {relative}")
            finally:
                target.write_bytes(original)


if __name__ == "__main__":
    main()
