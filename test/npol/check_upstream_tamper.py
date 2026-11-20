"""Negative checks for the reviewed provenance chain, on disposable copies."""
from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import shutil
import tempfile

from check_upstream import main as verify
from check_upstream import APP_SYMBOLS_SNAPSHOT, suite_source_path


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    manifest = json.loads((root / "docs/reference/rmcdhf_upstream_snapshot.json").read_text())
    with tempfile.TemporaryDirectory(prefix="upstream-tamper-") as temporary:
        work = Path(temporary)
        migration = json.loads((root / "docs/reference/rmcdhf_descriptor_snapshot.json").read_text())
        app_symbols = json.loads((root / APP_SYMBOLS_SNAPSHOT).read_text())
        renamed_paths = {name: record["suite_path"] for name, record in app_symbols["files"].items()}
        paths = [
            suite_source_path(root, name, renamed_paths)
            for name in manifest["files"] if name not in migration["removed"]
        ]
        paths += list((root / "docs/reference").glob("rmcdhf_*.json"))
        for source in paths:
            target = work / source.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        with contextlib.redirect_stdout(io.StringIO()):
            verify(work)
        for relative in ("src/apps/rmcdhf/basis_csf_two_j_plus_one.f90", "src/apps/rmcdhf/load_csf_basis_data.f90", "src/apps/rmcdhf/run_scf_iterations.f90", "src/apps/rmcdhf/solve_radial_orbital.f90",
                         "src/apps/rmcdhf/direct_angular_coefficient.f90",
                         "src/apps/rmcdhf/serial_mpi_context.f90",
                         "src/apps/rmcdhf/rmcdhf_main.f90",
                         "src/apps/rmcdhf/diagonalize_hamiltonian_block.f90",
                         "src/apps/rmcdhf/CMakeLists.txt",
                         "src/apps/rmcdhf_mpi/Makefile",
                         "src/apps/rmcdhf_mpi/read_al_average_options_mpi.f90",
                         "src/apps/rmcdhf/write_input_summary.f90",
                         "docs/reference/rmcdhf_io_units_snapshot.json",
                         "docs/reference/rmcdhf_app_variable_names_snapshot.json",
                         "docs/reference/rmcdhf_app_symbols_snapshot.json",
                         "docs/reference/rmcdhf_grasp_mpi_naming_snapshot.json",
                         "docs/reference/rmcdhf_grasp_shared_naming_snapshot.json",
                         "docs/reference/rmcdhf_davidson_naming_snapshot.json",
                         "docs/reference/rmcdhf_kind_alias_cleanup_snapshot.json",
                         "docs/reference/rmcdhf_csf_queries_snapshot.json",
                         "docs/reference/rmcdhf_descriptor_snapshot.json",
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
