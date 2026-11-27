"""Verify imported sources against the reviewed upstream and rename snapshots.

This check is independent of the external grasp checkout and never executes SCF.
It protects SOLVE/IMPROV/SCF/orthogonalization and the other imported routines from
silently reverting to the older orbopt implementation.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
import sys


APP_SYMBOLS_SNAPSHOT = Path("docs/reference/rmcdhf_app_symbols_snapshot.json")
APP_MODULES_SNAPSHOT = Path("docs/reference/rmcdhf_app_modules_snapshot.json")
SEMANTIC_NAMES_SNAPSHOT = Path("docs/reference/rmcdhf_semantic_names_snapshot.json")
AUDITED_NAMES_SNAPSHOT = Path("docs/reference/rmcdhf_audited_names_snapshot.json")
DIAGNOSTIC_LABELS_SNAPSHOT = Path("docs/reference/rmcdhf_diagnostic_labels_snapshot.json")
V3_J_LABELS_SNAPSHOT = Path("docs/reference/rmcdhf_v3_j_labels_snapshot.json")
APP_COMMENTS_SNAPSHOT = Path("docs/reference/rmcdhf_app_comments_snapshot.json")
MPI_RUN_LOGGING_SNAPSHOT = Path("docs/reference/rmcdhf_mpi_run_logging_snapshot.json")
MPI_MCP_IO_SNAPSHOT = Path("docs/reference/rmcdhf_mpi_mcp_io_snapshot.json")
RUN_TIMER_SNAPSHOT = Path("docs/reference/rmcdhf_run_timer_snapshot.json")
REQUIRED_ASF_ROOTS_SNAPSHOT = Path("docs/reference/rmcdhf_required_asf_roots_snapshot.json")
RADIAL_STREAM_SNAPSHOT = Path("docs/reference/rmcdhf_radial_stream_snapshot.json")
MIXING_STREAM_SNAPSHOT = Path("docs/reference/rmcdhf_mixing_stream_snapshot.json")
MPI_F08_SNAPSHOT = Path("docs/reference/rmcdhf_mpi_f08_snapshot.json")
SPARSE_STORAGE_SNAPSHOT = Path("docs/reference/rmcdhf_sparse_storage_snapshot.json")
BASIS_OWNERSHIP_SNAPSHOT = Path("docs/reference/rmcdhf_basis_ownership_snapshot.json")
HOST_MEMORY_SNAPSHOT = Path("docs/reference/rmcdhf_optional_host_memory_snapshot.json")
STATE_INPUT_SNAPSHOT = Path("docs/reference/rmcdhf_required_state_input_snapshot.json")


def suite_source_path(root: Path, upstream_path: str, renamed_paths: dict[str, str]) -> Path:
    # Historical snapshot keys keep the original source directory and filenames.
    return root / renamed_paths.get(upstream_path, upstream_path.replace("src/appl/", "src/apps/", 1))


def reverse_reviewed_edits(text: str, path: Path, edits: list[dict[str, str | int]]) -> str:
    for edit in reversed(edits):
        before, after, count = edit["before"], edit["after"], edit["count"]
        assert isinstance(before, str) and isinstance(after, str) and isinstance(count, int)
        pattern = "(?m)^" + re.escape(after)
        assert len(re.findall(pattern, text)) == count, f"Reviewed edit missing or ambiguous: {path}"
        text = re.sub(pattern, lambda _: before, text)
    return text


def original_bytes(path: Path, edits: list[dict[str, str | int]] | None = None) -> bytes:
    text = reverse_reviewed_edits(path.read_text(), path, edits or [])
    if path.name in ("getscd.f90", "getscdmpi.f90", "getscd_I.f90", "getscdmpi_I.f90"):
        text = text.replace("ISOFILE, RWFFILE, NPFILE)", "ISOFILE, RWFFILE)")
        text = text.replace("      CHARACTER (LEN = *), INTENT(IN) :: NPFILE\n", "")
        text = text.replace("      CHARACTER(*), INTENT(IN) :: NPFILE\n", "")
        text = text.replace("h, hp, n, rnt, r", "h, hp, n, rnt")
        imports = "      USE npot_C, ONLY: ZZ\n      USE npol_C, ONLY: NPOL_MESSAGE_LEN\n"
        for module in ("      USE npol_initialization, ONLY: NPOL_INIT\n",
                       "      USE npol_mpi, ONLY: NPOL_INIT_MPI, NPOL_ABORT_MPI\n"):
            text = text.replace(imports + module, "      !USE npot_C\n")
        text = text.replace("      INTEGER :: NP_STATUS\n      CHARACTER(NPOL_MESSAGE_LEN) :: NP_MESSAGE\n", "")
        text = text.replace(
            "      CALL NPOL_INIT(NPFILE, Z, N, R, ZZ, NP_STATUS, NP_MESSAGE)\n"
            "      IF (NP_STATUS /= 0) THEN\n         WRITE(ISTDE,'(A)') TRIM(NP_MESSAGE)\n"
            "         ERROR STOP 1\n      END IF\n", "")
        text = text.replace(
            "      CALL NPOL_INIT_MPI(NPFILE, Z, N, R, ZZ, MPI_COMM_WORLD, NP_STATUS, NP_MESSAGE)\n"
            "      IF (NP_STATUS /= 0) CALL NPOL_ABORT_MPI(MPI_COMM_WORLD, NP_MESSAGE)\n", "")
    elif path.name == "rscfvu.f90":
        text = text.replace("'isodata', 'rwfn.inp', 'npol.inp')", "'isodata', 'rwfn.inp')")
    elif path.name == "rscfmpivu.f90":
        text = text.replace("      USE npol_mpi, ONLY: NPOL_ABORT_MPI\n", "")
        text = text.replace("      CHARACTER(:), ALLOCATABLE :: npfile\n", "")
        text = text.replace(
            "      npfile = ''\n      IF (myid == 0) THEN\n"
            "         IF (lenperm == LEN(permdir)) CALL NPOL_ABORT_MPI(MPI_COMM_WORLD, &\n"
            "              'npol: permanent input directory reaches the legacy path limit')\n"
            "         npfile = TRIM(permdir) // '/npol.inp'\n      END IF\n"
            "      CALL GETSCDmpi (EOL, idblk, file1, file2, npfile)\n",
            "      CALL GETSCDmpi (EOL, idblk, file1, file2 )\n")
    elif path.name == "dacon.f90":
        text = text.replace("      USE dpbdt_I\n", "")
    else:
        return text.encode() if edits else path.read_bytes()
    return text.encode()


def main(root: Path) -> None:
    upstream_path = root / "docs/reference/rmcdhf_upstream_snapshot.json"
    compatibility_path = root / "docs/reference/rmcdhf_f2018_compatibility.json"
    manifest = json.loads(upstream_path.read_text())
    expected = manifest["files"]
    compatibility = json.loads(compatibility_path.read_text())["files"]
    assert set(compatibility).issubset(expected), "F2018 edits must refer to imported files"
    # Each before-suite identity was checked through original_bytes before the
    # case-only conversion. Pin the exact resulting bytes, not a relaxed hash
    # that would also ignore changes to string literals or comments.
    lowercase_path = root / "docs/reference/rmcdhf_lowercase_snapshot.json"
    lowercase = json.loads(lowercase_path.read_text())
    assert lowercase["schema_version"] == 1, "Unknown lowercase snapshot schema"
    assert hashlib.sha256(upstream_path.read_bytes()).hexdigest() == lowercase["upstream_manifest_sha256"], \
        "Upstream manifest changed; review lowercase provenance"
    assert hashlib.sha256(compatibility_path.read_bytes()).hexdigest() == lowercase["f2018_compatibility_sha256"], \
        "F2018 edits changed; review lowercase provenance"
    assert set(lowercase["files"]) == {name for name in expected if name.endswith(".f90")}, \
        "Lowercase snapshot must cover exactly the imported Fortran files"
    renamed_path = root / "docs/reference/rmcdhf_variable_rename_snapshot.json"
    renamed = json.loads(renamed_path.read_text())
    assert renamed["schema_version"] == 1, "Unknown variable-rename snapshot schema"
    assert hashlib.sha256(lowercase_path.read_bytes()).hexdigest() == renamed["previous_manifest_sha256"], \
        "Lowercase snapshot changed; review variable-rename provenance"
    assert renamed["files"].keys() == lowercase["files"].keys(), "Variable-rename snapshot coverage changed"
    for name, record in renamed["files"].items():
        assert record["before_suite_sha256"] == lowercase["files"][name]["after_suite_sha256"], \
            f"Variable-rename snapshot has a different predecessor: {name}"
    safety_path = root / "docs/reference/rmcdhf_input_safety_snapshot.json"
    safety = json.loads(safety_path.read_text())
    assert safety["schema_version"] == 1, "Unknown input-safety snapshot schema"
    assert hashlib.sha256(renamed_path.read_bytes()).hexdigest() == safety["previous_manifest_sha256"], \
        "Variable-rename snapshot changed; review input-safety provenance"
    assert set(safety["files"]).issubset(renamed["files"]), "Input-safety snapshot contains unrelated files"
    for name, record in safety["files"].items():
        assert record["before_suite_sha256"] == renamed["files"][name]["after_suite_sha256"], \
            f"Input-safety snapshot has a different predecessor: {name}"
    declarations = json.loads((root / "docs/reference/rmcdhf_declarations_snapshot.json").read_text())
    assert declarations["schema_version"] == 1, "Unknown declaration snapshot schema"
    assert hashlib.sha256(safety_path.read_bytes()).hexdigest() == declarations["previous_manifest_sha256"], \
        "Input-safety snapshot changed; review declaration provenance"
    prior_fortran = renamed["files"] | safety["files"]
    assert set(declarations["files"]).issubset(prior_fortran), "Declaration snapshot contains unrelated files"
    for name, record in declarations["files"].items():
        assert record["before_suite_sha256"] == prior_fortran[name]["after_suite_sha256"], \
            f"Declaration snapshot has a different predecessor: {name}"
    reviewed_fortran = prior_fortran | declarations["files"]
    numeric_path = root / "docs/reference/rmcdhf_numeric_kind_snapshot.json"
    numeric = json.loads(numeric_path.read_text())
    assert numeric["schema_version"] == 1, "Unknown numeric-kind snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_declarations_snapshot.json").read_bytes()).hexdigest() == \
        numeric["previous_manifest_sha256"], "Declaration snapshot changed; review numeric-kind provenance"
    numeric_scope = {f"src/appl/{name}/{routine}.f90" for name, routine in (
        ("rmcdhf", "scf"), ("rmcdhf_mem", "scf"),
        ("rmcdhf_mpi", "scfmpi"), ("rmcdhf_mem_mpi", "scfmpi"))}
    assert set(numeric["files"]) == numeric_scope, "Numeric-kind snapshot coverage changed"
    for name, record in numeric["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"Numeric-kind snapshot has a different predecessor: {name}"
    reviewed_fortran |= numeric["files"]
    spacing = json.loads((root / "docs/reference/rmcdhf_spacing_snapshot.json").read_text())
    assert spacing["schema_version"] == 1, "Unknown spacing snapshot schema"
    assert hashlib.sha256(numeric_path.read_bytes()).hexdigest() == spacing["previous_manifest_sha256"], \
        "Numeric-kind snapshot changed; review spacing provenance"
    assert set(spacing["files"]).issubset(reviewed_fortran), "Spacing snapshot contains unrelated files"
    for name, record in spacing["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"Spacing snapshot has a different predecessor: {name}"
    reviewed_fortran |= spacing["files"]
    quote_end = json.loads((root / "docs/reference/rmcdhf_quote_end_snapshot.json").read_text())
    assert quote_end["schema_version"] == 1, "Unknown quote/end snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_spacing_snapshot.json").read_bytes()).hexdigest() == \
        quote_end["previous_manifest_sha256"], "Spacing snapshot changed; review quote/end provenance"
    assert set(quote_end["files"]).issubset(reviewed_fortran), "Quote/end snapshot contains unrelated files"
    for name, record in quote_end["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"Quote/end snapshot has a different predecessor: {name}"
    reviewed_fortran |= quote_end["files"]
    decimal = json.loads((root / "docs/reference/rmcdhf_decimal_snapshot.json").read_text())
    assert decimal["schema_version"] == 1, "Unknown decimal snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_quote_end_snapshot.json").read_bytes()).hexdigest() == \
        decimal["previous_manifest_sha256"], "Quote/end snapshot changed; review decimal provenance"
    assert set(decimal["files"]).issubset(reviewed_fortran), "Decimal snapshot contains unrelated files"
    for name, record in decimal["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"Decimal snapshot has a different predecessor: {name}"
    reviewed_fortran |= decimal["files"]
    parser_style = json.loads((root / "docs/reference/rmcdhf_parser_style_snapshot.json").read_text())
    assert parser_style["schema_version"] == 1, "Unknown parser-style snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_decimal_snapshot.json").read_bytes()).hexdigest() == \
        parser_style["previous_manifest_sha256"], "Decimal snapshot changed; review parser-style provenance"
    assert set(parser_style["files"]).issubset(reviewed_fortran), "Parser-style snapshot contains unrelated files"
    for name, record in parser_style["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"Parser-style snapshot has a different predecessor: {name}"
    reviewed_fortran |= parser_style["files"]
    operator_style = json.loads((root / "docs/reference/rmcdhf_operator_style_snapshot.json").read_text())
    assert operator_style["schema_version"] == 1, "Unknown operator-style snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_parser_style_snapshot.json").read_bytes()).hexdigest() == \
        operator_style["previous_manifest_sha256"], "Parser-style snapshot changed; review operator-style provenance"
    assert set(operator_style["files"]).issubset(reviewed_fortran), "Operator-style snapshot contains unrelated files"
    for name, record in operator_style["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"Operator-style snapshot has a different predecessor: {name}"
    reviewed_fortran |= operator_style["files"]
    use_only = json.loads((root / "docs/reference/rmcdhf_use_only_snapshot.json").read_text())
    assert use_only["schema_version"] == 1, "Unknown USE ONLY snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_operator_style_snapshot.json").read_bytes()).hexdigest() == \
        use_only["previous_manifest_sha256"], "Operator-style snapshot changed; review USE ONLY provenance"
    assert set(use_only["files"]).issubset(reviewed_fortran), "USE ONLY snapshot contains unrelated files"
    for name, record in use_only["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"USE ONLY snapshot has a different predecessor: {name}"
    reviewed_fortran |= use_only["files"]
    use_state = json.loads((root / "docs/reference/rmcdhf_use_only_state_snapshot.json").read_text())
    assert use_state["schema_version"] == 1, "Unknown state USE ONLY snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_snapshot.json").read_bytes()).hexdigest() == \
        use_state["previous_manifest_sha256"], "USE ONLY snapshot changed; review state import provenance"
    assert set(use_state["files"]).issubset(reviewed_fortran), "State USE ONLY snapshot contains unrelated files"
    for name, record in use_state["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"State USE ONLY snapshot has a different predecessor: {name}"
    reviewed_fortran |= use_state["files"]
    use_state2 = json.loads((root / "docs/reference/rmcdhf_use_only_state2_snapshot.json").read_text())
    assert use_state2["schema_version"] == 1, "Unknown state2 USE ONLY snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_state_snapshot.json").read_bytes()).hexdigest() == \
        use_state2["previous_manifest_sha256"], "State USE ONLY snapshot changed; review state2 import provenance"
    assert set(use_state2["files"]).issubset(reviewed_fortran), "State2 USE ONLY snapshot contains unrelated files"
    for name, record in use_state2["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"State2 USE ONLY snapshot has a different predecessor: {name}"
    reviewed_fortran |= use_state2["files"]
    use_state3 = json.loads((root / "docs/reference/rmcdhf_use_only_state3_snapshot.json").read_text())
    assert use_state3["schema_version"] == 1, "Unknown state3 USE ONLY snapshot schema"
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_state2_snapshot.json").read_bytes()).hexdigest() == \
        use_state3["previous_manifest_sha256"], "State2 USE ONLY snapshot changed; review state3 import provenance"
    assert set(use_state3["files"]).issubset(reviewed_fortran), "State3 USE ONLY snapshot contains unrelated files"
    for name, record in use_state3["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], \
            f"State3 USE ONLY snapshot has a different predecessor: {name}"
    reviewed_fortran |= use_state3["files"]
    use_debug = json.loads((root / "docs/reference/rmcdhf_use_only_debug_snapshot.json").read_text())
    assert use_debug["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_state3_snapshot.json").read_bytes()).hexdigest() == use_debug["previous_manifest_sha256"]
    assert set(use_debug["files"]).issubset(reviewed_fortran)
    for name, record in use_debug["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_debug["files"]
    use_state4 = json.loads((root / "docs/reference/rmcdhf_use_only_state4_snapshot.json").read_text())
    assert use_state4["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_debug_snapshot.json").read_bytes()).hexdigest() == use_state4["previous_manifest_sha256"]
    assert set(use_state4["files"]).issubset(reviewed_fortran)
    for name, record in use_state4["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_state4["files"]
    use_pack = json.loads((root / "docs/reference/rmcdhf_use_only_pack_snapshot.json").read_text())
    assert use_pack["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_state4_snapshot.json").read_bytes()).hexdigest() == use_pack["previous_manifest_sha256"]
    assert set(use_pack["files"]).issubset(reviewed_fortran)
    for name, record in use_pack["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_pack["files"]
    use_pos = json.loads((root / "docs/reference/rmcdhf_use_only_pos_snapshot.json").read_text())
    assert use_pos["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_pack_snapshot.json").read_bytes()).hexdigest() == use_pos["previous_manifest_sha256"]
    assert set(use_pos["files"]).issubset(reviewed_fortran)
    for name, record in use_pos["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_pos["files"]
    use_eigv = json.loads((root / "docs/reference/rmcdhf_use_only_eigv_snapshot.json").read_text())
    assert use_eigv["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_pos_snapshot.json").read_bytes()).hexdigest() == use_eigv["previous_manifest_sha256"]
    assert set(use_eigv["files"]).issubset(reviewed_fortran)
    for name, record in use_eigv["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_eigv["files"]
    use_radial = json.loads((root / "docs/reference/rmcdhf_use_only_radial_snapshot.json").read_text())
    assert use_radial["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_eigv_snapshot.json").read_bytes()).hexdigest() == use_radial["previous_manifest_sha256"]
    assert set(use_radial["files"]).issubset(reviewed_fortran)
    for name, record in use_radial["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_radial["files"]
    use_scf = json.loads((root / "docs/reference/rmcdhf_use_only_scf_snapshot.json").read_text())
    assert use_scf["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_radial_snapshot.json").read_bytes()).hexdigest() == use_scf["previous_manifest_sha256"]
    assert set(use_scf["files"]).issubset(reviewed_fortran)
    for name, record in use_scf["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_scf["files"]
    use_memory = json.loads((root / "docs/reference/rmcdhf_use_only_memory_snapshot.json").read_text())
    assert use_memory["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_scf_snapshot.json").read_bytes()).hexdigest() == use_memory["previous_manifest_sha256"]
    assert set(use_memory["files"]).issubset(reviewed_fortran)
    for name, record in use_memory["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_memory["files"]
    use_def_orb = json.loads((root / "docs/reference/rmcdhf_use_only_def_orb_snapshot.json").read_text())
    assert use_def_orb["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_memory_snapshot.json").read_bytes()).hexdigest() == use_def_orb["previous_manifest_sha256"]
    assert set(use_def_orb["files"]).issubset(reviewed_fortran)
    for name, record in use_def_orb["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_def_orb["files"]
    use_orb_label = json.loads((root / "docs/reference/rmcdhf_use_only_orb_label_snapshot.json").read_text())
    assert use_orb_label["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_def_orb_snapshot.json").read_bytes()).hexdigest() == use_orb_label["previous_manifest_sha256"]
    assert set(use_orb_label["files"]).issubset(reviewed_fortran)
    for name, record in use_orb_label["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_orb_label["files"]
    use_mpi = json.loads((root / "docs/reference/rmcdhf_use_only_mpi_snapshot.json").read_text())
    assert use_mpi["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_orb_label_snapshot.json").read_bytes()).hexdigest() == use_mpi["previous_manifest_sha256"]
    assert set(use_mpi["files"]).issubset(reviewed_fortran)
    for name, record in use_mpi["files"].items():
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"], name
    reviewed_fortran |= use_mpi["files"]
    migration_path = root / "docs/reference/rmcdhf_descriptor_snapshot.json"
    migration = json.loads(migration_path.read_text())
    assert migration["schema_version"] == 1
    assert hashlib.sha256((root / "docs/reference/rmcdhf_use_only_mpi_snapshot.json").read_bytes()).hexdigest() == migration["previous_manifest_sha256"]
    for name, record in migration["removed"].items():
        assert name in expected
        assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"]
        assert not (root / name).exists(), f"Retired text loader restored: {name}"
        del expected[name]
    for name, record in migration["files"].items():
        assert name in expected
        if name in reviewed_fortran:
            assert record["before_suite_sha256"] == reviewed_fortran[name]["after_suite_sha256"]
    reviewed_files = reviewed_fortran | migration["files"]
    csf_queries_path = root / "docs/reference/rmcdhf_csf_queries_snapshot.json"
    csf_queries = json.loads(csf_queries_path.read_text())
    assert csf_queries["schema_version"] == 1, "Unknown CSF queries snapshot schema"
    assert hashlib.sha256(migration_path.read_bytes()).hexdigest() == csf_queries["previous_manifest_sha256"], \
        "Descriptor snapshot changed; review CSF queries provenance"
    assert set(csf_queries["files"]).issubset(expected.keys() & reviewed_fortran.keys())
    for name, record in csf_queries["files"].items():
        assert record["before_suite_sha256"] == reviewed_files[name]["after_suite_sha256"], \
            f"CSF queries snapshot has a different predecessor: {name}"
    reviewed_files |= csf_queries["files"]
    kind_cleanup_path = root / "docs/reference/rmcdhf_kind_alias_cleanup_snapshot.json"
    kind_cleanup = json.loads(kind_cleanup_path.read_text())
    assert kind_cleanup["schema_version"] == 1, "Unknown kind alias cleanup snapshot schema"
    assert hashlib.sha256(csf_queries_path.read_bytes()).hexdigest() == kind_cleanup["previous_manifest_sha256"], \
        "CSF queries snapshot changed; review kind alias cleanup provenance"
    assert set(kind_cleanup["files"]).issubset(expected.keys() & reviewed_fortran.keys())
    for name, record in kind_cleanup["files"].items():
        assert record["before_suite_sha256"] == reviewed_files[name]["after_suite_sha256"], \
            f"Kind alias cleanup snapshot has a different predecessor: {name}"
    reviewed_files |= kind_cleanup["files"]
    davidson_naming_path = root / "docs/reference/rmcdhf_davidson_naming_snapshot.json"
    davidson_naming = json.loads(davidson_naming_path.read_text())
    assert davidson_naming["schema_version"] == 1, "Unknown Davidson naming snapshot schema"
    assert hashlib.sha256(kind_cleanup_path.read_bytes()).hexdigest() == davidson_naming["previous_manifest_sha256"], \
        "Kind alias cleanup snapshot changed; review Davidson naming provenance"
    assert set(davidson_naming["files"]).issubset(expected)
    for name, record in davidson_naming["files"].items():
        if name in reviewed_files:
            previous_digest = reviewed_files[name]["after_suite_sha256"]
        else:
            # These build files matched upstream byte for byte before renaming.
            previous_digest = expected[name]["sha256"]
        assert record["before_suite_sha256"] == previous_digest, \
            f"Davidson naming snapshot has a different predecessor: {name}"
    reviewed_files |= davidson_naming["files"]
    core_naming_path = root / "docs/reference/rmcdhf_core_names_snapshot.json"
    core_naming = json.loads(core_naming_path.read_text())
    assert core_naming["schema_version"] == 1, "Unknown GRASP core naming snapshot schema"
    assert hashlib.sha256(davidson_naming_path.read_bytes()).hexdigest() == core_naming["previous_manifest_sha256"], \
        "Davidson naming snapshot changed; review GRASP core naming provenance"
    assert set(core_naming["files"]).issubset(expected)
    for name, record in core_naming["files"].items():
        if name in reviewed_files:
            previous_digest = reviewed_files[name]["after_suite_sha256"]
        else:
            previous_digest = expected[name]["sha256"]
        assert record["before_suite_sha256"] == previous_digest, \
            f"GRASP core naming snapshot has a different predecessor: {name}"
    reviewed_files |= core_naming["files"]
    shared_naming_path = root / "docs/reference/rmcdhf_grasp_shared_naming_snapshot.json"
    shared_naming = json.loads(shared_naming_path.read_text())
    assert shared_naming["schema_version"] == 1, "Unknown GRASP shared naming snapshot schema"
    assert hashlib.sha256(core_naming_path.read_bytes()).hexdigest() == shared_naming["previous_manifest_sha256"], \
        "Core naming snapshot changed; review GRASP shared naming provenance"
    assert set(shared_naming["files"]).issubset(expected.keys() & reviewed_files.keys())
    for name, record in shared_naming["files"].items():
        assert record["before_suite_sha256"] == reviewed_files[name]["after_suite_sha256"], \
            f"GRASP shared naming snapshot has a different predecessor: {name}"
    reviewed_files |= shared_naming["files"]
    mpi_naming_path = root / "docs/reference/rmcdhf_grasp_mpi_naming_snapshot.json"
    mpi_naming = json.loads(mpi_naming_path.read_text())
    assert mpi_naming["schema_version"] == 1, "Unknown GRASP MPI naming snapshot schema"
    assert hashlib.sha256(shared_naming_path.read_bytes()).hexdigest() == mpi_naming["previous_manifest_sha256"], \
        "Shared naming snapshot changed; review GRASP MPI naming provenance"
    assert set(mpi_naming["files"]).issubset(expected.keys() & reviewed_files.keys())
    inline_path = root / "docs/reference/rmcdhf_inline_output_formats_snapshot.json"
    inline_formats = json.loads(inline_path.read_text())
    assert inline_formats["schema_version"] == 1, "Unknown inline output formats snapshot schema"
    assert hashlib.sha256(mpi_naming_path.read_bytes()).hexdigest() == inline_formats["previous_manifest_sha256"], \
        "GRASP MPI naming snapshot changed; review inline output format provenance"
    assert set(inline_formats["files"]).issubset(expected.keys() & reviewed_files.keys())
    for name, record in inline_formats["files"].items():
        assert record["before_suite_sha256"] == reviewed_files[name]["after_suite_sha256"], \
            f"Inline output format snapshot has a different predecessor: {name}"
    reviewed_files |= inline_formats["files"]
    io_units_path = root / "docs/reference/rmcdhf_io_units_snapshot.json"
    io_units = json.loads(io_units_path.read_text())
    assert io_units["schema_version"] == 1, "Unknown I/O unit snapshot schema"
    assert hashlib.sha256(inline_path.read_bytes()).hexdigest() == io_units["previous_manifest_sha256"], \
        "Inline output formats snapshot changed; review I/O unit provenance"
    assert set(io_units["files"]).issubset(expected.keys() & reviewed_files.keys())
    for name, record in io_units["files"].items():
        assert record["before_suite_sha256"] == reviewed_files[name]["after_suite_sha256"], \
            f"I/O unit snapshot has a different predecessor: {name}"
    state_files_path = root / "docs/reference/rmcdhf_state_files_snapshot.json"
    state_files = json.loads(state_files_path.read_text())
    assert state_files["schema_version"] == 1
    assert hashlib.sha256(io_units_path.read_bytes()).hexdigest() == state_files["previous_manifest_sha256"], \
        "I/O unit snapshot changed; review state-name provenance"
    for name, record in state_files["files"].items():
        assert record["before_suite_sha256"] == io_units["files"][name]["after_suite_sha256"], name
    reviewed_io = io_units["files"] | state_files["files"]
    app_names_path = root / "docs/reference/rmcdhf_app_variable_names_snapshot.json"
    app_names = json.loads(app_names_path.read_text())
    assert app_names["schema_version"] == 1, "Unknown app variable names snapshot schema"
    assert hashlib.sha256(state_files_path.read_bytes()).hexdigest() == app_names["previous_manifest_sha256"], \
        "State-name snapshot changed; review app variable naming provenance"
    previous = reviewed_files | reviewed_io
    assert set(app_names["files"]).issubset(name for name in previous if name.endswith(".f90"))
    for name, record in app_names["files"].items():
        assert record["before_suite_sha256"] == previous[name]["after_suite_sha256"], \
            f"App variable names snapshot has a different predecessor: {name}"
    reviewed_migration = reviewed_io | app_names["files"]
    app_symbols_path = root / APP_SYMBOLS_SNAPSHOT
    app_symbols = json.loads(app_symbols_path.read_text())
    assert app_symbols["schema_version"] == 1, "Unknown app symbols snapshot schema"
    assert hashlib.sha256(app_names_path.read_bytes()).hexdigest() == app_symbols["previous_manifest_sha256"], \
        "App variable names snapshot changed; review app symbol naming provenance"
    assert app_symbols["files"].keys() <= expected.keys()
    renamed_paths = {}
    previous = reviewed_files | reviewed_migration
    for name, record in app_symbols["files"].items():
        prior = previous[name]["after_suite_sha256"] if name in previous else expected[name]["sha256"]
        assert record["before_suite_sha256"] == prior, f"App symbols snapshot has a different predecessor: {name}"
        suite_path = Path(record["suite_path"])
        previous_path = Path(name.replace("src/appl/", "src/apps/", 1))
        assert suite_path.parent == previous_path.parent, f"Unexpected app sources directory: {name}"
        renamed_paths[name] = suite_path.as_posix()
    reviewed_migration |= app_symbols["files"]

    app_style_path = root / "docs/reference/rmcdhf_app_source_style_snapshot.json"
    app_style = json.loads(app_style_path.read_text())
    assert app_style["schema_version"] == 1, "Unknown app source style snapshot schema"
    assert hashlib.sha256(app_symbols_path.read_bytes()).hexdigest() == app_style["previous_manifest_sha256"], \
        "App symbols snapshot changed; review app source style provenance"
    assert set(app_style["files"]).issubset(
        name for name in app_symbols["files"] if name.endswith(".f90")
    )
    for name, record in app_style["files"].items():
        assert record["before_suite_sha256"] == app_symbols["files"][name]["after_suite_sha256"], \
            f"App source style snapshot has a different predecessor: {name}"
    reviewed_migration |= app_style["files"]

    app_modules = json.loads((root / APP_MODULES_SNAPSHOT).read_text())
    assert app_modules["schema_version"] == 1, "Unknown app modules snapshot schema"
    assert hashlib.sha256(app_style_path.read_bytes()).hexdigest() == app_modules["previous_manifest_sha256"], \
        "App source style snapshot changed; review app module provenance"
    merged_interfaces = app_modules["merged_interfaces"]
    assert merged_interfaces.keys() <= expected.keys()
    assert set(merged_interfaces.values()) <= app_modules["files"].keys()
    assert app_modules["files"].keys() <= expected.keys() - merged_interfaces.keys()
    for interface, implementation in merged_interfaces.items():
        former = Path(renamed_paths[interface])
        retained = Path(renamed_paths[implementation])
        assert former.name.endswith("_interface.f90") and former.parent == retained.parent, interface
    for name, record in app_modules["files"].items():
        assert record["before_suite_sha256"] == reviewed_migration[name]["after_suite_sha256"], \
            f"App modules snapshot has a different predecessor: {name}"
    reviewed_migration |= app_modules["files"]
    expected = {name: record for name, record in expected.items() if name not in merged_interfaces}

    semantic_names = json.loads((root / SEMANTIC_NAMES_SNAPSHOT).read_text())
    assert semantic_names["schema_version"] == 1, "Unknown semantic names snapshot schema"
    assert semantic_names["previous_manifest_sha256"] == hashlib.sha256((root / APP_MODULES_SNAPSHOT).read_bytes()).hexdigest(), \
        "App modules snapshot changed; review semantic naming provenance"
    previous = reviewed_files | reviewed_migration
    for name, record in semantic_names["files"].items():
        assert name in expected and name.endswith(".f90") and name in previous, name
        assert record["before_suite_sha256"] == previous[name]["after_suite_sha256"], \
            f"Semantic names snapshot has a different predecessor: {name}"
    reviewed_migration |= semantic_names["files"]

    audited_names = json.loads((root / AUDITED_NAMES_SNAPSHOT).read_text())
    assert audited_names["schema_version"] == 1, "Unknown audited names snapshot schema"
    assert audited_names["previous_manifest_sha256"] == hashlib.sha256((root / SEMANTIC_NAMES_SNAPSHOT).read_bytes()).hexdigest(), \
        "Semantic names snapshot changed; review audited naming provenance"
    previous = reviewed_files | reviewed_migration
    for name, record in audited_names["files"].items():
        assert name in expected and name in previous, name
        assert record["before_suite_sha256"] == previous[name]["after_suite_sha256"], \
            f"Audited names snapshot has a different predecessor: {name}"
        if "suite_path" in record:
            former = suite_source_path(root, name, renamed_paths)
            updated = root / record["suite_path"]
            assert former.parent == updated.parent and former.suffix == updated.suffix == ".f90", record["suite_path"]
            renamed_paths[name] = record["suite_path"]
    reviewed_migration |= audited_names["files"]

    mcp_snapshot = json.loads((root / MPI_MCP_IO_SNAPSHOT).read_text())
    added = {}
    predecessor = AUDITED_NAMES_SNAPSHOT
    for snapshot in (DIAGNOSTIC_LABELS_SNAPSHOT, V3_J_LABELS_SNAPSHOT, APP_COMMENTS_SNAPSHOT, MPI_RUN_LOGGING_SNAPSHOT,
                     MPI_MCP_IO_SNAPSHOT, RUN_TIMER_SNAPSHOT, REQUIRED_ASF_ROOTS_SNAPSHOT,
                     RADIAL_STREAM_SNAPSHOT, MIXING_STREAM_SNAPSHOT, MPI_F08_SNAPSHOT, SPARSE_STORAGE_SNAPSHOT,
                     BASIS_OWNERSHIP_SNAPSHOT, HOST_MEMORY_SNAPSHOT, STATE_INPUT_SNAPSHOT):
        app_update = json.loads((root / snapshot).read_text())
        assert app_update["schema_version"] == 1, f"Unknown application snapshot schema: {snapshot}"
        assert app_update["previous_manifest_sha256"] == hashlib.sha256(
            (root / predecessor).read_bytes()
        ).hexdigest(), f"Previous manifest changed: {predecessor}"
        previous = reviewed_files | reviewed_migration
        for name, record in app_update["files"].items():
            allowed_build_file = snapshot == MPI_MCP_IO_SNAPSHOT and Path(name).name in {
                "CMakeLists.txt", "Makefile", "BUILDCONF.sh",
            }
            assert (name in expected or name in added) and (name.endswith(".f90") or allowed_build_file), name
            if name in previous:
                previous_digest = previous[name]["after_suite_sha256"]
            elif name in added:
                previous_digest = added[name]
            else:
                previous_digest = expected[name]["sha256"]
            assert record["before_suite_sha256"] == previous_digest, name
        reviewed_migration |= app_update["files"]
        if snapshot == MPI_MCP_IO_SNAPSHOT:
            added = mcp_snapshot["added"]
        predecessor = snapshot

    for name, record in mcp_snapshot["removed"].items():
        previous = reviewed_files | reviewed_migration
        previous_digest = previous[name]["after_suite_sha256"] if name in previous else expected[name]["sha256"]
        assert record["before_suite_sha256"] == previous_digest, name
        path = suite_source_path(root, name, renamed_paths)
        assert not path.exists(), f"Retired MPI memory probe restored: {path}"
        del expected[name]

    actual = {
        path
        for name in ("rmcdhf", "rmcdhf_mem", "rmcdhf_mpi", "rmcdhf_mem_mpi")
        for path in (root / "src/apps" / name).iterdir()
        if path.suffix in (".f90", ".sh", ".txt") or path.name == "Makefile"
    }
    expected_paths = {name: suite_source_path(root, name, renamed_paths) for name in expected}
    paths = set(expected_paths.values())
    assert len(paths) == len(expected), "App source paths must be unique"
    assert paths.isdisjoint(root / name for name in added)
    paths.update(root / name for name in added)
    assert actual == paths, (actual - paths, paths - actual)
    for name, digest in added.items():
        reviewed = reviewed_migration.get(name, {}).get("after_suite_sha256", digest)
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == reviewed, f"Unreviewed new app source: {name}"
    for relative, record in expected.items():
        path = expected_paths[relative]
        if relative in reviewed_migration:
            assert reviewed_migration[relative]["after_suite_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest(), \
                f"Unreviewed deviation from app migration snapshot: {relative}"
            continue
        if relative in reviewed_files:
            source = path.read_bytes()
            reviewed = reviewed_files[relative]
            digests = {hashlib.sha256(source).hexdigest()}
            if relative in mpi_naming["files"]:
                # Check each renamed import or build reference, then compare the
                # original spelling with the reviewed source digest.
                before_mpi_names = reverse_reviewed_edits(source.decode(), path, mpi_naming["files"][relative])
                digests.add(hashlib.sha256(before_mpi_names.encode()).hexdigest())
            assert reviewed["after_suite_sha256"] in digests, \
                f"Unreviewed deviation from reviewed GRASP snapshot: {relative}"
            continue
        edits = compatibility.get(relative, {}).get("edits", [])
        digest = hashlib.sha256(original_bytes(path, edits)).hexdigest()
        assert digest == record["sha256"], f"Unreviewed deviation from GRASP source: {relative}"
    print(f"PASS: {len(expected)} imported files match reviewed snapshots "
          f"({len(expected.keys() & reviewed_files.keys())} exact suite file hashes)")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
