"""Verify imported sources against reviewed upstream and lowercase snapshots.

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


def original_bytes(path: Path, edits: list[dict[str, str | int]] | None = None) -> bytes:
    text = path.read_text()
    for edit in reversed(edits or []):
        before, after, count = edit["before"], edit["after"], edit["count"]
        assert isinstance(before, str) and isinstance(after, str) and isinstance(count, int)
        pattern = "(?m)^" + re.escape(after)
        assert len(re.findall(pattern, text)) == count, f"Reviewed F2018 edit missing or ambiguous: {path}"
        text = re.sub(pattern, lambda _: before, text)
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
    lowercase = json.loads((root / "docs/reference/rmcdhf_lowercase_snapshot.json").read_text())
    assert lowercase["schema_version"] == 1, "Unknown lowercase snapshot schema"
    assert hashlib.sha256(upstream_path.read_bytes()).hexdigest() == lowercase["upstream_manifest_sha256"], \
        "Upstream manifest changed; review lowercase provenance"
    assert hashlib.sha256(compatibility_path.read_bytes()).hexdigest() == lowercase["f2018_compatibility_sha256"], \
        "F2018 edits changed; review lowercase provenance"
    assert set(lowercase["files"]) == {name for name in expected if name.endswith(".f90")}, \
        "Lowercase snapshot must cover exactly the imported Fortran files"
    actual = {
        str(path.relative_to(root))
        for name in ("rmcdhf", "rmcdhf_mem", "rmcdhf_mpi", "rmcdhf_mem_mpi")
        for path in (root / "src/appl" / name).iterdir()
        if path.suffix in (".f90", ".sh", ".txt") or path.name == "Makefile"
    }
    assert actual == set(expected), (actual - set(expected), set(expected) - actual)
    for relative, record in expected.items():
        if relative in lowercase["files"]:
            digest = hashlib.sha256((root / relative).read_bytes()).hexdigest()
            assert digest == lowercase["files"][relative]["after_suite_sha256"], \
                f"Unreviewed deviation from lowercase GRASP snapshot: {relative}"
            continue
        edits = compatibility.get(relative, {}).get("edits", [])
        digest = hashlib.sha256(original_bytes(root / relative, edits)).hexdigest()
        assert digest == record["sha256"], f"Unreviewed deviation from GRASP source: {relative}"
    print(f"PASS: {len(expected)} imported files match reviewed snapshots "
          f"({len(lowercase['files'])} exact lowercase Fortran identities)")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
