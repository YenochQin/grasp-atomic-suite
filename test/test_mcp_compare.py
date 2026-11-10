"""Independent binary fixtures and CLI regression tests for the MCP comparator."""

from __future__ import annotations

import argparse
import io
import json
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import mcp_compare as mcp
except ModuleNotFoundError as exc:
    if exc.name and exc.name.split(".")[0] == "numpy":
        print("SKIP: MCP parser tests require NumPy", file=sys.stderr)
        raise SystemExit(77) from exc
    raise

KEY = 215
T = KEY + 2
V = ((KEY + 1) * KEY + 2) * KEY + 2
# (block, unit, label, column CSF, row CSF, value)
TERMS = [
    (1, 31, T, 2, 1, 0.5),
    (1, 31, T, 2, 1, 0.25),
    (1, 31, T, 3, 2, -0.75),
    (1, 32, V, 3, 1, 0.125),
    (1, 33, V, 3, 3, -1.25),
    (2, 32, V, 1, 1, 0.375),
]


def fixture(
    root: Path,
    *,
    ranks: int = 1,
    endian: str = "<",
    marker: int = 4,
    integer: int = 4,
    reverse: bool = False,
    split: bool = False,
    terms: list | None = None,
    omit_pair: tuple | None = None,
    bad_index: bool = False,
    bad_column: bool = False,
    bad_count: bool = False,
    bad_endpoint: bool = False,
) -> None:
    root.mkdir()
    (root / "rcsf.inp").write_text("same CSFs and orbitals\n", encoding="ascii")
    terms = TERMS if terms is None else terms
    intcode = "i" if integer == 4 else "q"
    markcode = "i" if marker == 4 else "q"

    def ints(*items: int) -> bytes:
        return struct.pack(endian + intcode * len(items), *items)

    def record(payload: bytes, continuation: bool = False) -> bytes:
        parts = [payload]
        if split and continuation and len(payload) > 1:
            # Deliberately split inside a number, as the Fortran runtime may do.
            parts = [payload[:5], payload[5:]]
        out = b""
        for i, part in enumerate(parts):
            leading = -len(part) if i < len(parts) - 1 else len(part)
            trailing = -len(part) if i > 0 else len(part)
            out += (
                struct.pack(endian + markcode, leading)
                + part
                + struct.pack(endian + markcode, trailing)
            )
        return out

    for rank in range(ranks):
        stem = "mcp" if ranks == 1 else f"mcp{rank:03d}"
        common = (
            record(b"MCP" + ints(2, rank, ranks))
            + record(ints(2, 4, 4))
            + record(ints(0, 0, 0))
        )
        structure = (
            record(ints(0, 2, 1))
            + record(ints(3, 1))
            + record(b"   0+      1+   ")
            + common
        )
        outputs = {unit: common for unit in (31, 32, 33)}
        for block, ncf in ((1, 3), (2, 1)):
            pairs, ends = [], []
            for column in range(rank + 1, ncf + 1, ranks):
                rows = list(range(1, column))
                if reverse:
                    rows.reverse()
                rows.append(column)
                pairs.extend(
                    (column, row) for row in rows if (block, column, row) != omit_pair
                )
                ends.append(len(pairs))
            pair_index = {pair: i for i, pair in enumerate(pairs, 1)}
            structure += record(b"MCP" + ints(block, ncf)) + record(ints(len(pairs)))
            stored_ends = list(ends)
            if bad_endpoint and block == 1:
                stored_ends[-1] += 1
            structure += record(ints(*stored_ends, *(row for _, row in pairs)), True)
            for unit in (31, 32, 33):
                selected = [
                    t
                    for t in terms
                    if t[0] == block and t[1] == unit and (t[3] - 1) % ranks == rank
                ]
                outputs[unit] += record(
                    b"MCP" + ints(block, ncf, len(selected) + int(bad_count))
                )
                for label in sorted({t[2] for t in selected}):
                    values = [t for t in selected if t[2] == label]
                    if reverse:
                        values.reverse()
                    outputs[unit] += record(ints(label, len(values)))
                    payload = b""
                    for _, _, _, col, row, value in values:
                        index = 0 if bad_index else pair_index[(col, row)]
                        payload += ints(1 if bad_column else col, index) + struct.pack(
                            endian + "d", value
                        )
                    outputs[unit] += record(payload, True)
                outputs[unit] += record(ints(0, 0))
        (root / f"{stem}.30").write_bytes(structure)
        for unit, content in outputs.items():
            (root / f"{stem}.{unit}").write_bytes(content)


class MCPTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.left, self.right = self.root / "reference", self.root / "candidate"
        fixture(self.left)

    def run_cli(self, *args: str) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = mcp.main(list(args))
        return code, out.getvalue(), err.getvalue()

    def compare(self, **kwargs) -> dict:
        return mcp.compare(mcp.Dataset(self.left), mcp.Dataset(self.right), **kwargs)

    def test_formats_mpi_and_permutation(self) -> None:
        for endian in ("<", ">"):
            for marker in (4, 8):
                for integer in (4, 8):
                    with self.subTest(endian=endian, marker=marker, integer=integer):
                        fixture(
                            self.right,
                            ranks=2,
                            endian=endian,
                            marker=marker,
                            integer=integer,
                            reverse=True,
                            split=True,
                        )
                        report = self.compare()
                        self.assertEqual(report["status"], "equal")
                        self.assertEqual(
                            sum(s["compared"] for s in report["channels"]), len(TERMS)
                        )
                        shutil.rmtree(self.right)

    def test_empty_mpi_rank(self) -> None:
        fixture(self.right, ranks=4)
        self.assertEqual(self.compare()["status"], "equal")

    def test_difference_location_and_tolerance(self) -> None:
        terms = list(TERMS)
        terms[3] = (*terms[3][:-1], terms[3][-1] + 1e-9)
        fixture(self.right, terms=terms)
        report = self.compare()
        self.assertEqual(report["difference_count"], 1)
        d = report["differences"][0]
        self.assertEqual(
            (d["unit"], d["block"], d["column_csf"], d["row_csf"]), (32, 1, 3, 1)
        )
        self.assertEqual(d["orbitals"], (1, 2, 1, 2))
        self.assertEqual(self.compare(atol=2e-9)["status"], "equal")
        self.assertEqual(self.compare(rtol=2e-8)["status"], "equal")
        self.assertEqual(self.compare(max_diffs=0)["differences"], [])

    def test_duplicate_multiset_and_optional_sum(self) -> None:
        fixture(self.right, terms=[(1, 31, T, 2, 1, 0.75), *TERMS[2:]])
        self.assertEqual(self.compare()["status"], "different")
        self.assertEqual(self.compare(sum_duplicates=True)["status"], "equal")

    def test_missing_group_and_empty_channel(self) -> None:
        fixture(self.right, terms=[t for t in TERMS if t[1] != 31])
        report = self.compare()
        self.assertEqual(report["status"], "different")
        self.assertEqual(report["channels"][0]["missing"], 3)

    def test_missing_contribution(self) -> None:
        fixture(self.right, terms=TERMS[1:])
        self.assertEqual(self.compare()["channels"][0]["missing"], 1)

    def test_structure_difference_even_without_coefficients(self) -> None:
        fixture(self.right, omit_pair=(1, 2, 1), terms=[t for t in TERMS if t[1] != 31])
        self.assertTrue(
            any(d["type"] == "sparse_structure" for d in self.compare()["differences"])
        )

    def test_bad_indices_counts_endpoints_and_values(self) -> None:
        cases = [
            {"bad_index": True},
            {"bad_column": True},
            {"bad_count": True},
            {"bad_endpoint": True},
        ]
        cases += [
            {"terms": [(*TERMS[0][:-1], value), *TERMS[1:]]}
            for value in (float("nan"), float("inf"))
        ]
        for kwargs in cases:
            with self.subTest(kwargs=kwargs):
                fixture(self.right, **kwargs)
                with self.assertRaises(mcp.MCPError):
                    mcp.inspect(mcp.Dataset(self.right), KEY)
                shutil.rmtree(self.right)

    def test_truncation_and_trailing_bytes(self) -> None:
        fixture(self.right)
        for unit in range(30, 34):
            path = self.right / f"mcp.{unit}"
            original = path.read_bytes()
            for bad in (original[:-1], original + b"extra"):
                path.write_bytes(bad)
                with self.assertRaises(mcp.MCPError):
                    mcp.inspect(mcp.Dataset(self.right), KEY)
            path.write_bytes(original)

    def test_continuation_sign_corruption(self) -> None:
        path = self.left / "mcp.31"
        data = bytearray(path.read_bytes())
        struct.pack_into("<i", data, 19, -15)
        path.write_bytes(data)
        with self.assertRaises(mcp.MCPError):
            mcp.inspect(mcp.Dataset(self.left), KEY)

    def test_file_completeness(self) -> None:
        fixture(self.right, ranks=2)
        (self.right / "mcp001.33").unlink()
        with self.assertRaisesRegex(mcp.MCPError, "file set"):
            mcp.Dataset(self.right)
        for p in self.right.glob("mcp001.*"):
            p.unlink()
        with self.assertRaisesRegex(mcp.MCPError, "MPI ranks"):
            mcp.Dataset(self.right)

    def test_extra_and_mixed_files(self) -> None:
        fixture(self.right)
        extra = self.right / "mcp.34"
        extra.write_bytes(b"stale")
        with self.assertRaisesRegex(mcp.MCPError, "file set"):
            mcp.Dataset(self.right)
        extra.unlink()
        shutil.copyfile(self.right / "mcp.30", self.right / "mcp000.30")
        with self.assertRaisesRegex(mcp.MCPError, "mixed"):
            mcp.Dataset(self.right)

    def test_identity_is_checked(self) -> None:
        fixture(self.right)
        code, out, _ = self.run_cli("compare", str(self.left), str(self.right))
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(out)["csf_identity"]["status"], "sha256_identical")
        (self.right / "rcsf.inp").unlink()
        self.assertEqual(self.run_cli("compare", str(self.left), str(self.right))[0], 2)
        self.assertEqual(
            self.run_cli(
                "compare", str(self.left), str(self.right), "--assume-same-csfs"
            )[0],
            0,
        )
        (self.right / "rcsf.inp").write_text("changed", encoding="ascii")
        self.assertEqual(
            self.run_cli(
                "compare", str(self.left), str(self.right), "--assume-same-csfs"
            )[0],
            2,
        )

    def test_cli_reports_and_dump(self) -> None:
        fixture(self.right, terms=TERMS[:-1])
        report = self.root / "report.json"
        self.assertEqual(
            self.run_cli(
                "compare", str(self.left), str(self.right), "--json", str(report)
            )[0],
            1,
        )
        self.assertEqual(json.loads(report.read_text())["status"], "different")
        code, out, _ = self.run_cli(
            "dump", str(self.left), "--unit", "31", "--limit", "1"
        )
        self.assertEqual(code, 0)
        self.assertEqual(len(out.splitlines()), 1)
        self.assertEqual(json.loads(out)["coefficient"], 0.5)
        self.assertEqual(
            self.run_cli(
                "inspect", str(self.left), "--json", str(self.left / "mcp.31")
            )[0],
            2,
        )
        self.assertEqual(self.run_cli("inspect", str(self.left), "--keyorb", "2")[0], 2)

    def test_metadata_mismatch(self) -> None:
        fixture(self.right)
        for path in self.right.glob("mcp.*"):
            data = path.read_bytes().replace(
                struct.pack("<iii", 2, 4, 4), struct.pack("<iii", 3, 4, 4), 1
            )
            path.write_bytes(data)
        self.assertEqual(self.compare()["differences"][0]["type"], "metadata")
        path = self.right / "mcp.33"
        path.write_bytes(path.read_bytes()[:-1])
        with self.assertRaises(mcp.MCPError):
            self.compare()

    def test_tolerance_validation(self) -> None:
        for value in ("nan", "inf", "-1"):
            with self.assertRaises(argparse.ArgumentTypeError):
                mcp.nonnegative(value)

    def test_extreme_finite_coefficients(self) -> None:
        fixture(self.right, terms=[(*t[:-1], -1e308) for t in TERMS])
        shutil.rmtree(self.left)
        fixture(self.left, terms=[(*t[:-1], 1e308) for t in TERMS])
        report = self.compare()
        self.assertEqual(report["difference_count"], len(TERMS))
        self.assertEqual(report["channels"][0]["max_rel_error"], 2.0)
        self.assertIsInstance(report["channels"][0]["max_abs_error"], str)
        self.assertEqual(self.compare(rtol=2)["status"], "equal")
        json.dumps(report, allow_nan=False)

    def test_fortran_generated_records(self) -> None:
        compiler = shutil.which("gfortran-15")
        if compiler is None:
            self.skipTest("GNU Fortran 15 unavailable")
        source = self.root / "writer.f90"
        source.write_text(
            """program writer
  implicit none
  integer :: u
  character(2) :: suffix
  do u=30,33
    write(suffix,'(i2)') u
    open(u,file='mcp.'//suffix,form='unformatted',status='replace')
    if(u==30) then
      write(u) 0,1,1
      write(u) 2
      write(u) '   0+   '
    endif
    write(u) 'MCP',1,0,1
    write(u) 2,2,2
    write(u) .false.,0,.false.
    if(u==30) then
      write(u) 'MCP',1,2
      write(u) 3
      write(u) 1,3,1,1,2
    else if(u==31) then
      write(u) 'MCP',1,2,1
      write(u) 217,1
      write(u) 2,2,0.5d0
      write(u) 0,0
    else
      write(u) 'MCP',1,2,0
      write(u) 0,0
    endif
    close(u)
  enddo
end program
""",
            encoding="ascii",
        )
        # GNU runtime creates real chained subrecords for the 20-byte payload.
        subprocess.run(
            [
                compiler,
                "-std=f2018",
                "-fmax-subrecord-length=16",
                str(source),
                "-o",
                str(self.root / "writer"),
            ],
            check=True,
            capture_output=True,
        )
        self.right.mkdir()
        subprocess.run(
            [str(self.root / "writer")], cwd=self.right, check=True, capture_output=True
        )
        dataset = mcp.Dataset(self.right)
        groups = list(dataset.groups(31, KEY))
        self.assertEqual(groups[0][3].tolist(), [0.5])
        self.assertEqual(groups[0][2].tolist(), [7])
        self.assertEqual(mcp.inspect(dataset, KEY)["status"], "valid")


if __name__ == "__main__":
    unittest.main()
