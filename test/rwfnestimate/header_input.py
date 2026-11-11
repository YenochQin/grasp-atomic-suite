"""Exercise the native TOML reader against the ordinary text input path."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from high_l import header, validate

BINARY = Path()
FIXTURES = Path(__file__).parent / "fixtures"


def encode_header(lines: list[str], electrons: object = 28) -> str:
    return (
        "[header_info]\nheader_lines = " + json.dumps(lines) + "\n"
        + "electron_count = " + json.dumps(electrons) + "\n"
    )


class HeaderInput(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rwfn-header-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.counter = 0
        self.csf, self.expected = header(28)
        self.lines = self.csf.splitlines()[:5]

    def run_case(
        self, *, toml: str | None = None, csf: str | None = None,
        stdin: str = "y\n2\n*\n", arguments: list[str] | None = None,
        filename: str = "header.toml",
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        self.counter += 1
        directory = self.root / str(self.counter)
        directory.mkdir()
        shutil.copyfile(FIXTURES / "Ni.isodata", directory / "isodata")
        command = [str(BINARY)]
        original = None
        if toml is not None:
            path = directory / filename
            path.write_text(toml, encoding="utf-8")
            original = hashlib.sha256(path.read_bytes()).hexdigest()
            command += ["--header", filename]
        if csf is not None:
            (directory / "rcsf.inp").write_text(csf)
        if arguments is not None:
            command += arguments
        result = subprocess.run(command, cwd=directory, input=stdin, text=True,
                                capture_output=True, timeout=30, check=False)
        (directory / "stdout.log").write_text(result.stdout)
        (directory / "stderr.log").write_text(result.stderr)
        if original is not None:
            self.assertEqual(original, hashlib.sha256(path.read_bytes()).hexdigest())
        return result, directory

    def assert_success(self, result: subprocess.CompletedProcess[str]) -> None:
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("RWFNESTIMATE: Execution complete.", result.stderr)

    def assert_rejected(self, toml: str, message: str) -> None:
        result, directory = self.run_case(toml=toml)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(message, result.stderr)
        self.assertFalse((directory / "rwfn.inp").exists())

    def test_text_and_header_are_byte_identical_for_all_estimation_modes(self) -> None:
        for mode, stdin in [("tf", "y\n2\n*\n"), ("hyd", "y\n3\n*\n"),
                            ("custom", "y\n4\n*\n2\n")]:
            with self.subTest(mode=mode):
                old, text_dir = self.run_case(csf=self.csf, stdin=stdin)
                new, toml_dir = self.run_case(toml=encode_header(self.lines), stdin=stdin)
                self.assert_success(old)
                self.assert_success(new)
                self.assertIn("core subshells=7, orbitals=88, total electrons=28", new.stdout)
                self.assertEqual((text_dir / "rwfn.inp").read_bytes(), (toml_dir / "rwfn.inp").read_bytes())
                self.assertFalse((toml_dir / "rcsf.inp").exists())
                if mode == "tf":
                    validate(toml_dir / "rwfn.inp", self.expected, 28)

    def test_no_core_and_closed_core_only(self) -> None:
        no_core = self.lines.copy()
        no_core[1] = ""
        no_core[3] = self.lines[1] + " " + self.lines[3]
        occupation = "  1s ( 2)  2s ( 2)  2p-( 2)  2p ( 4)  3s ( 2)  3p-( 2)  3p ( 4)"
        cases = [(no_core, 28, occupation + self.csf.splitlines()[5]),
                 (["Core subshells:", "  1s", "Peel subshells:", "", "CSF(s):"], 2, "")]
        for lines, electrons, first in cases:
            with self.subTest(electrons=electrons):
                text, text_dir = self.run_case(csf="\n".join([*lines, first]) + "\n")
                toml, toml_dir = self.run_case(toml=encode_header(lines, electrons))
                self.assert_success(text)
                self.assert_success(toml)
                self.assertEqual((text_dir / "rwfn.inp").read_bytes(), (toml_dir / "rwfn.inp").read_bytes())

    def test_toml_syntax_and_long_paths(self) -> None:
        baseline, baseline_dir = self.run_case(csf=self.csf)
        self.assert_success(baseline)
        long_lines = self.lines.copy()
        long_lines[3] = (" " * 80).join(long_lines[3].split())
        multiline = "[header_info]\n# total includes closed core\nelectron_count=2_8\nheader_lines=[\n"
        multiline += ", # orbital metadata\n".join(json.dumps(line) for line in long_lines)
        multiline += ",\n]\n[extra]\nlabel='中文 metadata'\n"
        escaped = encode_header(self.lines).replace("Core subshells:", r"Core\u0020subshells:")
        inline = "header_info = {electron_count=28, header_lines=" + json.dumps(self.lines) + "}\n"
        literal = "[header_info]\nelectron_count=28\nheader_lines=[" + ",".join("'" + s + "'" for s in self.lines) + "]\n"
        for encoding in [multiline, escaped, inline, literal]:
            with self.subTest(encoding=encoding[:50]):
                result, directory = self.run_case(toml=encoding, filename="输入 with spaces " + "x" * 140 + ".toml")
                self.assert_success(result)
                self.assertEqual((baseline_dir / "rwfn.inp").read_bytes(), (directory / "rwfn.inp").read_bytes())

    def test_explicit_header_does_not_read_csf_text(self) -> None:
        result, _ = self.run_case(toml=encode_header(self.lines), csf="deliberately invalid text\n")
        self.assert_success(result)

    def test_missing_and_invalid_electron_counts(self) -> None:
        valid = encode_header(self.lines)
        self.assert_rejected(valid.replace("electron_count = 28\n", ""), "regenerate old headers with rCSFs")
        for value in [True, "28", 28.0, [28]]:
            with self.subTest(value=value):
                self.assert_rejected(encode_header(self.lines, value), "electron_count must be an integer")
        for value in [0, -1, 2147483648]:
            with self.subTest(value=value):
                self.assert_rejected(encode_header(self.lines, value), "electron_count must be positive")
        self.assert_rejected(encode_header(self.lines, 17), "smaller than the closed-core population")
        self.assert_rejected(encode_header(self.lines, 9999), "exceeds the orbital table capacity")

    def test_bad_toml_and_schema(self) -> None:
        self.assert_rejected("[header_info\n", "cannot parse TOML header")
        self.assert_rejected(encode_header(self.lines) + "electron_count=29\n", "cannot parse TOML header")
        self.assert_rejected("[something_else]\nx=1\n", "[header_info]")
        self.assert_rejected("[header_info]\nelectron_count=28\n", "header_lines must be an array")
        self.assert_rejected(encode_header(self.lines[:4]), "five strings")
        wrong_type = self.lines.copy()
        wrong_type[2] = 3
        self.assert_rejected(encode_header(wrong_type), "strings only")
        for position in [0, 2, 4]:
            lines = self.lines.copy()
            lines[position] = "incorrect label"
            self.assert_rejected(encode_header(lines), "header label")
        lines = self.lines.copy()
        lines[3] += "\n  11s"
        self.assert_rejected(encode_header(lines), "must not contain newline")

    def test_invalid_orbital_tables(self) -> None:
        for label, diagnostic in [("1s", "duplicate orbital"), ("0s", "must exceed l"),
                                  ("1p", "must exceed l"), ("4z", "unsupported orbital letter"),
                                  ("4s-", "invalid relativistic suffix"), ("6h--", "invalid subshell suffix"),
                                  ("p", "invalid subshell label"), ("12", "invalid subshell label")]:
            with self.subTest(label=label):
                lines = self.lines.copy()
                lines[3] += " " + label
                self.assert_rejected(encode_header(lines), diagnostic)
        empty = ["Core subshells:", "", "Peel subshells:", "", "CSF(s):"]
        self.assert_rejected(encode_header(empty), "empty orbital table")
        empty[3] = " ".join(f"{n}s" for n in range(1, 129))
        self.assert_rejected(encode_header(empty), "NNNW capacity")

    def test_cli_errors_and_help(self) -> None:
        help_result, directory = self.run_case(arguments=["--help"], stdin="")
        self.assertEqual(help_result.returncode, 0)
        self.assertIn("--header", help_result.stdout)
        self.assertFalse((directory / "rwfn.inp").exists())
        for arguments, diagnostic in [(["--header"], "missing command argument"),
                                      (["--unknown"], "unknown option"),
                                      (["--csf-text", "a", "--header", "b"], "select only one"),
                                      (["--header", "a", "--header", "b"], "select only one"),
                                      (["--header", "missing.toml"], "cannot parse TOML header"),
                                      (["--header", "rwfn.inp"], "reserved output name")]:
            with self.subTest(arguments=arguments):
                result, _ = self.run_case(arguments=arguments)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(diagnostic, result.stderr)

    def test_explicit_text_filename(self) -> None:
        result, _ = self.run_case(csf=self.csf, arguments=["--csf-text", "rcsf.inp"])
        self.assert_success(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    options, remaining = parser.parse_known_args()
    BINARY = options.binary.resolve()
    unittest.main(argv=[__file__, *remaining])
