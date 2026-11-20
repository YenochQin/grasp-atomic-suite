"""Check indentation and label rules in the suite's free-form Fortran."""

import argparse
from collections.abc import Iterator
from pathlib import Path
import re


LABEL = re.compile(r"^(?:[0-9]{1,5}(?=\s|$)|[a-z][a-z0-9_]*\s*:(?!:))", re.I)
LABEL_BRANCH = re.compile(r"(?:^|\))\s*go\s*to(?=\s|\()", re.I)
LABELED_DO = re.compile(r"^(?:do|assign)\s+[0-9]{1,5}(?=\s|,|$)", re.I)
IO_STATEMENT = re.compile(r"(?:^|\))\s*(?:read|write|open|close|wait|rewind|backspace|endfile)\s*\(", re.I)
IO_LABEL = re.compile(r"\b(?:err|end|eor|fmt)\s*=\s*[0-9]{1,5}(?=\s*[,)]|$)", re.I)


def statements(source: str) -> Iterator[tuple[int, str]]:
    """Yield statements with comments and character literals masked."""
    parts: list[str] = []
    quote = ""
    start = 0
    for number, raw in enumerate(source.splitlines(), 1):
        line = raw.lstrip()
        if line.startswith(("!$ ", "!$\t")):
            line = line[2:].lstrip()
        elif not line or line.startswith(("!", "#")):
            continue
        if parts:
            line = line.removeprefix("&")
        else:
            start = number

        code: list[str] = []
        for char in line:
            if quote:
                if char == quote:
                    quote = ""
            elif char in "'\"":
                quote = char
                code.append(" ")
            elif char == "!":
                break
            else:
                code.append(char)

        fragment = "".join(code).rstrip()
        continued = bool(quote) or fragment.endswith("&")
        parts.append(fragment.removesuffix("&"))
        if not continued:
            for statement in "".join(parts).split(";"):
                if statement.strip():
                    yield start, statement.strip()
            parts.clear()


def violations(source: str) -> Iterator[tuple[int, str]]:
    for number, line in enumerate(source.splitlines(), 1):
        indent = line[: len(line) - len(line.lstrip())]
        if indent.strip(" ") or len(indent) % 4:
            yield number, "Indent with groups of four spaces"
    for number, statement in statements(source):
        if LABEL.match(statement):
            yield number, "Use an unlabeled statement or construct"
        elif LABEL_BRANCH.search(statement):
            yield number, "Use structured control flow instead of GOTO"
        elif LABELED_DO.match(statement) or (IO_STATEMENT.search(statement) and IO_LABEL.search(statement)):
            yield number, "Use an unlabeled DO or I/O statement"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", choices=("apps", "lib"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    directory = root / "src" / args.directory
    sources = sorted(path for path in directory.rglob("*") if path.suffix.lower() == ".f90")
    if not sources:
        raise SystemExit(f"No Fortran sources found in {directory}")
    errors = [
        f"{path.relative_to(root)}:{line}: {message}"
        for path in sources for line, message in violations(path.read_text())
    ]
    if errors:
        raise SystemExit("\n".join(errors))
    print(f"Checked indentation and label rules in {len(sources)} Fortran sources under src/{args.directory}")


if __name__ == "__main__":
    main()
