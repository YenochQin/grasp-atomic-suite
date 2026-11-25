"""Independent reader for the on-disk Suite mixing layout used by regressions."""

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from struct import Struct, unpack_from


@dataclass
class Block:
    two_j: int
    parity: int
    csf_count: int
    roots: list[int]
    energies: list[float]
    coefficients: bytes


def read_stream(path: Path, descriptor: Path | None = None, *, signature: bytes = b"CSF--MIX") -> list[Block]:
    raw = path.read_bytes()
    header = Struct("<8s5q64s")
    directory = Struct("<6q")
    level = Struct("<qd")
    magic, electrons, orbitals, blocks, csfs, roots, basis_id = header.unpack_from(raw)
    assert magic == signature and min(electrons, orbitals, blocks, csfs, roots) > 0
    if descriptor is not None:
        assert basis_id == sha256(descriptor.read_bytes()).hexdigest().encode()
    table = header.size + directory.size * blocks
    end = table + level.size * roots
    results = []
    total_csfs = total_roots = 0
    for index in range(blocks):
        two_j, parity, count, levels, row, offset = directory.unpack_from(raw, header.size + index * directory.size)
        assert parity in (-1, 1) and two_j >= 0 and offset == end and row == total_roots and 0 < levels <= count
        root_ids, energies = [], []
        for level_index in range(row, row + levels):
            root, energy = level.unpack_from(raw, table + level_index * level.size)
            assert 1 <= root <= count
            root_ids.append(root)
            energies.append(energy)
        end = offset + 8 * count * levels
        results.append(Block(two_j, parity, count, root_ids, energies, raw[offset:end]))
        total_csfs += count
        total_roots += levels
    assert (total_csfs, total_roots, end) == (csfs, roots, len(raw))
    return results


def mixing_blocks(path: Path) -> list[Block]:
    with path.open("rb") as stream:
        signature = stream.read(8)
    if signature in (b"CSF--MIX", b"GRMIX001"):
        return read_stream(path, signature=signature)
    from rdensity.regression import records

    data = records(path)
    assert data[0] == b"G92MIX"
    results = []
    for index in range(2, len(data), 4):
        _, csfs, levels, j_code, parity = unpack_from("=5i", data[index])
        roots = list(unpack_from(f"={levels}i", data[index + 1]))
        shift, *energies = unpack_from(f"={levels+1}d", data[index + 2])
        results.append(Block(j_code - 1, parity, csfs, roots, [shift + energy for energy in energies], data[index + 3]))
    return results
