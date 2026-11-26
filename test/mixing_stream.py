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


def write_stream(
    path: Path, electrons: int, orbitals: int, blocks: list[Block], descriptor: Path | None = None
) -> None:
    """Write a native mixing fixture, optionally bound to its V3 descriptor."""
    digest = sha256(descriptor.read_bytes()).hexdigest() if descriptor else "0" * 64
    total_levels = sum(len(block.roots) for block in blocks)
    total_csfs = sum(block.csf_count for block in blocks)
    header = Struct("<8s5q64s")
    directory = Struct("<6q")
    level = Struct("<qd")
    offset = header.size + directory.size * len(blocks) + level.size * total_levels
    with path.open("wb") as stream:
        stream.write(header.pack(b"CSF--MIX", electrons, orbitals, len(blocks), total_csfs, total_levels, digest.encode()))
        row = 0
        for block in blocks:
            stream.write(directory.pack(block.two_j, block.parity, block.csf_count, len(block.roots), row, offset))
            row += len(block.roots)
            offset += len(block.coefficients)
        for block in blocks:
            for root, energy in zip(block.roots, block.energies, strict=True):
                stream.write(level.pack(root, energy))
        for block in blocks:
            stream.write(block.coefficients)


def read_stream(path: Path, descriptor: Path | None = None) -> list[Block]:
    raw = path.read_bytes()
    header = Struct("<8s5q64s")
    directory = Struct("<6q")
    level = Struct("<qd")
    magic, electrons, orbitals, blocks, csfs, roots, basis_id = header.unpack_from(raw)
    if magic != b"CSF--MIX" or min(electrons, orbitals, blocks, csfs, roots) <= 0:
        raise ValueError(f"{path}: invalid mixing stream header")
    if descriptor is not None and basis_id != sha256(descriptor.read_bytes()).hexdigest().encode():
        raise ValueError(f"{path}: mixing stream basis digest does not match {descriptor.name}")
    table = header.size + directory.size * blocks
    end = table + level.size * roots
    results = []
    total_csfs = total_roots = 0
    for index in range(blocks):
        two_j, parity, count, levels, row, offset = directory.unpack_from(raw, header.size + index * directory.size)
        if not (parity in (-1, 1) and two_j >= 0 and offset == end and row == total_roots and 0 < levels <= count):
            raise ValueError(f"{path}: invalid mixing block directory entry {index + 1}")
        root_ids, energies = [], []
        for level_index in range(row, row + levels):
            root, energy = level.unpack_from(raw, table + level_index * level.size)
            if not 1 <= root <= count:
                raise ValueError(f"{path}: invalid mixing root {root} in block {index + 1}")
            root_ids.append(root)
            energies.append(energy)
        end = offset + 8 * count * levels
        results.append(Block(two_j, parity, count, root_ids, energies, raw[offset:end]))
        total_csfs += count
        total_roots += levels
    if (total_csfs, total_roots, end) != (csfs, roots, len(raw)):
        raise ValueError(f"{path}: mixing stream dimensions do not match its payload")
    return results


def mixing_blocks(path: Path) -> list[Block]:
    with path.open("rb") as stream:
        signature = stream.read(8)
    if signature == b"CSF--MIX":
        return read_stream(path)
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
