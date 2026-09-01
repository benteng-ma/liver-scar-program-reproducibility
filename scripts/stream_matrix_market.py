from __future__ import annotations

from pathlib import Path

import numpy as np


def read_header(path: Path) -> tuple[int, int, int, int]:
    """Return n_rows, n_cols, n_nonzero and byte offset of coordinate data."""
    with path.open("rb") as handle:
        banner = handle.readline().decode("ascii").strip()
        if banner != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"unsupported Matrix Market banner: {banner}")
        while True:
            line = handle.readline()
            if not line:
                raise ValueError("missing Matrix Market dimension line")
            if not line.startswith(b"%"):
                n_rows, n_cols, n_nonzero = map(int, line.split())
                return n_rows, n_cols, n_nonzero, handle.tell()


def iter_coordinate_blocks(
    path: Path,
    offset: int,
    block_bytes: int = 64 * 1024 * 1024,
):
    """Yield zero-based row, zero-based column and integer value arrays."""
    with path.open("rb") as handle:
        handle.seek(offset)
        remainder = b""
        while True:
            block = handle.read(block_bytes)
            if not block:
                payload = remainder
                remainder = b""
            else:
                payload = remainder + block
                boundary = payload.rfind(b"\n")
                if boundary < 0:
                    remainder = payload
                    continue
                remainder = payload[boundary + 1 :]
                payload = payload[: boundary + 1]
            if payload:
                values = np.fromstring(payload.decode("ascii"), sep=" ", dtype=np.int64)
                if values.size % 3:
                    raise ValueError(
                        f"coordinate block has {values.size} integers, not a multiple of 3"
                    )
                yield values[0::3] - 1, values[1::3] - 1, values[2::3]
            if not block:
                break
        if remainder:
            raise ValueError("unparsed Matrix Market trailing bytes")


def aggregate_by_cell_group(
    path: Path,
    cell_groups: np.ndarray,
    n_groups: int,
    block_bytes: int = 64 * 1024 * 1024,
) -> tuple[np.ndarray, int]:
    """Aggregate a gene-by-cell Matrix Market file to gene-by-group counts."""
    n_genes, n_cells, n_nonzero, offset = read_header(path)
    if cell_groups.shape != (n_cells,):
        raise ValueError(
            f"cell group length {cell_groups.size} does not match matrix columns {n_cells}"
        )
    if cell_groups.min() < 0 or cell_groups.max() >= n_groups:
        raise ValueError("cell group codes are outside the declared group range")
    flat = np.zeros(n_genes * n_groups, dtype=np.int64)
    observed_nonzero = 0
    for rows, cols, values in iter_coordinate_blocks(path, offset, block_bytes):
        if rows.min() < 0 or rows.max() >= n_genes:
            raise ValueError("gene index out of bounds")
        if cols.min() < 0 or cols.max() >= n_cells:
            raise ValueError("cell index out of bounds")
        index = rows * n_groups + cell_groups[cols]
        flat += np.bincount(index, weights=values, minlength=flat.size).astype(np.int64)
        observed_nonzero += values.size
    if observed_nonzero != n_nonzero:
        raise ValueError(
            f"observed {observed_nonzero} coordinates, header declares {n_nonzero}"
        )
    return flat.reshape(n_genes, n_groups), observed_nonzero
