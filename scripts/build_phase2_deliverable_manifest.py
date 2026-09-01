from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    roots = [
        repo / "config",
        repo / "literature",
        repo / "manuscript",
        repo / "metadata",
        repo / "reports",
        repo / "results",
        repo / "scripts",
        repo / "tests",
    ]
    output = repo / "results" / "logs" / "phase2_deliverable_manifest.csv"
    excluded = {
        output,
        repo / "metadata" / "file_manifest.csv",
    }
    rows: list[dict[str, object]] = []
    top_level_files = [
        repo / "STATUS.md",
        repo / "requirements-phase0.txt",
        repo / "requirements-phase2.txt",
    ]
    for path in top_level_files:
        rows.append(
            {
                "path": path.relative_to(repo).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    for root in roots:
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            if path in excluded or "__pycache__" in path.parts:
                continue
            rows.append(
                {
                    "path": path.relative_to(repo).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    pd.DataFrame(rows).to_csv(output, index=False)
    print(f"Wrote {len(rows)} Phase 2 artifact hashes to {output}")


if __name__ == "__main__":
    main()
