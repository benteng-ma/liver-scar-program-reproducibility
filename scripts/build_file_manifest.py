"""Create an auditable hash manifest for local Phase 0 source and smoke files."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
paths = list((ROOT / "data" / "raw" / "phase0_metadata").rglob("*"))
paths += list((ROOT / "data" / "external").rglob("*"))


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


rows = []
for path in sorted(p for p in paths if p.is_file()):
    rel = path.relative_to(ROOT).as_posix()
    if "_supplementary/" in rel and path.suffix.lower() in {".jpg", ".gif", ".tif"}:
        continue
    rows.append({
        "relative_path": rel,
        "bytes": path.stat().st_size,
        "sha256": digest(path),
        "acquired_or_generated": "2026-08-30",
        "phase0_role": "smoke object" if "/external/" in rel else "source metadata/supplement",
        "tracked_by_git": "no; immutable source data are represented by this manifest",
    })

target = ROOT / "metadata" / "file_manifest.csv"
with target.open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
print(f"files={len(rows)} bytes={sum(int(row['bytes']) for row in rows)}")

