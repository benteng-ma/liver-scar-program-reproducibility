from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stream one gzip layer to disk and record byte-level provenance"
    )
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--block-mib", type=int, default=8)
    args = parser.parse_args()

    source = args.source.resolve(strict=True)
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    block_bytes = args.block_mib * 1024 * 1024
    source_digest = hashlib.sha256()
    output_digest = hashlib.sha256()
    output_bytes = 0

    with source.open("rb") as raw:
        for block in iter(lambda: raw.read(block_bytes), b""):
            source_digest.update(block)

    with gzip.open(source, "rb") as compressed, output.open("wb") as unwrapped:
        for block in iter(lambda: compressed.read(block_bytes), b""):
            unwrapped.write(block)
            output_digest.update(block)
            output_bytes += len(block)

    manifest = {
        "source": str(source),
        "source_bytes": source.stat().st_size,
        "source_sha256": source_digest.hexdigest().upper(),
        "output": str(output),
        "output_bytes": output_bytes,
        "output_sha256": output_digest.hexdigest().upper(),
        "operation": "one-layer gzip decompression",
    }
    manifest_path = output.with_suffix(output.suffix + ".unwrap.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
