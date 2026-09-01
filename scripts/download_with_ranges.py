from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import time
from pathlib import Path

import requests


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resumable verified HTTP range downloader")
    parser.add_argument("url")
    parser.add_argument("output", type=Path)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--chunk-mib", type=int, default=4)
    parser.add_argument("--workers", type=int, default=24)
    parser.add_argument("--retries", type=int, default=8)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def download_part(
    url: str,
    part_path: Path,
    start: int,
    end: int,
    retries: int,
) -> tuple[int, int]:
    expected = end - start + 1
    if part_path.exists() and part_path.stat().st_size == expected:
        return start, expected
    for attempt in range(1, retries + 1):
        temporary = part_path.with_suffix(part_path.suffix + ".tmp")
        try:
            with requests.get(
                url,
                headers={"Range": f"bytes={start}-{end}"},
                stream=True,
                timeout=(30, 600),
            ) as response:
                response.raise_for_status()
                if response.status_code != 206:
                    raise RuntimeError(f"range request returned HTTP {response.status_code}")
                content_range = response.headers.get("Content-Range", "")
                if not content_range.startswith(f"bytes {start}-{end}/"):
                    raise RuntimeError(f"unexpected Content-Range: {content_range}")
                with temporary.open("wb") as handle:
                    for block in response.iter_content(chunk_size=1024 * 1024):
                        if block:
                            handle.write(block)
            if temporary.stat().st_size != expected:
                raise RuntimeError(
                    f"part length {temporary.stat().st_size}, expected {expected}"
                )
            os.replace(temporary, part_path)
            return start, expected
        except Exception:
            if temporary.exists():
                temporary.unlink()
            if attempt == retries:
                raise
            time.sleep(min(60, 2**attempt))
    raise AssertionError("unreachable")


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    part_dir = args.output.with_name(args.output.name + ".parts")
    part_dir.mkdir(parents=True, exist_ok=True)
    chunk = args.chunk_mib * 1024 * 1024
    ranges = []
    for start in range(0, args.size, chunk):
        end = min(args.size - 1, start + chunk - 1)
        part = part_dir / f"{start:012d}-{end:012d}.part"
        ranges.append((args.url, part, start, end, args.retries))

    completed_bytes = sum(
        path.stat().st_size
        for _, path, start, end, _ in ranges
        if path.exists() and path.stat().st_size == end - start + 1
    )
    print(
        f"parts={len(ranges)} workers={args.workers} "
        f"resume_bytes={completed_bytes} total_bytes={args.size}",
        flush=True,
    )
    started = time.time()
    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(download_part, *item) for item in ranges]
        for future in concurrent.futures.as_completed(futures):
            _, part_bytes = future.result()
            completed_bytes += part_bytes
            done += 1
            elapsed = max(time.time() - started, 0.001)
            print(
                f"completed_parts={done}/{len(ranges)} "
                f"downloaded_or_verified={completed_bytes}/{args.size} "
                f"rate_mib_s={completed_bytes / elapsed / 1024 / 1024:.3f}",
                flush=True,
            )

    assembly = args.output.with_suffix(args.output.suffix + ".assembling")
    with assembly.open("wb") as target:
        for _, part, _, _, _ in sorted(ranges, key=lambda item: item[2]):
            with part.open("rb") as source:
                for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
                    target.write(block)
    if assembly.stat().st_size != args.size:
        raise RuntimeError(
            f"assembled length {assembly.stat().st_size}, expected {args.size}"
        )
    os.replace(assembly, args.output)
    manifest = {
        "url": args.url,
        "output": str(args.output),
        "size": args.size,
        "sha256": sha256(args.output),
        "chunk_bytes": chunk,
        "parts": len(ranges),
        "workers": args.workers,
        "range_validation": "HTTP 206 and exact Content-Range/part length",
    }
    manifest_path = args.output.with_suffix(args.output.suffix + ".download.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
