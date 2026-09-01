from __future__ import annotations

import hashlib
import importlib
import json
import platform
import sys
from pathlib import Path


EXPECTED = {
    "anndata": "0.12.7",
    "h5py": "3.15.1",
    "matplotlib": "3.10.8",
    "numpy": "2.5.1",
    "openpyxl": "3.1.5",
    "pandas": "2.3.3",
    "yaml": "6.0.3",
    "requests": "2.32.5",
    "sklearn": "1.8.0",
    "scipy": "1.16.3",
    "seaborn": "0.13.2",
    "statsmodels": "0.14.6",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    cfg = Path(sys.prefix) / "pyvenv.cfg"
    cfg_text = cfg.read_text(encoding="utf-8") if cfg.exists() else ""
    isolated = "include-system-site-packages = false" in cfg_text.lower()
    observed = {}
    for module_name, expected_version in EXPECTED.items():
        module = importlib.import_module(module_name)
        observed_version = str(module.__version__)
        observed[module_name] = observed_version
        if observed_version != expected_version:
            raise RuntimeError(
                f"{module_name}: expected {expected_version}, observed {observed_version}"
            )
    if not isolated:
        raise RuntimeError("Phase 2 environment exposes system site packages")

    payload = {
        "verified": True,
        "implementation": "method-equivalent isolated Python",
        "python": sys.version,
        "executable": str(Path(sys.executable).resolve()),
        "executable_sha256": sha256(Path(sys.executable)),
        "platform": platform.platform(),
        "isolated_from_system_site_packages": isolated,
        "packages": observed,
        "requirements_sha256": sha256(repo / "requirements-phase2.txt"),
        "frozen_random_seed": 20260830,
        "frozen_random_modules_per_program_dataset_state": 1000,
    }
    output = repo / "results" / "logs" / "phase2_environment_verification.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
