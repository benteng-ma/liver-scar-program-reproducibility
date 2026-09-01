"""Verify the isolated Phase 0 Python runtime and write an audit log."""

from __future__ import annotations

import platform
from pathlib import Path

import h5py
import numpy
import openpyxl
import pandas
import scipy


root = Path(__file__).resolve().parents[1]
lines = [
    f"python {platform.python_version()}",
    f"numpy {numpy.__version__}",
    f"pandas {pandas.__version__}",
    f"scipy {scipy.__version__}",
    f"h5py {h5py.__version__}",
    f"openpyxl {openpyxl.__version__}",
    "phase0_disease_effects FALSE",
    "environment_scope PHASE0_METADATA_AND_SMOKE_ONLY",
]
target = root / "results" / "logs" / "environment_verification.txt"
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))

