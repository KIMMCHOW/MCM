"""Shared helpers of the 'basics' analysis scripts (AN0, AN1, AN2, AN4, AN16).

* ``meta(...)``: the top-level ``"meta"`` block of every results JSON (FINAL_SPEC §14): spec version, git SHA,
  whether ``src/`` differs from that commit, seed, solver and step.
* ``write_json(name, data)``: results/<name> with NumPy scalars/arrays converted and non-finite floats -> null.
* ``write_tex(name, text)``: paper/tables/<name>.
* ``load_proto()``: the reference prototype module ``docs/design/scripts/final_proto.py`` (regression targets).
"""
from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
TABLES = ROOT / "paper" / "tables"
SPEC = "FINAL_SPEC 2026-11-17"


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def git_info() -> dict:
    try:
        sha = _git("rev-parse", "HEAD")
        dirty = bool(_git("status", "--porcelain", "--", "src"))
    except Exception:  # not a git checkout
        sha, dirty = "unknown", None
    return dict(git_sha=sha, git_src_dirty=dirty)


def meta(script: str, seed=None, solver=None, h_s=None, **extra) -> dict:
    """Top-level meta block (FINAL_SPEC §14)."""
    import phonebatt
    d = dict(spec=SPEC, seed=seed, solver=solver, h_s=h_s, script=f"analysis/{script}",
             phonebatt_version=phonebatt.__version__, generated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                                                 time.gmtime()))
    d.update(git_info())
    d.update(extra)
    return d


def jsonable(o):
    """Recursively convert NumPy objects to JSON types; NaN/inf -> None."""
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        f = float(o)
        return f if math.isfinite(f) else None
    return o


def write_json(name: str, data: dict) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    out = RESULTS / name
    with open(out, "w") as f:
        json.dump(jsonable(data), f, indent=1, allow_nan=False)
        f.write("\n")
    return out


def write_tex(name: str, text: str) -> Path:
    TABLES.mkdir(parents=True, exist_ok=True)
    out = TABLES / name
    out.write_text(text)
    return out


def load_proto():
    """Import docs/design/scripts/final_proto.py (the [PROTO] reference implementation) as a module."""
    path = ROOT / "docs" / "design" / "scripts" / "final_proto.py"
    spec = importlib.util.spec_from_file_location("final_proto", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["final_proto"] = mod
    spec.loader.exec_module(mod)
    return mod


def f(x) -> float:
    """float() of a 0-d array / NumPy scalar."""
    return float(np.asarray(x))


def bib_keys() -> set[str]:
    """Entry keys of paper/refs.bib (used to check that every \\cite in a generated table resolves)."""
    import re
    txt = (ROOT / "paper" / "refs.bib").read_text()
    return set(re.findall(r"^@\w+\{([^,]+),", txt, flags=re.M))
