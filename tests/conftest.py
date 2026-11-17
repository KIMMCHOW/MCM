"""Pytest configuration: put src/ (and analysis/) on sys.path and share cached model runs."""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT / "src", ROOT / "analysis"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import pytest  # noqa: E402


@lru_cache(maxsize=None)
def lsoda(prim: str, Tamb_C: float = 25.0, q_loss: float = 0.0, s0: float = 1.0, **kw):
    """Cached LSODA run of a constant primitive (fresh or at SOH = 1 − q_loss, r_growth = γ_R q_loss)."""
    from phonebatt import RP63, PRIMITIVES, aged, run
    p = aged(RP63, q_loss) if q_loss else RP63
    if kw:
        p = p.replace(**kw)
    return run(PRIMITIVES[prim], Tamb_C, s0=s0, p=p)


@pytest.fixture(scope="session")
def run_const():
    return lsoda
