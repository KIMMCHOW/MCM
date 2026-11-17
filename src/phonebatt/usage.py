"""Usage inputs u(t): activity primitives, schedules, scenarios S1–S16 and the usage CTMC (FINAL_SPEC §1.4, §5).

``UsageState`` is one constant input (§1.4).  ``Schedule`` is a cyclic list of ``(duration_s, UsageState)``.
``CTMC`` is the hub-and-spoke × connectivity usage chain of §5.3; its rates are in h⁻¹.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Callable, Iterable, Sequence

import numpy as np

from .params import RP63, Params, aged


# ------------------------------------------------------------------------------------------- usage state
@dataclass(frozen=True, eq=False)
class UsageState:
    """Piecewise-constant input u (§1.4).  Numeric fields may be arrays over the sample axis (jitter)."""

    name: str = ""
    s_scr: float = 0.0  # screen on
    awake: float = 0.0  # wakelock flag a_flag; effective a = max(s_scr, awake)
    b: float = 0.4  # brightness slider
    APL: float = 0.6  # average picture level
    m_APL: float = 1.0  # dark-theme multiplier
    f_r: float = 60.0  # refresh rate (Hz)
    x_c: float = 0.0  # SoC throughput demand
    w: float = 0.0  # Wi-Fi associated
    c: str = "none"  # cellular camping RAT: 'none' | 'lte' | 'nr'
    rho: str = "none"  # data path: 'wifi' | 'cell' | 'none'
    kappa: str = "burst"  # traffic class: 'burst' | 'rt'
    lam: float = 0.0  # burst / wake-up rate (1/s)
    d_dl: float = 0.0  # Mbit/s
    d_ul: float = 0.0  # Mbit/s
    RSRP: float = -85.0  # dBm
    RSSI: float | None = None  # dBm; Wi-Fi cliff multiplier (S11b) when set
    g: float = 0.0  # GNSS duty
    c_r: float = 0.0  # rear camera
    c_f: float = 0.0  # front camera
    au: float = 0.0  # audio
    bt: float = 0.0  # Bluetooth audio
    P_sens: float = 0.0  # W
    dPb: float = 0.0  # burst-envelope amplitude before m_b (W)
    P_extra: float = 0.0  # W, explicit extra load (tests, structural variants)
    T_amb_C: float | None = None  # ambient override for this segment (°C)

    def with_(self, **kw) -> "UsageState":
        return replace(self, **kw)

    @property
    def a(self):
        return np.maximum(self.s_scr, self.awake)


def _prim() -> dict[str, UsageState]:
    """Activity primitives of FINAL_SPEC §5.1."""
    U: dict[str, UsageState] = {}
    sb = dict(s_scr=0, awake=0, kappa="burst", lam=1 / 120, d_dl=0.0, dPb=1.0, c="lte")
    U["SB-W"] = UsageState("SB-W", **sb, w=1, rho="wifi")
    U["SB-L"] = UsageState("SB-L", **sb, w=0, rho="cell")
    U["SB-5"] = UsageState("SB-5", **dict(sb, lam=1 / 60, c="nr"), w=0, rho="cell")
    U["RD"] = UsageState("RD", s_scr=1, awake=1, b=0.4, APL=0.65, f_r=60, x_c=0.05, w=1, c="lte", rho="wifi",
                         lam=1 / 30, d_dl=0.05, dPb=3.0)
    U["SO"] = UsageState("SO", s_scr=1, awake=1, b=0.4, APL=0.55, f_r=120, x_c=0.25, w=1, c="lte", rho="wifi",
                         lam=2.0, d_dl=2.0, dPb=3.0)
    U["SO-L"] = replace(U["SO"], name="SO-L", w=0, rho="cell")
    U["VW"] = UsageState("VW", s_scr=1, awake=1, b=0.4, APL=0.30, f_r=60, x_c=0.10, w=1, c="lte", rho="wifi",
                         lam=0.25, d_dl=3.0, au=1, dPb=2.0)
    U["VL"] = replace(U["VW"], name="VL", w=0, rho="cell")
    U["V5"] = replace(U["VW"], name="V5", w=0, c="nr", rho="cell")
    U["NV"] = UsageState("NV", s_scr=1, awake=1, b=0.7, APL=0.55, f_r=60, x_c=0.25, w=0, c="lte", rho="cell",
                         lam=0.2, d_dl=0.1, g=1, P_sens=0.047, dPb=3.0)
    U["GM"] = UsageState("GM", s_scr=1, awake=1, b=0.5, APL=0.40, f_r=120, x_c=0.80, w=1, c="lte", rho="wifi",
                         lam=20.0, d_dl=0.2, d_ul=0.05, au=1, dPb=4.0)
    U["GM-L"] = replace(U["GM"], name="GM-L", w=0, rho="cell")
    U["VC"] = UsageState("VC", s_scr=1, awake=1, b=0.5, APL=0.45, f_r=60, x_c=0.30, w=0, c="lte", rho="cell",
                         kappa="rt", c_f=1, au=1, dPb=3.0)
    U["VC-W"] = replace(U["VC"], name="VC-W", w=1, rho="wifi")
    U["CM"] = UsageState("CM", s_scr=1, awake=1, b=0.6, APL=0.45, f_r=60, x_c=0.05, w=1, c="lte", rho="wifi",
                         lam=0.0, c_r=1, dPb=4.0)
    U["AU"] = UsageState("AU", s_scr=0, awake=1, w=1, c="lte", rho="wifi", lam=1 / 30, d_dl=0.256, au=1, bt=1,
                         dPb=1.0)
    U["AU-L"] = replace(U["AU"], name="AU-L", w=0, rho="cell")
    U["TRK"] = UsageState("TRK", s_scr=0, awake=1, w=0, c="lte", rho="cell", lam=1 / 60, d_dl=0.01, g=1,
                          P_sens=0.047, dPb=1.0)
    U["BU"] = UsageState("BU", s_scr=1, awake=1, b=1.0, APL=0.45, f_r=60, x_c=0.50, w=0, c="lte", rho="cell",
                         lam=1.0, d_ul=4.0, c_r=1, dPb=0.0)
    return U


PRIMITIVES: dict[str, UsageState] = _prim()

# §8 V1 Bonn configuration: airplane mode, OLED on and black, x_c = 0, no bursts.
BONN_IDLE = UsageState("BONN-idle", s_scr=1, b=0.4, APL=0.0, f_r=60, x_c=0.0, w=0, c="none", rho="none", dPb=0.0)
BONN_CASES = {"idle": BONN_IDLE, "cam": replace(BONN_IDLE, name="BONN-cam", c_r=1),
              "gps": replace(BONN_IDLE, name="BONN-gps", g=1),
              "cxl": replace(BONN_IDLE, name="BONN-cxl", c_r=1, g=1)}


# ------------------------------------------------------------------------------------------- schedules
class Schedule:
    """Cyclic deterministic schedule: a list of ``(duration_s, UsageState)`` repeated until the TTE event."""

    def __init__(self, segments: Sequence[tuple[float, UsageState]], name: str = ""):
        if not segments:
            raise ValueError("empty schedule")
        self.segments = [(float(d), u) for d, u in segments]
        self.name = name
        self.durations = np.array([d for d, _ in self.segments])
        self.ends = np.cumsum(self.durations)
        self.period = float(self.ends[-1])

    def __len__(self):
        return len(self.segments)

    def __iter__(self):
        return iter(self.segments)

    def __getitem__(self, k):
        return self.segments[k]

    @classmethod
    def constant(cls, u: UsageState | str, dt_s: float = 3600.0) -> "Schedule":
        u = PRIMITIVES[u] if isinstance(u, str) else u
        return cls([(dt_s, u)], name=u.name)

    @classmethod
    def from_minutes(cls, items: Iterable[tuple[float, str | UsageState]], name: str = "") -> "Schedule":
        return cls([(60.0 * m, PRIMITIVES[n] if isinstance(n, str) else n) for m, n in items], name=name)

    def map(self, fn: Callable[[UsageState], UsageState]) -> "Schedule":
        """Apply ``fn`` to every segment input (levers, structural variants)."""
        return Schedule([(d, fn(u)) for d, u in self.segments], name=self.name)

    def index_at(self, t_s):
        """Segment index active at time(s) t (cyclic; the segment [start, end) convention)."""
        tm = np.mod(np.asarray(t_s, dtype=float), self.period)
        return np.minimum(np.searchsorted(self.ends, tm, side="right"), len(self.segments) - 1)

    def mean(self, fn: Callable[[UsageState], float]) -> float:
        """Time average of ``fn(u)`` over one period."""
        return float(sum(d * fn(u) for d, u in self.segments) / self.period)


def s12_items() -> list[tuple[int, str]]:
    """S12 scheduled full day from 07:00 (minutes, primitive), FINAL_SPEC §5.2."""
    s = [(30, "RD"), (45, "NV")] + [(26, "SB-W"), (4, "SO")] * 7 + [(30, "VW"), (30, "SO")] + \
        [(26, "SB-W"), (4, "SO")] * 8 + [(45, "NV"), (30, "VC"), (60, "GM"), (120, "VW"), (60, "SO"), (60, "RD"),
                                         (480, "SB-W")]
    assert sum(m for m, _ in s) == 1440
    return s


def s12_schedule() -> Schedule:
    return Schedule.from_minutes(s12_items(), name="S12")


def s9_schedule(Tamb_C: float | None = None) -> Schedule:
    """S9 cold hike: periodic 600 s, TRK 570 s then BU 30 s."""
    trk, bu = PRIMITIVES["TRK"], PRIMITIVES["BU"]
    if Tamb_C is not None:
        trk, bu = replace(trk, T_amb_C=Tamb_C), replace(bu, T_amb_C=Tamb_C)
    return Schedule([(570.0, trk), (30.0, bu)], name="S9")


def wifi_cliff(RSSI):
    """Wi-Fi weak-signal multiplier m_W(RSSI) = 1 + 1.13·exp((−85 − RSSI)/2.54 dB) (S11b, DI13 [DERIVED])."""
    return 1.0 + 1.13 * np.exp((-85.0 - np.asarray(RSSI, dtype=float)) / 2.54)


# ------------------------------------------------------------------------------------------- scenarios (§5.2)
@dataclass
class Scenario:
    id: str
    schedule: "Schedule | CTMC"
    Tamb_C: float = 25.0
    T0_C: float | None = None
    s0: float = 1.0
    p: Params = field(default_factory=lambda: RP63)
    note: str = ""

    def run(self, **kw):
        from .model import run
        if isinstance(self.schedule, CTMC):
            raise TypeError("stochastic scenario: use model.run_batch")
        return run(self.schedule, self.Tamb_C, s0=self.s0, T0_C=self.T0_C, p=self.p, **kw)


_CONST = {"S1": "SB-W", "S2a": "SB-L", "S2b": "SB-5", "S3": "RD", "S4": "SO", "S5a": "VW", "S5b": "VL",
          "S5c": "V5", "S6": "NV", "S7": "GM", "S8a": "VC", "S8b": "VC-W", "S16": "AU", "CM": "CM"}


def scenario(id: str, Tamb_C: float | None = None, s0: float = 1.0, T0_C: float | None = None,
             p: Params = RP63, q_loss: float | None = None, soh: float | None = None, **opts) -> Scenario:
    """Scenario S1–S16 of FINAL_SPEC §5.2.

    Common options: ``Tamb_C`` (default 25 °C), ``s0`` (displayed initial SOC), ``T0_C`` (default = T_amb),
    ``p`` (base parameters), ``q_loss`` or ``soh`` (aged pack, r_growth = γ_R q_loss).
    Scenario-specific options: S7 ``throttling=False``; S10 ``base='S4'|'S12'``; S11 ``variant='a'|'b'|'c'``
    with ``RSRP`` (dBm) or ``RSSI`` (dBm, S11b); S13/S14 CTMC keyword options (``rat='nr'``, ...);
    S15 ``base`` in {S4, S5a, S7, S12, S13}.  Constant scenarios use a 3600 s repeated segment.
    """
    if soh is not None:
        q_loss = 1.0 - soh
    if q_loss is not None:
        p = aged(p, q_loss)
    Ta = 25.0 if Tamb_C is None else Tamb_C
    if id in _CONST:
        if id == "S7" and opts.get("throttling") is False:
            p = replace(p, throttling=False, burst=False)
        return Scenario(id, Schedule.constant(_CONST[id]), Ta, T0_C, s0, p)
    if id == "S9":
        T0 = 20.0 if T0_C is None else T0_C
        return Scenario(id, s9_schedule(Ta), Ta, T0, s0, p, note="T0 = 20 °C (carried out from indoors)")
    if id == "S10":
        base = opts.get("base", "S4")
        sc = scenario(base, Tamb_C=Ta, s0=s0, T0_C=T0_C, p=p)
        sc.id = f"S10({base})"
        return sc
    if id == "S11":
        var = opts.get("variant", "a")
        if var == "a":
            u = replace(PRIMITIVES["SO-L"], RSRP=float(opts.get("RSRP", -85.0)))
        elif var == "c":
            u = replace(PRIMITIVES["VC"], RSRP=float(opts.get("RSRP", -85.0)))
        elif var == "b":
            u = replace(PRIMITIVES["SO"], RSSI=float(opts.get("RSSI", -85.0)))
        else:
            raise ValueError(var)
        return Scenario(f"S11{var}", Schedule.constant(u), Ta, T0_C, s0, p)
    if id == "S12":
        return Scenario(id, s12_schedule(), Ta, T0_C, s0, p)
    if id in ("S13", "S14"):
        prof = "typical" if id == "S13" else "heavy"
        return Scenario(id, CTMC(prof, **opts), Ta, T0_C, s0, p)
    if id == "S15":
        return scenario(opts.pop("base", "S4"), Tamb_C=Ta, s0=s0, T0_C=T0_C, p=p, **opts)
    raise KeyError(id)


SCENARIO_IDS = ("S1", "S2a", "S2b", "S3", "S4", "S5a", "S5b", "S5c", "S6", "S7", "S8a", "S8b", "S9", "S10",
                "S11", "S12", "S13", "S14", "S15", "S16")


# ------------------------------------------------------------------------------------------- CTMC (§5.3)
ACTS = ("idle", "social", "video", "game", "nav", "vcall", "camera", "audio")
MAP_WIFI = dict(idle="SB-W", social="SO", video="VW", game="GM", nav="NV", vcall="VC-W", camera="CM", audio="AU")
MAP_LTE = dict(idle="SB-L", social="SO-L", video="VL", game="GM-L", nav="NV", vcall="VC", camera="CM", audio="AU-L")
PROFILES = {
    # (lam_on [1/h], p_j, d_j [min])
    "typical": (4.0, dict(social=.55, video=.12, game=.05, nav=.04, vcall=.04, camera=.08, audio=.12),
                dict(social=5, video=15, game=20, nav=20, vcall=10, camera=2, audio=30)),
    "heavy": (4.0, dict(social=.45, video=.20, game=.15, nav=.05, vcall=.05, camera=.05, audio=.05),
              dict(social=8, video=25, game=30, nav=25, vcall=15, camera=3, audio=30)),
}


@dataclass
class Paths:
    """Sampled CTMC paths: per path, jump times (s, increasing, padded with +inf) and the state after each jump.

    ``times[k, 0] = 0`` and ``states[k, 0] = x0``.
    """

    times: np.ndarray  # (n, M) seconds
    states: np.ndarray  # (n, M) int
    h: float

    @property
    def n(self) -> int:
        return self.times.shape[0]

    def state_at(self, t_s) -> np.ndarray:
        """State active at time t (scalar or per-path array) for every path."""
        t = np.broadcast_to(np.asarray(t_s, dtype=float), (self.n,))
        idx = (self.times <= t[:, None]).sum(axis=1) - 1
        return self.states[np.arange(self.n), np.maximum(idx, 0)]

    def state_for_step(self, k) -> np.ndarray:
        """State held during RK4 step k (scalar or per-path array): a jump at time t takes effect in step
        floor(t/h), i.e. the state after all jumps with t < (k+1) h (prototype convention)."""
        t = (np.asarray(k, dtype=float) + 1.0) * self.h
        t = np.broadcast_to(t, (self.n,))
        idx = (self.times < t[:, None]).sum(axis=1) - 1
        return self.states[np.arange(self.n), np.maximum(idx, 0)]

    def cursor(self) -> "PathCursor":
        return PathCursor(self)


class PathCursor:
    """Amortised O(n) per-step state lookup for monotonically increasing step indices."""

    def __init__(self, paths: Paths):
        self.paths = paths
        self.ptr = np.zeros(paths.n, dtype=np.int64)
        self.rows = np.arange(paths.n)
        self.M = paths.times.shape[1]

    def state_for_step(self, k: int) -> np.ndarray:
        t = (k + 1.0) * self.paths.h
        tm = self.paths.times
        while True:
            nxt = np.minimum(self.ptr + 1, self.M - 1)
            adv = (self.ptr + 1 < self.M) & (tm[self.rows, nxt] < t)
            if not adv.any():
                break
            self.ptr = np.where(adv, self.ptr + 1, self.ptr)
        return self.paths.states[self.rows, self.ptr]


class CTMC:
    """Usage CTMC: hub-and-spoke activity chain × Wi-Fi/cellular connectivity chain (§5.3).

    Parameters
    ----------
    profile : 'typical' | 'heavy' or a tuple ``(lam_on [1/h], {act: p}, {act: d_min})``.
    mean_w, mean_c : mean connectivity sojourns on Wi-Fi / cellular [h].
    rat : 'lte' (default) or 'nr' (5G variant: cellular states use NR primitives).
    lam_on, p, d_min : optional overrides of the profile values.
    """

    def __init__(self, profile="typical", mean_w: float = 2.5, mean_c: float = 2.5, rat: str = "lte",
                 lam_on: float | None = None, p: dict | None = None, d_min: dict | None = None,
                 state_map: Sequence[str] | None = None):
        if isinstance(profile, str):
            lam0, p0, d0 = PROFILES[profile]
            self.profile = profile
        else:
            lam0, p0, d0 = profile
            self.profile = "custom"
        self.lam_on = float(lam0 if lam_on is None else lam_on)
        self.p = dict(p0 if p is None else p)
        self.d_min = dict(d0 if d_min is None else d_min)
        self.mean_w, self.mean_c, self.rat = float(mean_w), float(mean_c), rat
        self.acts = ACTS
        if state_map is None:
            state_map = [m[a] for a in ACTS for m in (MAP_WIFI, MAP_LTE)]
        self.state_map = list(state_map)  # primitive name per product state (activity-major, Wi-Fi first)
        self.states = [(a, c) for a in ACTS for c in ("wifi", "cell")]

    # -------------------------------------------------------------------- generators (rates in 1/h)
    def activity_generator(self) -> np.ndarray:
        n = len(ACTS)
        L = np.zeros((n, n))
        for j, a in enumerate(ACTS[1:], 1):
            L[0, j] = self.lam_on * self.p[a]
            L[j, 0] = 60.0 / self.d_min[a]
        np.fill_diagonal(L, -L.sum(1))
        return L

    def connectivity_generator(self) -> np.ndarray:
        return np.array([[-1 / self.mean_w, 1 / self.mean_w], [1 / self.mean_c, -1 / self.mean_c]])

    def generator(self) -> np.ndarray:
        """Product generator Λ = Λ_A ⊗ I₂ + I₈ ⊗ Λ_C (16 states, 1/h)."""
        La, Lc = self.activity_generator(), self.connectivity_generator()
        return np.kron(La, np.eye(2)) + np.kron(np.eye(len(ACTS)), Lc)

    def primitives(self) -> list[UsageState]:
        """Inputs of the 16 product states.  With ``rat='nr'`` every LTE camping/data RAT becomes 5G NSA."""
        us = [PRIMITIVES[n] for n in self.state_map]
        if self.rat == "nr":
            us = [replace(u, c="nr") if u.c == "lte" else u for u in us]
        return us

    def state_powers(self, p: Params = RP63, Tamb_C: float = 25.0, capped: bool = True, map_fn=None) -> np.ndarray:
        """Quasi-steady state powers (F5 fixed points if ``capped``) for the analytic formulas [W]."""
        from .power import demand_power, steady_power
        us = self.primitives()
        if map_fn is not None:
            us = [map_fn(n, u) for n, u in zip(self.state_map, us)]
        f = (lambda u: steady_power(u, p, Tamb_C)) if capped else (lambda u: demand_power(u, p))
        return np.array([f(u) for u in us])

    # -------------------------------------------------------------------- path sampling
    def sample_paths(self, n: int, t_max: float, h: float = 10.0, seed: int = 7, x0: int = 0) -> Paths:
        """Sample ``n`` independent paths on [0, t_max] (t_max in hours); jump times in seconds.

        Sojourns are exponential with the diagonal rates; the next state follows the embedded jump chain.
        Vectorised over paths (one draw per path per jump).
        """
        rng = np.random.default_rng(seed)
        L = self.generator()
        rates = -np.diag(L)
        J = (L - np.diag(np.diag(L))) / rates[:, None]
        cum = np.cumsum(J, axis=1)
        cum[:, -1] = 1.0
        t_end = t_max * 3600.0
        s = np.full(n, x0, dtype=np.int64)
        t = np.zeros(n)
        times, states = [t.copy()], [s.copy()]
        alive = np.ones(n, bool)
        while alive.any():
            tau = rng.exponential(1.0, n) / rates[s] * 3600.0
            u = rng.random(n)
            t_new = t + tau
            s_new = (u[:, None] > cum[s]).sum(axis=1)
            t = np.where(alive, t_new, t)
            s = np.where(alive, s_new, s)
            times.append(np.where(alive, t, np.inf))
            states.append(s.copy())
            alive &= t < t_end
        return Paths(np.stack(times, axis=1), np.stack(states, axis=1).astype(np.int64), float(h))


# ------------------------------------------------------------------------------------------- levers (§12.1)
UI_STATES = ("RD", "SO", "SO-L", "NV")


def _scr(fn):
    return lambda n, u: fn(u) if np.all(np.asarray(u.s_scr) > 0) else u


# id -> (label, map(primitive name, UsageState) -> UsageState).  Applied to S12 segments and CTMC states.
LEVERS: dict[str, tuple[str, Callable]] = {
    "skip_game": ("Skip 3D gaming (GM -> standby)",
                  lambda n, u: PRIMITIVES["SB-W"] if n in ("GM", "GM-L") else u),
    "dark_theme": ("Dark theme (APL x0.25 on RD, SO, SO-L, NV)",
                   lambda n, u: replace(u, m_APL=0.25) if n in UI_STATES else u),
    "voice_call": ("Voice call instead of video call",
                   lambda n, u: replace(u, s_scr=0, awake=1, c_f=0, x_c=0.05) if n in ("VC", "VC-W") else u),
    "brightness_m01": ("Brightness -0.1 on all screen-on states",
                       _scr(lambda u: replace(u, b=np.maximum(u.b - 0.1, 0.0)))),
    "refresh_60": ("Cap refresh at 60 Hz", _scr(lambda u: replace(u, f_r=60.0))),
    "vcall_wifi": ("Video call over Wi-Fi instead of LTE", lambda n, u: PRIMITIVES["VC-W"] if n == "VC" else u),
    "half_bg": ("Halve background wake-ups", lambda n, u: replace(u, lam=u.lam / 2) if n.startswith("SB") else u),
    "standby_5g": ("Avoid: standby on 5G NSA without Wi-Fi",
                   lambda n, u: replace(PRIMITIVES["SB-5"], lam=u.lam) if n.startswith("SB") else u),
    "weak_signal": ("Avoid: weak signal (-110 dBm) all day", lambda n, u: replace(u, RSRP=-110.0)),
    "brightness_p02": ("Avoid: brightness +0.2 (outdoors)",
                       _scr(lambda u: replace(u, b=np.minimum(u.b + 0.2, 1.0)))),
    "os_coalesce": ("OS: tail-aware coalescing of background traffic (x5)",
                    lambda n, u: replace(u, lam=u.lam / 5) if n.startswith("SB") else u),
    "os_frame_cap": ("OS: thermal-aware frame cap (game x_c 0.80 -> 0.64)",
                     lambda n, u: replace(u, x_c=0.64) if n in ("GM", "GM-L") else u),
}


__all__ = [
    "UsageState", "PRIMITIVES", "BONN_IDLE", "BONN_CASES", "Schedule", "s12_items", "s12_schedule",
    "s9_schedule", "wifi_cliff", "Scenario", "scenario", "SCENARIO_IDS", "ACTS", "MAP_WIFI", "MAP_LTE", "PROFILES",
    "Paths", "PathCursor", "CTMC", "LEVERS", "UI_STATES",
]
