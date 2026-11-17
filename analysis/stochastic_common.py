"""Shared helpers of the 'stochastic' analysis group (AN5 ``run_ctmc.py``, AN11 ``run_usage_fluct.py``).

FINAL_SPEC §5.3 (usage CTMC), §7 F6–F10 (closed forms), §10.5 (usage fluctuations).

* ``scaled_ctmc(profile, kappa)``: the S13/S14 CTMC with every sojourn multiplied by κ (π fixed, §10.5 item 1).
* ``UsageSim``: a vectorised event-driven simulator of the activity × connectivity process with general session
  laws (exponential or lognormal with a given CV; idle and connectivity sojourns stay exponential, F9/X15) and
  optional per-visit jitter (one jitter draw per activity sojourn, kept across connectivity switches, §4.6).
  It returns exact L0 time-to-empty samples (energy E_u depleted at the state power, no time step) or stores the
  jump paths as ``phonebatt.usage.Paths`` so that ``model.run_batch`` integrates the full L2 model on them.
* ``ExpandedCTMC``: a CTMC whose ``primitives()`` is an arbitrary list of inputs (e.g. 16 states × K jitter
  variants); ``run_batch`` only needs ``primitives()`` when paths are supplied.
* ``jitter_bank`` / ``jittered_inputs``: the ±30 % per-visit multipliers on b, APL, λ, d_dl (U(0.7, 1.3) each,
  scrambled Sobol', fixed seed).
* ``l0_from_paths``: exact L0 TTE on stored paths (pairs with the L2 result on the same path).
* ``summ``: moments / quantiles / Pr(TTE < 24 h) of a TTE sample with standard errors.

Nothing here changes ``src/phonebatt`` (frozen core); the CTMC conventions follow ``phonebatt.usage.CTMC``:
product state index = 2·activity + connectivity (activity order ``ACTS``, Wi-Fi first), rates in 1/h, path times
in seconds.
"""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np

from scipy.special import ndtri

from phonebatt import RP63, analytic, cell
from phonebatt.params import unit_samples
from phonebatt.power import steady_power
from phonebatt.usage import ACTS, CTMC, PROFILES, Paths, UsageState

EU_WH = float(cell.usable_energy(RP63))
SESSIONS = ACTS[1:]


# ------------------------------------------------------------------------------------------- CTMC variants
def scaled_ctmc(profile: str, kappa: float = 1.0, **kw) -> CTMC:
    """CTMC of ``profile`` with every sojourn (sessions, idle, connectivity) multiplied by ``kappa``.

    Keyword overrides (lam_on, p, d_min, mean_w, mean_c) are applied before scaling.
    """
    lam0, p0, d0 = PROFILES[profile]
    lam = kw.pop("lam_on", lam0)
    d = kw.pop("d_min", d0)
    mw = kw.pop("mean_w", 2.5)
    mc = kw.pop("mean_c", 2.5)
    return CTMC(profile, mean_w=mw * kappa, mean_c=mc * kappa, lam_on=lam / kappa,
                d_min={a: v * kappa for a, v in d.items()}, **kw)


class ExpandedCTMC(CTMC):
    """CTMC whose inputs are an explicit list (paths must be supplied to ``run_batch``)."""

    def __init__(self, base: CTMC, inputs: list[UsageState]):
        self.__dict__.update(base.__dict__)
        self._inputs = list(inputs)

    def primitives(self) -> list[UsageState]:
        return self._inputs


def state_powers(ctmc: CTMC, p=RP63, Tamb_C: float = 25.0) -> np.ndarray:
    """Quasi-steady (F5) state powers of the 16 product states [W]."""
    return ctmc.state_powers(p, Tamb_C, capped=True)


# ------------------------------------------------------------------------------------------- per-visit jitter
JITTER_FIELDS = ("b", "APL", "lam", "d_dl")


def jitter_bank(K: int = 256, lo: float = 0.7, hi: float = 1.3, seed: int = 4242) -> np.ndarray:
    """(K, 4) multipliers U(lo, hi) on (b, APL, λ, d_dl), scrambled Sobol' (§4.6 per-visit jitter)."""
    return lo + (hi - lo) * unit_samples(len(JITTER_FIELDS), K, method="sobol", seed=seed)


def jittered_inputs(u: UsageState, mult: np.ndarray) -> UsageState:
    """Array-valued input: ``u`` with b, APL, λ, d_dl multiplied by the columns of ``mult`` (K, 4)."""
    return replace(u, b=u.b * mult[:, 0], APL=u.APL * mult[:, 1], lam=u.lam * mult[:, 2], d_dl=u.d_dl * mult[:, 3])


def jitter_power_bank(ctmc: CTMC, mult: np.ndarray, p=RP63, Tamb_C: float = 25.0) -> np.ndarray:
    """(16, K) quasi-steady powers of every product state under each jitter vector [W]."""
    return np.stack([np.broadcast_to(steady_power(jittered_inputs(u, mult), p, Tamb_C), (mult.shape[0],))
                     for u in ctmc.primitives()])


def jitter_expanded_inputs(ctmc: CTMC, mult: np.ndarray) -> list[UsageState]:
    """16·K scalar inputs, index = state·K + k (for ``ExpandedCTMC`` and ``UsageSim(jitter_K=K)``)."""
    out = []
    for u in ctmc.primitives():
        for m in mult:
            out.append(replace(u, name=f"{u.name}*", b=u.b * m[0], APL=u.APL * m[1], lam=u.lam * m[2],
                               d_dl=u.d_dl * m[3]))
    return out


# ------------------------------------------------------------------------------------------- simulator
class UsageSim:
    """Event-driven activity × connectivity simulator with general session laws (vectorised over paths).

    Parameters
    ----------
    ctmc : the CTMC whose rates define the process (``lam_on``, ``p``, ``d_min``, ``mean_w``, ``mean_c``).
    law : 'exp' (Markov; reproduces the CTMC) or 'lognormal' (sessions only; idle stays exponential).
    cv : coefficient of variation of the lognormal session durations.
    jitter_K : number of jitter variants (0 = none); a variant index k ~ U{0..K−1} is drawn at every
        activity jump and kept across connectivity jumps.  Expanded state index = (2a + c)·K + k.
    connectivity : False freezes the connectivity on Wi-Fi (hub-and-spoke chain alone, F8/F9 exact).
    """

    def __init__(self, ctmc: CTMC, law: str = "exp", cv: float = 1.0, jitter_K: int = 0,
                 connectivity: bool = True):
        if law not in ("exp", "lognormal"):
            raise ValueError(law)
        self.lam_on = float(ctmc.lam_on)  # 1/h
        self.p = np.array([ctmc.p[a] for a in SESSIONS], dtype=float)
        self.p = self.p / self.p.sum()
        self.cum_p = np.cumsum(self.p)
        self.cum_p[-1] = 1.0
        self.d = np.array([ctmc.d_min[a] / 60.0 for a in SESSIONS], dtype=float)  # h
        self.mean_c = np.array([ctmc.mean_w, ctmc.mean_c], dtype=float)  # h, sojourn on Wi-Fi / cellular
        self.law, self.cv, self.K, self.conn = law, float(cv), int(jitter_K), bool(connectivity)
        s2 = math.log(1.0 + self.cv ** 2)
        self._ln_sig = math.sqrt(s2)
        self._ln_mu = np.log(self.d) - 0.5 * s2

    # -- draws.  Common random numbers: a session duration is a monotone transform of one uniform for every
    # law, so runs with the same seed but different laws are coupled path by path; the jitter indices use a
    # separate stream, so the activity/connectivity paths do not depend on jitter_K.
    def _session(self, rng, j):
        u = rng.random(j.size)
        if self.law == "exp":
            return -np.log1p(-u) * self.d[j]
        return np.exp(self._ln_mu[j] + self._ln_sig * ndtri(u))

    def _idle(self, rng, m):
        return rng.exponential(1.0 / self.lam_on, m)

    def _conn(self, rng, c):
        if not self.conn:
            return np.full(c.size, np.inf)
        return rng.exponential(1.0, c.size) * self.mean_c[c]

    def _kidx(self, rng, m):
        return self._rng_j.integers(0, self.K, m) if self.K else np.zeros(m, dtype=np.int64)

    def _index(self, a, c, k):
        return (2 * a + c) * max(self.K, 1) + k

    def run(self, n: int, seed: int, P: np.ndarray | None = None, Eu: float = EU_WH, store: bool = False,
            t_max_h: float = 400.0, x0: int = 0):
        """Simulate ``n`` paths from product state ``x0`` (activity idle and Wi-Fi by default).

        With ``store=False`` returns ``(tte_h, state_at_depletion)``: exact L0 TTE with energy ``Eu`` [Wh]
        depleted at ``P[expanded state]`` [W]; paths without depletion before ``t_max_h`` get NaN.
        With ``store=True`` returns ``Paths`` (jump times in s up to ``t_max_h``) for ``run_batch``.
        """
        rng = np.random.default_rng(seed)
        self._rng_j = np.random.default_rng([seed, 7919])
        a = np.full(n, x0 // 2, dtype=np.int64)
        c = np.full(n, x0 % 2, dtype=np.int64)
        if not self.conn:
            c[:] = 0
        k = self._kidx(rng, n)
        ra = np.where(a == 0, self._idle(rng, n), 0.0)
        if (a > 0).any():
            ra[a > 0] = self._session(rng, a[a > 0] - 1)
        rc = self._conn(rng, c)
        t = np.zeros(n)
        idx = np.arange(n)  # live paths (compacted)
        if store:
            times, states = [np.zeros(n)], [self._index(a, c, k)]
        else:
            P = np.asarray(P, dtype=float)
            E = np.full(n, float(Eu))
            tte = np.full(n, np.nan)
            xdep = np.full(n, -1, dtype=np.int64)
        while idx.size:
            dt = np.minimum(ra, rc)
            if not store:
                s = self._index(a, c, k)
                e = P[s] * dt
                done = e >= E
                if done.any():
                    tte[idx[done]] = t[done] + E[done] / P[s[done]]
                    xdep[idx[done]] = s[done]
                E = E - e
            t = t + dt
            act = ra <= rc
            # activity jumps
            ia = np.nonzero(act)[0]
            if ia.size:
                was_idle = a[ia] == 0
                st = ia[was_idle]
                if st.size:
                    j = np.searchsorted(self.cum_p, rng.random(st.size), side="right")
                    j = np.minimum(j, self.p.size - 1)
                    a[st] = j + 1
                    ra[st] = self._session(rng, j)
                en = ia[~was_idle]
                if en.size:
                    a[en] = 0
                    ra[en] = self._idle(rng, en.size)
                k[ia] = self._kidx(rng, ia.size)
                rc[ia] -= dt[ia]
            # connectivity jumps
            ic = np.nonzero(~act)[0]
            if ic.size:
                c[ic] = 1 - c[ic]
                rc[ic] = self._conn(rng, c[ic])
                ra[ic] -= dt[ic]
            if store:
                tt = np.full(n, np.inf)
                tt[idx] = t
                ss = np.zeros(n, dtype=np.int64)
                ss[idx] = self._index(a, c, k)
                times.append(tt * 3600.0)
                states.append(ss)
                keep = t < t_max_h
            else:
                keep = (E > 0) & (t < t_max_h)
            if not keep.all():
                idx, a, c, k, ra, rc, t = (v[keep] for v in (idx, a, c, k, ra, rc, t))
                if not store:
                    E = E[keep]
        if store:
            T = np.stack(times, axis=1)
            S = np.stack(states, axis=1)
            # Paths convention: after a path's horizon the times are +inf (the state there is irrelevant)
            return Paths(T, S, 10.0)
        return tte, xdep


# ------------------------------------------------------------------------------------------- path utilities
def l0_from_paths(paths: Paths, P: np.ndarray, Eu: float = EU_WH) -> np.ndarray:
    """Exact L0 TTE [h] on stored paths (``P`` indexed by the path state); NaN if the horizon is too short."""
    t = paths.times / 3600.0
    st = paths.states
    with np.errstate(invalid="ignore"):  # inf - inf in the padding
        dt = np.diff(t, axis=1)
    dt = np.where(np.isfinite(dt), dt, 0.0)
    e = np.asarray(P, float)[st[:, :-1]] * dt
    cum = np.cumsum(e, axis=1)
    hit = cum >= Eu
    any_hit = hit.any(axis=1)
    i = np.argmax(hit, axis=1)
    rows = np.arange(st.shape[0])
    before = np.where(i > 0, cum[rows, np.maximum(i - 1, 0)], 0.0)
    out = t[rows, i] + (Eu - before) / np.asarray(P, float)[st[rows, i]]
    return np.where(any_hit, out, np.nan)


def subset_paths(paths: Paths, rows) -> Paths:
    rows = np.asarray(rows)
    return Paths(paths.times[rows].copy(), paths.states[rows].copy(), paths.h)


# ------------------------------------------------------------------------------------------- statistics
def summ(x, t_star: float = 24.0) -> dict:
    """Moments, quantiles and Pr(TTE < t*) of a TTE sample [h], with standard errors (NaN/inf dropped).

    Also Pr(TTE < 16 h) (§9 item 5: a day started at 07:00 that dies before 23:00).
    """
    x = np.asarray(x, dtype=float)
    n_all = x.size
    x = x[np.isfinite(x)]
    n = x.size
    m, sd = float(x.mean()), float(x.std(ddof=1))
    p = float(np.mean(x < t_star))
    sk = float(np.mean(((x - m) / sd) ** 3)) if sd > 0 else math.nan
    return dict(n=n, n_censored=int(n_all - n), mean_h=m, sd_h=sd, se_mean_h=sd / math.sqrt(n),
                se_sd_h=sd * math.sqrt(max(kurt(x) - 1.0, 0.0) / (4.0 * n)) if n > 3 else math.nan,
                cv=sd / m, skew=sk, q05_q50_q95_h=[float(v) for v in np.percentile(x, [5, 50, 95])],
                p_lt_24h=p, se_p_lt_24h=math.sqrt(max(p * (1 - p), 1.0 / n) / n),
                p_lt_16h=float(np.mean(x < 16.0)))


def kurt(x) -> float:
    """Kurtosis (non-excess) of a sample; used for the SE of the SD: Var(s) ≈ σ²(κ − 1)/(4n)."""
    x = np.asarray(x, float)
    s = x.std()
    return float(np.mean(((x - x.mean()) / s) ** 4)) if s > 0 else 3.0


def ig_summary(mu: float, lam: float) -> dict:
    """IG(mean μ, shape λ) quantiles, Pr(T < 24 h), SD and skewness (F7)."""
    q = analytic.ig_quantiles(mu, lam)
    return dict(mean_h=mu, sd_h=math.sqrt(mu ** 3 / lam), q05_q50_q95_h=[float(v) for v in q],
                p_lt_24h=float(analytic.ig_cdf(24.0, mu, lam)), skew=3.0 * math.sqrt(mu / lam))


def boot_ci(x, fn, n_boot: int = 1000, seed: int = 0, alpha: float = 0.05):
    """Percentile bootstrap CI of a statistic of a 1-D sample."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    rng = np.random.default_rng(seed)
    vals = np.array([fn(x[rng.integers(0, x.size, x.size)]) for _ in range(n_boot)])
    return [float(np.percentile(vals, 100 * alpha / 2)), float(np.percentile(vals, 100 * (1 - alpha / 2)))]


# ------------------------------------------------------------------------------------------- closed forms
def hub_inputs(ctmc: CTMC, P16: np.ndarray, average_connectivity: bool = False):
    """(lam_on, p_j, d_j [h], P_idle, P_j) of the hub-and-spoke chain: Wi-Fi powers, or π_c-averaged powers."""
    Pa = P16.reshape(-1, 2)
    if average_connectivity:
        w = np.array([1.0 / ctmc.mean_w, 1.0 / ctmc.mean_c])[::-1]
        w = w / w.sum()  # π_wifi = mean_w/(mean_w + mean_c)
        Pa = Pa @ w
    else:
        Pa = Pa[:, 0]
    p = np.array([ctmc.p[a] for a in SESSIONS])
    d = np.array([ctmc.d_min[a] / 60.0 for a in SESSIONS])
    return float(ctmc.lam_on), p, d, float(Pa[0]), Pa[1:]


def f9_full_sigma2(ctmc: CTMC, P16: np.ndarray, cv2: float) -> dict:
    """Approximate CLT variance rate of the full 16-state chain with general session laws (F9 increment).

    σ²_full(CV) ≈ σ²_F7,full + [σ²_F9,hub(CV) − σ²_F8,hub] with the hub at connectivity-averaged powers:
    the activity part changes exactly as F9 says; the connectivity part is kept at its Markov value.
    """
    m = analytic.ctmc_moments(ctmc.generator(), P16, EU_WH)
    lam, p, d, P0, Pj = hub_inputs(ctmc, P16, average_connectivity=True)
    s_exp = analytic.semi_markov_sigma2(lam, p, d, P0, Pj, 1.0, 1.0)
    s_cv = analytic.semi_markov_sigma2(lam, p, d, P0, Pj, cv2, 1.0)
    s2 = m["sigma2"] + (s_cv - s_exp)
    return dict(sigma2=s2, sigma2_exp=m["sigma2"], hub_exp=s_exp, hub_cv=s_cv, Pbar=m["Pbar"],
                sd_clt_h=math.sqrt(s2 * EU_WH / m["Pbar"] ** 3), mu_dynkin_exp_h=m["mu_dynkin"])


def json_list(x, nd: int = 6):
    return [None if not np.isfinite(v) else round(float(v), nd) for v in np.asarray(x, float).ravel()]
