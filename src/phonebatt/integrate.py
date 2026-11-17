"""Vectorized fixed-step integrators with terminal-event detection.

Used for batch runs (Monte Carlo, Sobol) where thousands of parameter samples
are integrated simultaneously.  Each sample is an independent ODE system

    dx/dt = f(t, x; theta),    x in R^n,

stacked as an array of shape (n_state, n_samples).  Integration of a sample
stops when its event function g(t, x) crosses zero from above; the crossing time
is located by linear interpolation inside the step, which is second-order
accurate and adequate for TTE (h ~ seconds, TTE ~ hours).
"""
from __future__ import annotations

from typing import Callable

import numpy as np

RHS = Callable[[float, np.ndarray], np.ndarray]
EVENT = Callable[[float, np.ndarray], np.ndarray]


def rk4_batch(
    f: RHS,
    x0: np.ndarray,
    t_end: float,
    h: float,
    event: EVENT | None = None,
    record_every: int = 0,
):
    """Integrate a batch of ODE systems with classical RK4.

    Parameters
    ----------
    f : callable (t, x) -> dx/dt, x has shape (n_state, n_samples).
    x0 : initial state, shape (n_state, n_samples).
    t_end : final time [s].
    h : step size [s].
    event : callable (t, x) -> g, shape (n_samples,).  A sample terminates
        when g changes sign from positive to non-positive.  NaN in g is treated
        as an immediate event (e.g. infeasible power demand).
    record_every : if > 0, store the state every ``record_every`` steps.

    Returns
    -------
    t_event : array (n_samples,), event time (np.inf if never reached).
    x_final : state at event (or t_end).
    history : (times, states) if record_every > 0 else None.
    """
    x = np.array(x0, dtype=float, copy=True)
    n = x.shape[1]
    alive = np.ones(n, dtype=bool)
    t_event = np.full(n, np.inf)
    x_final = x.copy()
    times, states = [], []
    g_prev = event(0.0, x) if event is not None else None
    if g_prev is not None:
        dead0 = ~(g_prev > 0)
        t_event[dead0] = 0.0
        alive &= ~dead0
    t = 0.0
    k = 0
    n_steps = int(np.ceil(t_end / h))
    for k in range(n_steps):
        if not alive.any():
            break
        if record_every and k % record_every == 0:
            times.append(t)
            states.append(x.copy())
        k1 = f(t, x)
        k2 = f(t + 0.5 * h, x + 0.5 * h * k1)
        k3 = f(t + 0.5 * h, x + 0.5 * h * k2)
        k4 = f(t + h, x + h * k3)
        x_new = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        x_new[:, ~alive] = x[:, ~alive]
        t_new = t + h
        if event is not None:
            g_new = event(t_new, x_new)
            bad = ~np.isfinite(g_new)
            crossed = alive & ((g_new <= 0) | bad)
            if crossed.any():
                gp = g_prev[crossed]
                gn = np.where(bad[crossed], 0.0, g_new[crossed])
                denom = gp - gn
                frac = np.where(denom > 0, gp / np.where(denom > 0, denom, 1.0), 1.0)
                frac = np.clip(frac, 0.0, 1.0)
                t_event[crossed] = t + frac * h
                x_final[:, crossed] = x[:, crossed] + frac * (x_new[:, crossed] - x[:, crossed])
                alive &= ~crossed
            g_prev = np.where(alive, g_new, g_prev)
        x = x_new
        t = t_new
    x_final[:, alive] = x[:, alive]
    hist = (np.array(times), np.array(states)) if record_every else None
    return t_event, x_final, hist



def solve_ivp_segmented(segments, x0, t_max: float, problem, method: str = "LSODA", **solver_kw):
    """Thin segment-by-segment wrapper around ``scipy.integrate.solve_ivp`` (FINAL_SPEC §6.2).

    The input is piecewise constant: ``segments`` is a sequence of ``(duration_s, payload)`` repeated
    cyclically.  At each switch ``problem(j, payload, t, x)`` is called with the cyclic segment index j and the
    carried state; it returns ``(rhs, events)`` for solve_ivp, or ``None`` to stop before integrating (e.g. an
    event is already violated at the switching instant).  The solver is restarted at every switch.

    Yields ``(k, j, t_start, x_start, sol)`` per segment (``sol`` is None when ``problem`` stopped).  The caller
    stops iterating once it has found its terminal event; iteration also ends at ``t_max``.
    """
    from scipy.integrate import solve_ivp

    x = np.array(x0, dtype=float, copy=True)
    t, k, m = 0.0, 0, len(segments)
    while t < t_max:
        j = k % m
        dur, payload = segments[j]
        prob = problem(j, payload, t, x)
        if prob is None:
            yield k, j, t, x.copy(), None
            return
        rhs, events = prob
        sol = solve_ivp(rhs, (t, min(t + dur, t_max)), x, method=method, events=events, **solver_kw)
        if not sol.success:
            raise RuntimeError(f"solve_ivp failed in segment {k}: {sol.message}")
        yield k, j, t, x.copy(), sol
        x, t = sol.y[:, -1].copy(), float(sol.t[-1])
        k += 1
