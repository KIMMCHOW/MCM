"""Numerical acceptance tests T1–T5 (FINAL_SPEC §6.4)."""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from phonebatt import PRIMITIVES, RP63, CTMC, aged, cell, run, run_batch
from phonebatt.analytic import tte_L1
from phonebatt.model import Inputs, Segment, fast_state, rhs_from_state
from phonebatt.params import K0
from phonebatt.usage import Schedule, UsageState, s9_schedule


# ------------------------------------------------------------------------------------------- T1
@pytest.mark.slow
def test_T1_rk4_convergence_against_lsoda():
    """|TTE_RK4(h=2 s) − TTE_LSODA|/TTE < 0.2 % and h = 2 s vs 1 s < 0.1 % on S4, S7, S9, S10a."""
    so = Schedule.constant(PRIMITIVES["SO"])
    cases = []
    # S4 (25 °C fresh) and S10a (−10 °C, SOH 0.8) share the SO schedule: one batch of two samples
    pb_theta = dict(q_loss=np.array([0.0, 0.2]))
    for h in (2.0, 1.0):
        b = run_batch(so, np.array([25.0, -10.0]), pb_theta, h=h)
        cases.append(("S4/S10a", h, b["tte_h"]))
    ref = np.array([run(so, 25.0).tte_h, run(so, -10.0, p=aged(RP63, 0.2)).tte_h])
    b2, b1 = cases[0][2], cases[1][2]
    assert np.all(np.abs(b2 - ref) / ref < 2e-3), (b2, ref)
    assert np.all(np.abs(b2 - b1) / b1 < 1e-3), (b2, b1)
    # S7 gaming at 25 °C
    gm = Schedule.constant(PRIMITIVES["GM"])
    ref = run(gm, 25.0).tte_h
    t2 = run_batch(gm, 25.0, None, h=2.0)["tte_h"][0]
    t1 = run_batch(gm, 25.0, None, h=1.0)["tte_h"][0]
    assert abs(t2 - ref) / ref < 2e-3 and abs(t2 - t1) / t1 < 1e-3, (t2, t1, ref)
    # S9 cold hike at −10 °C, SOH 0.8 (T0 = 20 °C)
    s9 = s9_schedule(-10.0)
    p8 = aged(RP63, 0.2)
    ref = run(s9, -10.0, T0_C=20.0, p=p8).tte_h
    t2 = run_batch(s9, -10.0, None, h=2.0, T0_C=20.0, p=p8)["tte_h"][0]
    t1 = run_batch(s9, -10.0, None, h=1.0, T0_C=20.0, p=p8)["tte_h"][0]
    assert abs(t2 - ref) / ref < 2e-3 and abs(t2 - t1) / t1 < 1e-3, (t2, t1, ref)


@pytest.mark.slow
def test_T1b_ctmc_step_convergence():
    """CTMC: mean TTE on 100 paths at h = 10 s vs 5 s (same continuous-time paths) < 0.5 %."""
    c = CTMC("typical")
    paths = c.sample_paths(100, 120.0, h=10.0, seed=5)
    m10 = np.nanmean(run_batch(c, 25.0, None, h=10.0, paths=paths, t_max_h=120.0)["tte_h"])
    m5 = np.nanmean(run_batch(c, 25.0, None, h=5.0, paths=paths, t_max_h=120.0)["tte_h"])
    assert abs(m10 - m5) / m5 < 5e-3, (m10, m5)


# ------------------------------------------------------------------------------------------- T2
def test_T2_electrical_energy_closure():
    """∫I·U dt = ∫P dt + ∫I²R0 dt + ∫v1²/R1 dt + ½ C1 v1(t_end)², C1 = τ1/R1, isothermal, v1(0) = 0.

    ∫I·U dt is computed independently as 3600 Q_max [F(z0) − F(z_end)] with the closed-form OCV integral F
    (T = T_ref, so the entropic term vanishes); the loss integrals are quadrature states of the ODE.
    """
    p = RP63.replace(C_th=1e12, P_susp=0.0)
    T0 = 25.0 + K0
    # a two-level load (constant power steps) so that v1 relaxes repeatedly
    loads = [(900.0, UsageState("hi", P_extra=2.5)), (900.0, UsageState("lo", P_extra=0.4))]
    ze = float(cell.z_empty(p))
    x = np.array([1.0 - 0.0 * ze, 0.0, T0, 0.0, 0.0, 0.0])
    z0 = x[0]
    t = 0.0
    for k in range(9):  # 2.25 h, ends mid-"hi" segment's successor so that v1(t_end) > 0
        dur, u = loads[k % 2]
        seg = Segment(u, T0, p)

        def f(tt, y, seg=seg):
            st = fast_state(y[:3], seg.inp, p)
            dz, dv1, dT = rhs_from_state(y[:3], st, seg.inp, p)
            return [dz, dv1, dT, st["P"], st["I"] ** 2 * st["R0"], y[1] ** 2 / st["R1"]]

        sol = solve_ivp(f, (t, t + dur), x, method="DOP853", rtol=1e-11, atol=1e-12)
        x, t = sol.y[:, -1], sol.t[-1]
    z_end, v1_end = x[0], x[1]
    R1 = float(cell.R1(T0, p))
    C1 = float(cell.tau1(T0, p)) / R1
    lhs = 3600.0 * float(p.Q_max) * float(cell.ocv_integral(z0, p) - cell.ocv_integral(z_end, p))
    rhs = x[3] + x[4] + x[5] + 0.5 * C1 * v1_end ** 2
    assert v1_end > 1e-3
    assert abs(lhs - rhs) / lhs < 1e-4, (lhs, rhs, (lhs - rhs) / lhs)


# ------------------------------------------------------------------------------------------- T3
@pytest.mark.parametrize("P", [0.5, 1.06, 1.49, 2.52, 3.5])
def test_T3_isothermal_L2_vs_F2(P):
    """Isothermal (C_th = 1e12), no bursts, constant P: L2 (LSODA) vs F2 quadrature < 0.05 % (PROTO ≤ 0.001 %)."""
    p = RP63.replace(C_th=1e12, P_susp=0.0, burst=False)
    r = run(UsageState("const", P_extra=P), 25.0, p=p)
    f2 = float(tte_L1(P, p, 25.0))
    assert r.event == "cutoff"
    assert abs(r.tte_h - f2) / f2 < 5e-4, (r.tte_h, f2)
    assert abs(r.tte_h - f2) / f2 < 2e-5  # the [PROTO] level (≤ 0.001 %, allowing 0.002 %)


# ------------------------------------------------------------------------------------------- T4
def test_T4_sanity_no_nan_monotone_z_thermal_steady_state():
    # (a) no NaN on feasible paths and z monotone (RK4 batch with recorded history, several parameter draws)
    from phonebatt.params import sample
    th = sample(["R0_ref", "Ea_R0", "C_th", "R_th", "P_SoC_max", "m_b"], 16, seed=11)
    b = run_batch(Schedule.constant(PRIMITIVES["SO"]), 0.0, th, h=10.0, t_max_h=8.0, record=True)
    X = b["traj"]["x"]  # (steps, 3, n)
    assert np.all(np.isfinite(X))
    assert np.all(np.diff(X[:, 0, :], axis=0) <= 0.0)
    # LSODA trajectory of the scheduled day: monotone z, finite values
    r = run(Schedule.from_minutes([(30, "RD"), (30, "GM"), (60, "SB-W")]), 25.0, record=True)
    for k in ("z", "V", "I", "T", "P"):
        assert np.all(np.isfinite(r.traj[k])), k
    assert np.all(np.diff(r.traj["z"]) <= 0.0)
    # (b) thermal steady state T -> T_amb + k_env R_th q (1e-3 K): constant load held long after τ_th.  A very
    # large capacity freezes the slow drift of q with U(z), so the state is a true steady state (with the
    # RP-63 capacity the quasi-steady lag τ_th·dT_ss/dt is ≈ 1.3 mK at 1 W).
    p = RP63.replace(k_env=1.2, Q_n=1000.0)
    u = UsageState("const", P_extra=1.0)
    rr = run(u, 10.0, p=p, record=True, t_max_h=6.0)
    tr = rr.traj
    j = np.argmin(np.abs(tr["t_h"] - 5.0))
    x = np.array([tr["z"][j], tr["v1"][j], tr["T"][j]])
    seg = Segment(u, 10.0 + K0, p)
    st = seg.state(x)
    q = st["I"] * st["V"] + st["I"] ** 2 * st["R0"] + x[1] ** 2 / st["R1"] - st["I"] * x[2] * st["dUdT"]
    T_ss = 10.0 + K0 + float(p.k_env * p.R_th) * q
    assert abs(x[2] - T_ss) < 1e-3, x[2] - T_ss


# ------------------------------------------------------------------------------------------- T5
def test_T5_forced_20W_step_is_infeasible():
    """A forced 20 W demand step on a cold pack (P_max < 20 W) ends with event 'infeasible' at the switch."""
    sch = Schedule([(600.0, PRIMITIVES["SO"]), (600.0, UsageState("step20W", s_scr=1, awake=1, P_extra=20.0))])
    p = RP63
    # check the premise: P_max = Ũ²/(4 R0) < 20 W at the switch
    r = run(sch, -20.0, s0=0.5, p=p)
    assert r.event == "infeasible" and abs(r.tte_h * 3600 - 600.0) < 1e-6
    b = run_batch(sch, -20.0, None, s0=0.5, h=2.0, p=p)
    assert b["event"][0] == "infeasible" and abs(b["tte_h"][0] * 3600 - 600.0) < 1e-6
    # NaN path inside rk4_batch: the combined event function is NaN at the switch
    from phonebatt.model import event_batch
    from phonebatt.power import power_parts
    Pns, dem, dPb, _ = power_parts(sch[1][1], p)
    x = np.array([[b["z_end"][0]], [0.0], [b["T_end_C"][0] + K0]])
    g = event_batch(x, Inputs(np.array([Pns]), np.array([dem]), np.array([dPb]), np.array([253.15])), p)
    assert np.isnan(g[0])
