"""Unit and consistency tests of the building blocks (params, cell, power, usage, model internals)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from phonebatt import cell, power
from phonebatt.model import (Inputs, Segment, _scalar_constants, _scalar_state, classify, event_batch,
                             event_components, fast_state)
from phonebatt.params import RP63, UQ_SPEC, USAGE_UQ_SPEC, K0, sample, with_theta
from phonebatt.usage import CTMC, PRIMITIVES, Schedule, s12_schedule


def _rand_states(n, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.uniform(-0.02, 1.0, n)
    v1 = rng.uniform(0.0, 0.08, n)
    T = rng.uniform(250.0, 330.0, n)
    return np.stack([z, v1, T])


def test_fast_state_matches_separate_functions():
    p = RP63.replace(q_loss=0.1)
    x = _rand_states(500)
    inp = Inputs(np.full(500, 0.7), np.linspace(0, 3, 500), np.full(500, 3.0), np.full(500, 268.15))
    st = fast_state(x, inp, p)
    z, v1, T = x
    np.testing.assert_allclose(st["U"], cell.ocv(z, T, p), rtol=0, atol=1e-13)
    np.testing.assert_allclose(st["R0"], cell.R0(z, T, p), rtol=1e-13)
    np.testing.assert_allclose(st["R1"], cell.R1(T, p), rtol=1e-13)
    np.testing.assert_allclose(st["tau1"], cell.tau1(T, p), rtol=1e-13)
    Ts = power.skin_temperature(T, inp.Tamb, p)
    np.testing.assert_allclose(st["Psoc"], power.throttle(inp.dem, Ts, p), rtol=1e-13)
    I, V, D = cell.current(st["Ue"], st["R0"], st["P"])
    np.testing.assert_allclose(st["I"], I, rtol=1e-13)
    np.testing.assert_allclose(st["Vpk"], cell.peak_voltage(st["Ue"], st["R0"], st["P"] + inp.dPb), rtol=1e-13)


def test_scalar_state_matches_vector_state():
    p = RP63.replace(q_loss=0.2)
    c = _scalar_constants(p)
    x = _rand_states(200, seed=3)
    for k in range(x.shape[1]):
        inp = Inputs(0.6, 0.2 + 0.02 * k, 3.0 * (k % 2), 263.15 + 0.2 * k)
        vs = fast_state(x[:, k], inp, p)
        ss = _scalar_state(x[0, k], x[1, k], x[2, k], inp, p, c)
        for key in ("U", "R0", "R1", "tau1", "Psoc", "P", "I", "V", "D", "Vpk"):
            assert math.isclose(ss[key], float(vs[key]), rel_tol=1e-12, abs_tol=1e-14), key


def test_current_root_is_cancellation_free_and_exact():
    Ue = np.array([4.0, 3.7, 3.5])
    R = np.array([0.06, 0.2, 0.3])
    for P in (1e-9, 0.04, 1.0, 5.0):
        I, V, D = cell.current(Ue, R, P)
        np.testing.assert_allclose(V * I, P, rtol=1e-12)
        np.testing.assert_allclose(V, Ue - R * I, rtol=1e-14)
        assert np.all(V > Ue / 2)


def test_display_soc_maps_and_inverse():
    p = RP63
    s = np.linspace(0, 1, 11)
    for mode in ("D1", "D2"):
        z = cell.soc_to_z(s, p, mode)
        np.testing.assert_allclose(cell.display_soc(z, p, mode), s, atol=1e-10)
    # §2.11 check values for the energy fraction at s = 0.5 / 0.2 / 0.1
    z = cell.soc_to_z(np.array([0.5, 0.2, 0.1]), p, "D1")
    np.testing.assert_allclose(cell.display_soc(z, p, "D2"), [0.477, 0.187, 0.092], atol=6e-4)
    # max |s − ς| = 0.024 at s ≈ 0.6
    sg = np.linspace(0, 1, 2001)
    d = sg - cell.display_soc(cell.soc_to_z(sg, p, "D1"), p, "D2")
    assert abs(np.max(np.abs(d)) - 0.024) < 6e-4 and abs(sg[np.argmax(np.abs(d))] - 0.6) < 0.06
    assert cell.ios_display(0.5) == 0.53 and cell.ios_display(0.95) == 1.0


def test_ocv_integral_closed_form_matches_quadrature():
    from scipy.integrate import quad
    p = RP63.replace(dU=0.01)
    for z in (0.05, 0.3, 0.77, 1.0):
        q = quad(lambda zz: float(cell.ocv_ref(zz, p)), 0.0, z, limit=200)[0]
        assert abs(float(cell.ocv_integral(z, p)) - q) < 1e-11


def test_sample_distributions_and_dirichlet():
    th = sample(None, 1024, seed=2026)
    for k, (dist, lo, mode, hi) in UQ_SPEC.items():
        v = th[k]
        assert v.shape == (1024,)
        assert np.all(v >= lo - 1e-12) and np.all(v <= hi + 1e-12), k
    # log-uniform: median is the geometric mean
    lo, hi = UQ_SPEC["R0_ref"][1], UQ_SPEC["R0_ref"][3]
    assert abs(np.median(th["R0_ref"]) / math.sqrt(lo * hi) - 1) < 0.02
    # triangular: mean (a + m + b)/3
    a, m, b = UQ_SPEC["C_th"][1:]
    assert abs(th["C_th"].mean() - (a + m + b) / 3) < 1.0
    u = sample(["p_typical", "lam_on_mul"], 512, seed=1)
    assert u["p_typical"].shape == (7, 512)
    np.testing.assert_allclose(u["p_typical"].sum(axis=0), 1.0, atol=1e-12)
    base = np.array(USAGE_UQ_SPEC["p_typical"][2])
    np.testing.assert_allclose(u["p_typical"].mean(axis=1), base, atol=0.01)
    # reproducible
    np.testing.assert_array_equal(sample(["Q_n"], 16)["Q_n"], sample(["Q_n"], 16)["Q_n"])


def test_params_broadcast_through_power_and_cell():
    th = sample(["P_SoC_max", "R0_ref", "Ea_R0"], 64)
    pb = with_theta(RP63, th)
    Pns, dem, dPb, br = power.power_parts(PRIMITIVES["GM"], pb)
    assert np.shape(dem) == (64,)
    for k in range(0, 64, 13):
        pk = RP63.replace(**{n: float(v[k]) for n, v in th.items()})
        assert math.isclose(float(dem[k]), float(power.power_parts(PRIMITIVES["GM"], pk)[1]), rel_tol=1e-13)
        assert math.isclose(float(cell.R0(0.3, 260.0, pb)[k]), float(cell.R0(0.3, 260.0, pk)), rel_tol=1e-13)


def test_radio_poisson_fractions_and_weak_signal():
    p = RP63
    # standby: tail and promotion terms of the exact Poisson time average (E14)
    u = PRIMITIVES["SB-L"]
    r = p.radio("lte")
    lam = u.lam
    phi = ((r["P_tail"] - r["P_idle"]) * (1 - math.exp(-lam * r["T_tail"])) + (r["P_act"] - r["P_tail"])
           * (1 - math.exp(-lam * r["T_act"])) + lam * math.exp(-lam * r["T_tail"]) * r["E_pro"])
    assert math.isclose(float(power.radio_power(u, p)), r["P_idle"] + phi, rel_tol=1e-12)
    assert math.isclose(float(power.m_sig(-105.0, p)), 1.52, rel_tol=1e-12)
    assert float(power.m_sig(-300.0, p)) == 3.5
    # §3.4 check values: video radio power over Wi-Fi / LTE / 5G NSA
    for name, val in (("VW", 0.074), ("VL", 0.220), ("V5", 0.363)):
        assert abs(float(power.radio_power(PRIMITIVES[name], p)) - val) < 5e-4


def test_soc_demand_check_values_and_inverse():
    xs = np.array([0.05, 0.10, 0.25, 0.30, 0.50, 0.80])
    np.testing.assert_allclose(power.soc_demand(xs, RP63), [0.056, 0.126, 0.427, 0.562, 1.307, 3.163], atol=5e-4)
    np.testing.assert_allclose(power.x_eff(power.soc_demand(xs, RP63), RP63), xs, atol=1e-12)


def test_display_power_check_values():
    # full-white maximum brightness 0.12 + 0.25*4/0.7 = 1.549 W (§3.2)
    u = PRIMITIVES["RD"].with_(b=1.0, APL=1.0, f_r=60.0)
    assert abs(float(power.display_power(u, RP63)) - (0.12 + 0.25 * 4 / 0.7)) < 1e-12
    # dark-mode savings in the S4 state (§3.2 [PROTO]): 5.1 / 8.6 / 12.6 / 33.2 % at b = 0.3 / 0.4 / 0.5 / 1.0
    for b, sav in ((0.3, 5.1), (0.4, 8.6), (0.5, 12.6), (1.0, 33.2)):
        base = PRIMITIVES["SO"].with_(b=b)
        dark = base.with_(APL=0.15)
        P0, P1 = power.demand_power(base, RP63), power.demand_power(dark, RP63)
        assert abs(100 * (1 - P1 / P0) - sav) < 0.06


def test_schedule_lookup_and_s12():
    s = s12_schedule()
    assert s.period == 86400.0
    assert s.index_at(0.0) == 0 and s.index_at(30 * 60 - 1e-6) == 0 and s.index_at(30 * 60) == 1
    assert s.index_at(86400.0 + 1.0) == 0
    assert all(d % 10 == 0 for d in s.durations)


def test_ctmc_generator_and_paths():
    c = CTMC("typical")
    L = c.generator()
    assert L.shape == (16, 16)
    np.testing.assert_allclose(L.sum(axis=1), 0.0, atol=1e-12)
    paths = c.sample_paths(200, 48.0, h=10.0, seed=3)
    tm = paths.times
    fin = np.isfinite(tm[:, 1:])
    assert np.all((tm[:, 1:][fin] - tm[:, :-1][fin]) > 0)
    cur = paths.cursor()
    for k in (0, 5, 100, 1000, 17279):
        np.testing.assert_array_equal(cur.state_for_step(k), paths.state_for_step(k))
    # empirical occupation matches π
    from phonebatt.analytic import stationary
    pi = stationary(L)
    tt = np.arange(0, 48 * 3600, 60.0)
    occ = np.zeros(16)
    for t in tt[::10]:
        occ += np.bincount(paths.state_at(t), minlength=16)
    occ /= occ.sum()
    assert np.max(np.abs(occ - pi)) < 0.03


def test_classify_priorities():
    # rows: charge, cutoff, brownout, thermal (normalised), D/Ũ²; one column per sample
    comps = np.array([[5.0, -1.0, -3.0, 9.0, 1.0],   # cutoff and brownout violated -> cutoff (priority)
                      [5.0, 0.5, -0.1, 9.0, 1.0],    # brownout only
                      [5.0, 0.01, 0.002, 9.0, 1.0],  # interpolated crossing, none <= 0 -> argmin
                      [-0.1, 2.0, 2.0, 9.0, 1.0],    # charge
                      [5.0, -1.0, -1.0, -1.0, -1e-3]]).T  # D <= 0 -> infeasible
    assert list(classify(comps)) == ["cutoff", "brownout", "brownout", "charge", "infeasible"]


def test_event_batch_nan_when_infeasible():
    x = np.array([[0.5], [0.0], [253.15]])
    inp = Inputs(np.array([25.0]), np.array([0.0]), np.array([0.0]), np.array([253.15]))
    assert np.isnan(event_batch(x, inp, RP63))[0]
    seg = Segment(PRIMITIVES["SO"].with_(P_extra=25.0), 253.15, RP63)
    assert seg.first_violated(x[:, 0]) == "infeasible"
