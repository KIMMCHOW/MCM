"""Regression tests against the [PROTO] numbers of FINAL_SPEC (table §15.3 and the tables it points to).

[PROTO] numbers are regression targets computed by docs/design/scripts/final_proto.py; they are not validation
results.  Tolerances are those of §15.3 (relative unless stated).
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from phonebatt import (CTMC, PRIMITIVES, RP63, aged, analytic, cell, power, run, run_batch, s12_schedule,
                       s9_schedule)
from phonebatt import aging as G
from phonebatt import validation as V
from phonebatt.params import K0
from phonebatt.usage import ACTS, MAP_WIFI, PROFILES, Schedule

from conftest import lsoda


def rel(a, b):
    return abs(a - b) / abs(b)


# ------------------------------------------------------------------------------------------- cell
def test_ocv_and_empty_point():
    p = RP63
    assert abs(float(cell.ocv(1.0, 298.15, p)) - 4.4496) < 1e-4
    assert abs(float(cell.mean_ocv(p)) - 3.8699) < 1e-4
    assert abs(float(cell.z_empty(p)) - 0.0263) < 1e-4
    for z, U in ((0.1, 3.644), (0.2, 3.694), (0.3, 3.728), (0.5, 3.815), (0.8, 4.088), (1.0, 4.450)):
        assert abs(float(cell.ocv(z, 298.15, p)) - U) < 6e-4


def test_energies():
    assert rel(float(cell.nominal_energy(RP63)), 18.189) < 1e-3
    assert rel(float(cell.usable_energy(RP63)), 17.788) < 1e-3
    assert rel(float(cell.usable_energy(V.iphone14pro())), 12.565) < 1e-3


@pytest.mark.parametrize("Tc,R0m,R1m,tau", [(40, 42.4, 5.6, 20.0), (25, 60.0, 10.0, 30.0), (10, 88.1, 19.0, 47.0),
                                            (0, 116.6, 30.3, 65.1), (-10, 157.6, 50.0, 92.6),
                                            (-20, 218.1, 86.0, 135.3)])
def test_resistance_table(Tc, R0m, R1m, tau):
    T = Tc + K0
    tol = 1e-3 if Tc == -20 else 2e-3  # §15.3 lists −20 °C at 1e-3; other rows are printed to 0.1
    assert rel(1e3 * float(cell.R0(0.5, T, RP63)), R0m) < tol + 0.05 / R0m
    assert rel(1e3 * float(cell.R1(T, RP63)), R1m) < tol + 0.05 / R1m
    assert rel(float(cell.tau1(T, RP63)), tau) < tol + 0.05 / tau


def test_lee20_plausibility_ratio():
    T5, T25 = -5 + K0, 25 + K0
    r = (cell.R0(0.5, T5, RP63) + cell.R1(T5, RP63)) / (cell.R0(0.5, T25, RP63) + cell.R1(T25, RP63))
    assert 2.4 < float(r) < 2.6  # "≈ 2.5"


# ------------------------------------------------------------------------------------------- power
PRIM_P = {"SB-W": 0.040, "SB-L": 0.047, "SB-5": 0.102, "RD": 0.456, "SO": 1.060, "SO-L": 1.218, "VW": 0.548,
          "VL": 0.694, "V5": 0.837, "NV": 1.486, "GM": 3.906, "GM-L": 4.302, "VC": 2.518, "VC-W": 1.622,
          "CM": 1.687, "AU": 0.224, "AU-L": 0.261, "TRK": 0.464, "BU": 4.095}


@pytest.mark.parametrize("name", list(PRIM_P))
def test_primitive_powers(name):
    assert abs(float(power.demand_power(PRIMITIVES[name], RP63)) - PRIM_P[name]) <= 0.5e-3 + 1e-9


def test_steady_game_power():
    assert abs(float(power.steady_power("GM", RP63, 25.0)) - 3.116) <= 1e-3
    assert abs(float(power.steady_power("GM", RP63, 35.0)) - 1.728) <= 1e-3
    d = analytic.thermal_cap("GM", RP63, 25.0)
    assert abs(float(d["rho_perf"]) - 0.864) < 6e-4 and abs(float(d["tte_cap_h"]) - 5.71) < 6e-3
    assert abs(float(d["sharp_cap_W"]) - 3.27) < 6e-3 and abs(float(d["P_nonsoc"]) - 0.743) < 6e-4
    d = analytic.thermal_cap("GM", RP63, 35.0)
    assert abs(float(d["rho_perf"]) - 0.531) < 6e-4 and abs(float(d["tte_cap_h"]) - 10.30) < 6e-3
    assert abs(float(d["sharp_cap_W"]) - 1.45) < 6e-3


# ------------------------------------------------------------------------------------------- §5.4 main table
MAIN = [  # (label, primitive, T_amb, TTE h, event, s_off %, T_end, Ts_end)
    ("S1", "SB-W", 25, 442.7, "cutoff", 0.02, 25.1, 25.2),
    ("S2a", "SB-L", 25, 378.3, "cutoff", 0.02, 25.1, 25.3),
    ("S2b", "SB-5", 25, 173.8, "cutoff", 0.04, 25.3, 25.6),
    ("S3", "RD", 25, 38.88, "cutoff", 0.18, 26.4, 27.7),
    ("S4", "SO", 25, 16.64, "cutoff", 0.40, 28.4, 31.2),
    ("S5a", "VW", 25, 32.32, "cutoff", 0.21, 26.7, 28.2),
    ("S5b", "VL", 25, 25.49, "cutoff", 0.27, 27.2, 29.0),
    ("S5c", "V5", 25, 21.10, "cutoff", 0.32, 27.7, 29.9),
    ("S6", "NV", 25, 11.83, "cutoff", 0.56, 29.8, 33.7),
    ("S7-25", "GM", 25, 5.68, "cutoff", 1.07, 34.5, 42.4),
    ("S7-35", "GM", 35, 10.29, "cutoff", 0.57, 40.3, 44.7),
    ("S7-15", "GM", 15, 4.40, "brownout", 1.89, 27.6, 38.2),
    ("S8a", "VC", 25, 6.93, "cutoff", 0.93, 33.1, 39.9),
    ("S8b", "VC-W", 25, 10.83, "cutoff", 0.61, 30.2, 34.5),
    ("CM", "CM", 25, 10.40, "brownout", 0.68, 30.4, 34.9),
    ("S16", "AU", 25, 79.20, "cutoff", 0.09, 25.7, 26.3),
]


@pytest.mark.parametrize("row", MAIN, ids=[r[0] for r in MAIN])
def test_main_table(row):
    label, prim, Ta, tte, ev, soff, Tend, Tsend = row
    r = lsoda(prim, float(Ta))
    assert rel(r.tte_h, tte) < 0.01
    assert r.event == ev
    assert abs(100 * r.s_off - soff) < 0.02
    assert abs(r.T_end_C - Tend) < 0.06 and abs(r.Ts_end_C - Tsend) < 0.06


def test_gaming_without_throttling():
    r = run(PRIMITIVES["GM"], 25.0, p=RP63.replace(burst=False, throttling=False))
    assert rel(r.tte_h, 4.43) < 0.01


@pytest.mark.parametrize("s0,tte", [(1.0, 34.92), (0.5, 11.79), (0.2, 8.70)])
def test_s12_day(s0, tte):
    r = run(s12_schedule(), 25.0, s0=s0)
    assert rel(r.tte_h, tte) < 0.01 and r.event == "cutoff"
    if s0 == 1.0:
        assert abs(100 * r.s_off - 0.93) < 0.02


def test_s12_mean_demand():
    P = s12_schedule().mean(lambda u: float(power.demand_power(u, RP63)))
    assert abs(P - 0.528) < 5e-4


# ------------------------------------------------------------------------------------------- S9 cold hike
S9 = [(25, 1.0, 26.99, 1.7), (25, 0.8, 21.32, 2.7), (-10, 1.0, 25.33, 7.2), (-10, 0.8, 17.00, 22.2),
      (-20, 1.0, 22.67, 16.9), (-20, 0.8, 12.33, 44.2)]


@pytest.mark.parametrize("Ta,soh,tte,soff", S9, ids=[f"{a}C-soh{b}" for a, b, _, _ in S9])
def test_s9_cold_hike(Ta, soh, tte, soff):
    r = run(s9_schedule(Ta), Ta, T0_C=20.0, p=aged(RP63, 1 - soh))
    assert r.event == "cutoff"
    assert rel(r.tte_h, tte) < 0.02
    assert abs(100 * r.s_off - soff) < 1.0


def test_s9_cold_soak_does_not_change_aged_result():
    r = run(s9_schedule(-10.0), -10.0, T0_C=-10.0, p=aged(RP63, 0.2))
    assert rel(r.tte_h, 17.00) < 0.02


# ------------------------------------------------------------------------------------------- stranded map
MAP = {25: [(0.4, 16.64), (0.5, 14.94), (0.9, 13.22)], 10: [(0.7, 16.57), (1.3, 14.80), (2.1, 13.05)],
       0: [(1.4, 16.43), (2.4, 14.62), (3.8, 12.80)], -10: [(2.8, 16.18), (4.9, 14.22), (8.5, 12.15)],
       -20: [(6.4, 15.54), (13.1, 12.99), (23.3, 10.23)]}
MAP_CASES = [(Ta, q, *MAP[Ta][i]) for Ta in MAP for i, q in enumerate((0.0, 0.1, 0.2))]


@pytest.mark.parametrize("Ta,q,soff,tte", MAP_CASES, ids=[f"{a}C-q{b}" for a, b, _, _ in MAP_CASES])
def test_stranded_charge_map(Ta, q, soff, tte):
    r = lsoda("SO", float(Ta), q)
    assert abs(100 * r.s_off - soff) <= 0.3
    assert rel(r.tte_h, tte) < 0.01
    assert r.event == ("cutoff" if (Ta == 25 and q == 0.0) else "brownout")


@pytest.mark.parametrize("m_b,tte,soff", [(1.0, 12.15, 8.5), (0.5, 12.88, 2.6), (0.0, 13.00, 1.5)])
def test_governor_O1(m_b, tte, soff):
    r = run(PRIMITIVES["SO"], -10.0, p=aged(RP63, 0.2).replace(m_b=m_b))
    assert rel(r.tte_h, tte) < 0.01 and abs(100 * r.s_off - soff) <= 0.3


# ------------------------------------------------------------------------------------------- other §5.4 numbers
@pytest.mark.parametrize("Ta,soh,tte,ev,soff", [(25, 0.8, 29.19, "cutoff", 0.3), (0, 1.0, 34.50, "brownout", 1.9),
                                                (0, 0.8, 24.78, "brownout", 4.7)])
def test_s10b_scheduled_day_aged(Ta, soh, tte, ev, soff):
    r = run(s12_schedule(), Ta, p=aged(RP63, 1 - soh))
    assert rel(r.tte_h, tte) < 0.01 and r.event == ev and abs(100 * r.s_off - soff) <= 0.3


@pytest.mark.parametrize("rsrp,t_sol,t_vc", [(-85, 14.46, 6.93), (-95, 13.72, 6.17), (-105, 13.05, 5.58),
                                             (-110, 12.73, 5.33), (-120, 12.15, 4.87)])
def test_s11_weak_signal(rsrp, t_sol, t_vc):
    r1 = run(PRIMITIVES["SO-L"].with_(RSRP=float(rsrp)), 25.0)
    r2 = run(PRIMITIVES["VC"].with_(RSRP=float(rsrp)), 25.0)
    assert rel(r1.tte_h, t_sol) < 0.01 and rel(r2.tte_h, t_vc) < 0.01


def test_navigation_cold_fresh():
    r = run(PRIMITIVES["NV"], -10.0)
    assert rel(r.tte_h, 11.43) < 0.01


# ------------------------------------------------------------------------------------------- V1 / V6
BONN = {("14 Pro", "idle"): 54.48, ("14 Pro", "cam"): 8.96, ("14 Pro", "gps"): 24.98, ("14 Pro", "cxl"): 7.48,
        ("13 Pro", "idle"): 51.85, ("13 Pro", "cam"): 8.52, ("13 Pro", "gps"): 23.77, ("13 Pro", "cxl"): 7.11}
BONN_BLIND = {("14 Pro", "idle"): -5.1, ("14 Pro", "cam"): 5.6, ("14 Pro", "gps"): -49.5, ("14 Pro", "cxl"): -11.3,
              ("13 Pro", "idle"): -2.2, ("13 Pro", "cam"): -11.2, ("13 Pro", "gps"): -21.7,
              ("13 Pro", "cxl"): -15.2}


def test_bonn_predictions():
    pred = V.bonn_predictions()
    for k, tte in BONN.items():
        assert rel(pred[k]["model_h"], tte) < 0.01, k
    assert abs(100 * V.bonn_additivity(pred) - 0.41) < 0.02
    obs = V.bonn_observations()
    assert abs(obs[("14 Pro", "idle")]["obs_h"] - 55.03) < 0.005


def test_bonn_blind_prior_and_calibration():
    pred = V.bonn_predictions(P_awake=0.10)
    for k, e in BONN_BLIND.items():
        assert abs(pred[k]["err_pct"] - e) < 0.6, (k, pred[k]["err_pct"])  # printed to 0.1 pp; cam 5.55 → 5.6
    assert abs(V.bonn_calibrate_awake() - 0.0877) < 1e-4


V6 = {"14P_local_video": 26.08, "14P_stream_video": 23.40, "14P_audio_BT": 55.90, "14PM_wifi_web_150": 15.18,
      "14PM_local_video_150": 33.01, "14PM_idle_minbright": 66.63, "14P_cell_web_150": 10.29,
      "Pixel8_wifi_web": 16.37, "Pixel8_local_video": 35.70, "Pixel8_idle_minbright": 70.52, "S24_wifi_web": 14.30,
      "S24_local_video": 31.20, "S24_idle_minbright": 61.64, "Pixel8a_S12_day": 34.67}


def test_v6_band():
    pred = V.v6_predictions()
    for k, tte in V6.items():
        assert rel(pred[k]["model_h"], tte) < 0.01, k
    n_pass = sum(v["pass_"] for v in pred.values())
    assert n_pass == 9  # 9 of these 14 model rows pass; the spec's '11 of 16' adds the GM-power and standby rows
    assert abs(100 / lsoda("SB-W", 25.0).tte_h - 0.23) < 0.005  # standby drain band row (S1, model run)


# ------------------------------------------------------------------------------------------- analytic
def test_F2_matches_isothermal_L2_quickly():
    p = RP63.replace(C_th=1e12, P_susp=0.0, burst=False)
    from phonebatt.usage import UsageState
    r = run(UsageState("c", P_extra=1.06), 25.0, p=p)
    assert rel(r.tte_h, float(analytic.tte_L1(1.06, p))) < 5e-4


def test_F3_a_coefficient():
    a = analytic.a_coefficient(RP63, 25.0)
    assert rel(float(a["total"]), 0.00818) < 0.02
    assert abs(float(a["ohm"]) - 0.00495) < 1e-5 and abs(float(a["cut"]) - 0.00322) < 1e-5
    assert abs(float(a["dUdz_ze"]) - 8.82) < 0.01


def test_F3_battery_loss_percentages():
    # L1 vs L0: −0.4 %, −0.9 %, −3.1 % at 0.5, 1.06, 3.5 W
    for P, loss in ((0.5, -0.4), (1.06, -0.9), (3.5, -3.1)):
        d = 100 * (float(analytic.tte_L1(P, RP63)) / float(analytic.tte_L0(P, RP63)) - 1)
        assert abs(d - loss) < 0.06


def test_F4_brownout_root():
    r = lsoda("SO", -20.0)
    assert abs(r.T_end_C - (-16.6)) < 0.06
    b = analytic.brownout_root(1.06, 3.0, r.T_end_C + K0, RP63)
    assert rel(float(b["z_star"]), 0.0881) < 0.01
    assert abs(100 * float(b["s_off"]) - 6.35) < 0.02 and abs(100 * r.s_off - 6.38) < 0.02


PCUT = {(25, 0.0): ((20.2, 12.5, 9.0), (28.3, 18.6, 13.9)), (-10, 0.0): ((6.9, 4.5, 3.4), (11.0, 7.4, 5.6)),
        (-10, 0.2): ((4.6, 3.0, 2.2), (7.3, 4.9, 3.7)), (-20, 0.2): ((3.2, 2.1, 1.6), (5.3, 3.6, 2.7))}


@pytest.mark.parametrize("key", list(PCUT))
def test_F4_pcut_table(key):
    Ta, q = key
    t = analytic.pcut_table(aged(RP63, q), Ta)
    np.testing.assert_allclose(t["P_cut"], PCUT[key][0], atol=0.051)
    np.testing.assert_allclose(t["P_cut_sys"], PCUT[key][1], atol=0.051)


def _typical():
    c = CTMC("typical")
    P = c.state_powers(RP63)
    return c, P, analytic.ctmc_moments(c.generator(), P, float(cell.usable_energy(RP63)))


def test_ctmc_typical_moments():
    c, P, m = _typical()
    assert rel(m["Pbar"], 0.4400) < 1e-3
    assert rel(m["sigma2"], 0.2492) < 1e-3
    assert rel(m["mu_dynkin"], 41.25) < 1e-3
    assert rel(m["sd_clt"], 7.21) < 1e-3
    assert abs(m["mu_fluid"] - 40.43) < 0.006 and abs(m["tau_c"] - 0.234) < 6e-4 and abs(m["cv_P"] - 1.659) < 6e-4
    assert abs(m["cv_tte"] - 0.175) < 6e-4
    q = analytic.ig_quantiles(m["mu_dynkin"], m["lam_IG"])
    np.testing.assert_allclose(q, [30.25, 40.59, 54.49], atol=0.006)
    assert abs(float(analytic.ig_cdf(24.0, m["mu_dynkin"], m["lam_IG"])) - 0.002) < 6e-4
    pia = m["pi"].reshape(-1, 2).sum(1)
    assert abs(24 * (1 - pia[0] - pia[-1]) - 6.5) < 0.06  # screen-on h/day
    share = {a: float(m["pi"].reshape(-1, 2)[i] @ P.reshape(-1, 2)[i]) / m["Pbar"] for i, a in enumerate(ACTS)}
    for a, v in dict(social=0.279, game=0.282, nav=0.106, video=0.100, audio=0.078, vcall=0.074, idle=0.058,
                     camera=0.024).items():
        assert abs(share[a] - v) < 6e-4, a


def test_ctmc_heavy_moments():
    c = CTMC("heavy")
    m = analytic.ctmc_moments(c.generator(), c.state_powers(RP63), float(cell.usable_energy(RP63)))
    assert rel(m["Pbar"], 0.8224) < 1e-3 and rel(m["sigma2"], 0.8388) < 1e-3
    assert abs(m["mu_dynkin"] - 22.50) < 0.006 and abs(m["sd_clt"] - 5.18) < 0.006
    assert abs(m["mu_fluid"] - 21.63) < 0.006 and abs(m["cv_tte"] - 0.230) < 6e-4


def test_ctmc_rate_scaling():
    c, P, m = _typical()
    L = c.generator()
    Eu = float(cell.usable_energy(RP63))
    for kap, mu, sd in ((0.25, 40.63, 3.61), (4.0, 43.71, 14.43)):
        mm = analytic.ctmc_moments(L / kap, P, Eu)
        assert abs(mm["mu_dynkin"] - mu) < 0.006 and abs(mm["sd_clt"] - sd) < 0.006
        assert abs(mm["sd_clt"] / m["sd_clt"] - math.sqrt(kap)) < 1e-10


@pytest.mark.parametrize("prof,s2", [("typical", 0.236495), ("heavy", 0.820905)])
def test_hub_closed_form_equals_deviation_matrix(prof, s2):
    lam, pj, dmin = PROFILES[prof]
    PaW = np.array([float(power.steady_power(MAP_WIFI[a], RP63)) for a in ACTS])
    hub = analytic.hub_closed_form(lam, [pj[a] for a in ACTS[1:]], [dmin[a] / 60 for a in ACTS[1:]], PaW[0],
                                   PaW[1:], ACTS[1:])
    mm = analytic.ctmc_moments(CTMC(prof).activity_generator(), PaW, 17.788)
    assert abs(hub["sigma2"] - mm["sigma2"]) < 1e-10
    assert abs(hub["sigma2"] - s2) < 1e-6
    # F9 with CV = 1 reduces to F8
    sm = analytic.semi_markov_sigma2(lam, [pj[a] for a in ACTS[1:]], [dmin[a] / 60 for a in ACTS[1:]], PaW[0],
                                     PaW[1:], 1.0, 1.0)
    assert abs(sm - hub["sigma2"]) < 1e-12
    if prof == "typical":
        for a, v in dict(game=0.807, nav=0.102, social=0.032, vcall=0.032, audio=0.021).items():
            assert abs(hub["share"][a] - v) < 6e-4
        for a, v in dict(game=-0.257, social=-0.169, audio=0.065).items():
            assert abs(hub["elast_d"][a] - v) < 6e-4
        assert abs(hub["elast_lam_on"] - (-0.531)) < 6e-4
        # lognormal sessions with CV 1.5 multiply the session terms by (1 + 2.25)/2 = 1.625
        sm15 = analytic.semi_markov_sigma2(lam, [pj[a] for a in ACTS[1:]], [dmin[a] / 60 for a in ACTS[1:]],
                                           PaW[0], PaW[1:], 2.25, 1.0)
        assert abs(sm15 / hub["sigma2"] - 1.625) < 1e-12
    else:
        assert abs(hub["share"]["game"] - 0.94) < 0.006


def test_two_state_formula():
    a, b, P1, P0 = 3.0, 0.7, 2.4, 0.05
    L = np.array([[-a, a], [b, -b]])  # off -> on at rate a, on -> off at rate b
    m = analytic.ctmc_moments(L, np.array([P0, P1]), 10.0)
    assert abs(m["sigma2"] - analytic.two_state_sigma2(a, b, P1, P0)) < 1e-10


def test_mc_L0_mean_vs_dynkin():
    c, P, m = _typical()
    x = analytic.mc_L0(c.generator(), P, float(cell.usable_energy(RP63)), n=20000, seed=1)
    assert rel(x.mean(), m["mu_dynkin"]) <= 0.01
    assert abs(x.std() / 7.04 - 1) < 0.05  # [PROTO] MC SD 7.04 h (seed-dependent)


@pytest.mark.slow
@pytest.mark.parametrize("prof,mean", [("typical", 39.53), ("heavy", 21.22)])
def test_mc_L2_mean(prof, mean):
    b = run_batch(CTMC(prof), 25.0, None, n=400, h=10.0, seed=7, t_max_h=120.0)
    assert rel(np.nanmean(b["tte_h"]), mean) < 0.03


# ------------------------------------------------------------------------------------------- aging
def test_g_stress_values():
    np.testing.assert_allclose(G.g_stress(np.array([0.1, 0.5, 0.7, 0.8, 0.9, 1.0]), RP63),
                               [0.702, 1.026, 1.724, 1.749, 2.043, 2.423], atol=6e-4)


def test_aging_calibration():
    kH, kL = G.calibrate(RP63)
    assert rel(kH, 2.697e-3) < 0.01 and rel(kL, 1.288e-3) < 0.01
    c = G.calibrate_details(RP63)
    assert abs(c["q_cal_N"] - 0.0782) < 6e-5 and abs(c["mean_g"] - 1.749) < 6e-4
    p = RP63.replace(k_H=c["k_H"], k_L=c["k_L"])
    st = G.AgingState()
    prof = G.day_profile("full", 25.0, 1.0, p, cal_mode=True)
    for _ in range(1000):
        st = G.advance_day(st, prof, p)
    assert abs(st.soh - 0.8) < 1e-9


AGING = {"A1": ((0.883, 0.834, 0.795), 2.86), "A2": ((0.897, 0.853, 0.818), math.inf),
         "A3": ((0.889, 0.841, 0.804), math.inf), "A4": ((0.796, 0.709, 0.641), 0.96),
         "A5": ((0.841, 0.774, 0.721), 1.58), "A6a": ((0.858, 0.796, 0.747), 1.93),
         "A6b": ((0.849, 0.783, 0.730), 1.72), "A6c": ((0.896, 0.853, 0.819), math.inf),
         "A7": ((0.872, 0.820, 0.779), 2.46)}


@pytest.mark.parametrize("sid", list(AGING))
def test_aging_scenarios(sid):
    r = G.run_years(years=3, p=RP63, **G.AGING_SCENARIOS[sid])
    soh = [r["soh"][365 * y - 1] for y in (1, 2, 3)]
    np.testing.assert_allclose(soh, AGING[sid][0], atol=5e-4 + 1e-9)
    if sid == "A1":
        assert abs(soh[1] - 0.834) <= 0.005  # §15.3 row
    yrs = AGING[sid][1]
    if math.isinf(yrs):
        assert math.isinf(r["years_to_80"])
    else:
        assert abs(r["years_to_80"] - yrs) < 6e-3


# ------------------------------------------------------------------------------------------- levers
LEVERS = {"skip_game": (34.9, 36.9), "dark_theme": (7.6, 4.6), "voice_call": (4.8, 4.0),
          "brightness_m01": (4.4, 3.1), "refresh_60": (4.1, 4.7), "vcall_wifi": (3.9, 1.6), "half_bg": (0.1, 0.4),
          "standby_5g": (-4.4, -4.5), "weak_signal": (-3.9, -5.7), "brightness_p02": (-10.3, -7.5),
          "os_coalesce": (0.1, 0.7), "os_frame_cap": (3.4, 2.7)}


def test_levers():
    L = analytic.lever_effects(RP63)
    assert abs(L["base"]["s12_tte_inf_h"] - 35.91) < 0.006 and abs(L["base"]["s12_Pbar_W"] - 0.4954) < 6e-5
    for k, (s12, ct) in LEVERS.items():
        assert abs(L[k]["s12_pct"] - s12) <= 0.3, k
        assert abs(L[k]["ctmc_pct"] - ct) <= 0.3, k
    assert abs(L["wifi_share_80"]["ctmc_pct"] - 3.7) <= 0.3
    assert abs(L["wifi_share_80"]["ctmc_mu_h"] - 42.78) < 0.006
    assert abs(float(power.steady_power(PRIMITIVES["GM"].with_(x_c=0.64), RP63)) - 2.722) < 6e-4
