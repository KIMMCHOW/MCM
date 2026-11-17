"""Model-level tests: batch vectorisation, events, scenarios, loaders."""
from __future__ import annotations

import math

import numpy as np
import pytest

from phonebatt import (CTMC, PRIMITIVES, RP63, Schedule, aged, run, run_batch, sample, scenario, with_theta)
from phonebatt import validation as V
from phonebatt.params import iphone14pro
from phonebatt.usage import SCENARIO_IDS, UsageState, wifi_cliff


def test_batch_samples_equal_individual_runs():
    """Per-sample θ in one batch gives the same TTE as separate single-sample batches (vectorisation)."""
    th = sample(["R0_ref", "Ea_R0", "C_th", "R_skin", "P_SoC_max", "m_b", "V_cut"], 6, seed=4)
    sch = Schedule.constant(PRIMITIVES["GM"])
    Ta = np.array([25.0, 10.0, 0.0, 35.0, -5.0, 20.0])
    b = run_batch(sch, Ta, th, h=10.0)
    for k in range(6):
        bk = run_batch(sch, Ta[k], {n: v[k:k + 1] for n, v in th.items()}, h=10.0)
        assert b["tte_h"][k] == pytest.approx(bk["tte_h"][0], rel=1e-12)
        assert b["event"][k] == bk["event"][0]
        assert b["s_off"][k] == pytest.approx(bk["s_off"][0], rel=1e-9)


def test_batch_rk4_matches_lsoda_on_short_runs():
    for prim, Ta, q in (("GM", 25.0, 0.0), ("VC", 0.0, 0.2)):
        p = aged(RP63, q)
        r = run(PRIMITIVES[prim], Ta, p=p)
        b = run_batch(Schedule.constant(PRIMITIVES[prim]), Ta, None, h=2.0, p=p)
        assert abs(b["tte_h"][0] - r.tte_h) / r.tte_h < 2e-3
        assert b["event"][0] == r.event
        assert abs(b["s_off"][0] - r.s_off) < 2e-3


def test_thermal_shutdown_event():
    p = RP63.replace(throttling=False, burst=False, k_env=1.5)
    r = run(PRIMITIVES["GM"], 50.0, p=p)
    assert r.event == "thermal" and abs(r.T_end_C - 68.0) < 1e-6
    b = run_batch(Schedule.constant(PRIMITIVES["GM"]), 50.0, None, h=2.0, p=p)
    assert b["event"][0] == "thermal" and abs(b["tte_h"][0] - r.tte_h) < 2 * 2.0 / 3600


def test_switch_time_event_at_camera_start_in_cold():
    """A burst segment that is infeasible for the cut-off at the switch ends the run at the switch time."""
    sch = Schedule([(3600.0 * 9, PRIMITIVES["TRK"]), (600.0, PRIMITIVES["BU"])])
    p = aged(RP63, 0.2)
    r = run(sch, -20.0, T0_C=-20.0, p=p)
    assert r.event == "cutoff"
    t = r.tte_h * 3600.0
    k = round((t - 9 * 3600.0) / sch.period)
    assert k >= 0 and abs(t - (9 * 3600.0 + k * sch.period)) < 1e-6  # exactly at a TRK -> BU switch
    assert r.s_off > 0.2  # a large stranded charge: the burst cannot be served
    b = run_batch(sch, -20.0, None, T0_C=-20.0, p=p, h=10.0)
    assert b["event"][0] == "cutoff" and abs(b["tte_h"][0] - r.tte_h) <= 10.0 / 3600 + 1e-9


def test_horizon_returns_infinite_tte():
    u = UsageState("tiny", P_extra=1e-4)
    r = run(u, 25.0, p=RP63.replace(P_susp=0.0), t_max_h=10.0)
    assert math.isinf(r.tte_h) and r.event == "none"
    b = run_batch(Schedule.constant(u), 25.0, None, p=RP63.replace(P_susp=0.0), h=10.0, t_max_h=10.0)
    assert math.isinf(b["tte_h"][0]) and b["event"][0] == "none"


def test_initial_soc_and_outputs():
    r = run(PRIMITIVES["SO"], 25.0, s0=0.5)
    r1 = run(PRIMITIVES["SO"], 25.0)
    assert 0.45 < r.tte_h / r1.tte_h < 0.52
    assert 0 < r.soe_off < r.s_off  # energy fraction below the coulombic display near empty
    assert abs(r.P_mean_W - 1.0599) < 1e-3
    assert math.isnan(r.t_disp0_h)  # the display never reaches 0 before the cut-off at 25 °C
    # a cold, low load: OCV is raised by the entropic term, the display reaches 0 % before shutdown
    rc = run(PRIMITIVES["SB-W"], -20.0, s0=0.05)
    assert rc.s_off < 0 and rc.t_disp0_h < rc.tte_h


@pytest.mark.parametrize("sid", SCENARIO_IDS)
def test_scenario_constructors(sid):
    sc = scenario(sid)
    assert sc.schedule is not None
    if sid in ("S13", "S14"):
        assert isinstance(sc.schedule, CTMC)
    else:
        assert len(sc.schedule) >= 1
    assert scenario("S9", Tamb_C=-10.0, soh=0.8).T0_C == 20.0
    assert scenario("S10", Tamb_C=-10.0, q_loss=0.2).p.rg == pytest.approx(0.5)


def test_wifi_cliff_and_s11b():
    assert wifi_cliff(-85.0) == pytest.approx(2.13)
    assert wifi_cliff(-90.0) == pytest.approx(1 + 1.13 * math.exp(5 / 2.54))
    r0 = run(PRIMITIVES["SO"], 25.0)
    r = scenario("S11", variant="b", RSSI=-90.0).run()
    assert r.tte_h < r0.tte_h


def test_ctmc_5g_variant_and_state_map():
    c = CTMC("typical", rat="nr")
    us = c.primitives()
    assert all(u.c in ("nr", "none") for u in us)
    assert c.state_map[0] == "SB-W" and c.state_map[1] == "SB-L"
    assert CTMC("typical").state_powers(RP63)[0] < c.state_powers(RP63)[1]


def test_loaders():
    obs = V.bonn_observations()
    assert len(obs) == 8 and obs[("13 Pro", "gps")]["obs_h"] == pytest.approx(29.76, abs=0.01)
    for T in V.HG2_TESTS:
        d = V.load_hg2(T, "Dis_0p5C")
        assert d["I_A"].mean() > 1.4 and d["V"][-1] < 2.9
    s = V.load_hg2_summary()
    q = {(r["T_amb_C"], r["kind"]): r["Q_Ah"] for r in s["results"]}
    assert V.load_hg2(25, "C20DisCh")["Ah"][-1] == pytest.approx(q[(25, "C20DisCh")], rel=2e-3)
    n = V.load_nasa()
    assert n["I_A"].max() > 1.9
    ox = V.load_oxford()
    assert ox["I_A"].max() == pytest.approx(5.0, abs=0.01)
    cal = V.hg2_calibration()
    assert cal["Q_Ah"] == pytest.approx(2.781, abs=1e-3) and cal["rmse_T_K"] < 1.0
    p = V.hg2_params(cal)
    assert p.ocv_table is not None and p.R0_ref + p.R1_ref == pytest.approx(0.038)
