"""Regression tests for defects found in the adversarial review of 2026-11-17 ([REVIEW 2026-11-17] in the code
and in FINAL_SPEC).  Each test failed on the reviewed commit 911187e and passes after the fix."""
from __future__ import annotations

import math

import numpy as np

from phonebatt import RP63, PRIMITIVES, aged, run, run_batch, sample
from phonebatt import validation as V
from phonebatt.params import hg2_cell, nasa_b0005
from phonebatt.usage import Schedule, UsageState, s9_schedule


def test_batch_detects_crossing_in_last_step_before_lighter_segment():
    """A smooth event crossing inside the last RK4 step of a heavy segment must be detected even though the
    next segment is lighter (the event function was only evaluated with the new input: +25 % TTE here)."""
    p = aged(RP63, 0.2)
    heavy, light = PRIMITIVES["SO"], PRIMITIVES["SB-W"]
    tc = run(heavy, -10.0, s0=0.3, p=p).tte_h * 3600.0
    h = 10.0
    d1 = math.ceil(tc / h) * h
    assert 0.0 < d1 - tc < h  # the crossing lies inside the last step of the heavy segment
    sch = Schedule([(d1, heavy), (3 * 3600.0, light)])
    r = run(sch, -10.0, s0=0.3, p=p)
    b = run_batch(sch, -10.0, None, s0=0.3, h=h, p=p)
    assert r.event == b["event"][0] == "brownout"
    assert abs(b["tte_h"][0] - r.tte_h) * 3600.0 < h
    assert abs(b["s_off"][0] - r.s_off) < 1e-3


def test_batch_s9_parameter_sample_matches_lsoda():
    """S9 (−10 °C, SOH 0.8) parameter sample whose cut-off crossing falls in the last 2 s of a BU burst: the
    batch result was one full 600 s period late (+2.4 %)."""
    th = sample(["R0_ref", "Ea_R0", "R1_ref", "tau1_ref", "V_cut", "Q_n", "C_th"], 48, seed=99)
    k = 45
    p = aged(RP63, 0.2).replace(**{n: float(v[k]) for n, v in th.items()})
    sch = s9_schedule(-10.0)
    r = run(sch, -10.0, T0_C=20.0, p=p)
    b = run_batch(sch, -10.0, None, h=2.0, T0_C=20.0, p=p)
    assert b["event"][0] == r.event == "cutoff"
    assert abs(b["tte_h"][0] - r.tte_h) / r.tte_h < 2e-3  # T1 tolerance


def test_batch_switch_event_is_at_switch_time():
    """An event caused by the input switch itself (camera-and-upload burst in the cold) happens at the switch
    instant (§6.1); the batch interpolated it up to one step earlier."""
    sch = Schedule([(3600.0 * 9, PRIMITIVES["TRK"]), (600.0, PRIMITIVES["BU"])])
    p = aged(RP63, 0.2)
    r = run(sch, -20.0, T0_C=-20.0, p=p)
    b = run_batch(sch, -20.0, None, T0_C=-20.0, p=p, h=10.0)
    assert b["event"][0] == r.event == "cutoff"
    assert abs(b["tte_h"][0] - r.tte_h) * 3600.0 < 1e-6


def test_cell_only_devices_do_not_heat_with_external_load():
    """A bare cell on a tester: the delivered power is dissipated outside the cell (A4 applies to phones only).
    With the phone heat balance the calibrated HG2 cell reached the 68 °C shutdown at 2C."""
    cal = V.hg2_calibration()
    p = V.hg2_params(cal)
    d = V.load_hg2(25, "Dis_2C")
    P = float(np.mean(d["I_A"] * d["V"]))
    r = run(UsageState("cc2C", P_extra=P), 25.0, p=p)
    assert r.event == "cutoff"
    rise_obs = float(d["T_C"].max() - d["T_C"][0])
    assert abs((r.T_end_C - 25.0) - rise_obs) < 2.0
    assert abs(r.tte_h * 3600.0 - d["t_s"][-1]) / d["t_s"][-1] < 0.05
    assert not hg2_cell().load_heat and not nasa_b0005().load_heat and RP63.load_heat
