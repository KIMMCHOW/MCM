"""AN4: analytic checks of the insight layer (FINAL_SPEC §7 F1–F5, §6.4 T3; §14).

Outputs
  results/analytic.json   F1 (E_u, E_n), F2_vs_L2_pct{P} (T3: isothermal, no bursts, constant P at 25 °C; plus 0 and
                          −20 °C as extra rows), L1_vs_L0_pct{P}, a_coeff{ohm, cut, total} (+ a vs T and SOH),
                          F4_check (S4 at −20 °C, plus the whole S10a stranded-charge map), safe_burst_W table,
                          F5_cap{Tamb} (vs the L2 ODE), Pcut_table (4 conditions × z = 0.5/0.2/0.1),
                          F2_steady_vs_L2_pct (F2 at the thermal steady state vs L2 for constant scenarios at 25 °C)
  paper/tables/tab_pcut.tex  cut-off- and burst-limited power table of F4 (optional, for §4 of the paper)

All L2 runs use the LSODA reference solver (§6.2).
"""
from __future__ import annotations

import time

import numpy as np

from basics_common import f, meta, write_json, write_tex
from phonebatt import RP63, aged, analytic, cell, power, run, scenario
from phonebatt.params import K0
from phonebatt.usage import PRIMITIVES, UsageState

P_T3 = (0.5, 1.06, 1.49, 2.52, 3.5)
PROTO = dict(E_u=17.788, E_n=18.189, a_ohm=0.00495, a_cut=0.00322, a_total=0.00818, dUdz_ze=8.82,
             L1_vs_L0={0.5: -0.4, 1.06: -0.9, 3.5: -3.1}, F4_z_star=0.0881, F4_s_off=6.35, F4_ode_s_off=6.38,
             F4_T_end=-16.6,
             F5={15: dict(P_star=3.902, ode_tte=4.40), 25: dict(P_star=3.116, tte_cap=5.71, ode_tte=5.68, rho=0.864,
                                                              sharp=3.27),
                 35: dict(P_star=1.728, tte_cap=10.30, ode_tte=10.29, rho=0.531, sharp=1.45)},
             P_nonsoc_GM=0.743,
             Pcut={"25|1.0": ([20.2, 12.5, 9.0], [28.3, 18.6, 13.9]), "-10|1.0": ([6.9, 4.5, 3.4], [11.0, 7.4, 5.6]),
                   "-10|0.8": ([4.6, 3.0, 2.2], [7.3, 4.9, 3.7]), "-20|0.8": ([3.2, 2.1, 1.6], [5.3, 3.6, 2.7])},
             stranded={  # §5.4 S10a map: s_off % (TTE h)
                 (25, 1.0): (0.4, 16.64), (25, 0.9): (0.5, 14.94), (25, 0.8): (0.9, 13.22),
                 (10, 1.0): (0.7, 16.57), (10, 0.9): (1.3, 14.80), (10, 0.8): (2.1, 13.05),
                 (0, 1.0): (1.4, 16.43), (0, 0.9): (2.4, 14.62), (0, 0.8): (3.8, 12.80),
                 (-10, 1.0): (2.8, 16.18), (-10, 0.9): (4.9, 14.22), (-10, 0.8): (8.5, 12.15),
                 (-20, 1.0): (6.4, 15.54), (-20, 0.9): (13.1, 12.99), (-20, 0.8): (23.3, 10.23)})


def iso_params(p):
    """T3 configuration: isothermal (C_th = 1e12), no bursts, no floor, so P is exactly the constant load."""
    return p.replace(C_th=1e12, P_susp=0.0, burst=False)


def main() -> None:
    t0 = time.time()
    out = dict(meta=meta("run_analytic.py", solver="LSODA rtol 1e-7 (L2); Gauss-Legendre n=200 + bisection (F2)",
                         h_s=None))
    Eu, En = f(cell.usable_energy(RP63)), f(cell.nominal_energy(RP63))
    out["F1"] = dict(E_u_Wh=Eu, E_n_Wh=En, z_e=f(cell.z_empty(RP63)), proto=dict(E_u_Wh=17.788, E_n_Wh=18.189))

    # ------------------------------------------------------------------------------- F2 vs L2 (T3)
    f2 = {}
    f2T = {}
    for Tc in (25.0, 0.0, -20.0):
        pi = iso_params(RP63)
        rows = {}
        for P in P_T3:
            r = run(UsageState("const", P_extra=P), Tc, p=pi)
            tl1 = f(analytic.tte_L1(P, pi, Tc))
            rows[f"{P}"] = dict(L2_h=r.tte_h, F2_h=tl1, diff_pct=100 * (tl1 / r.tte_h - 1), event=r.event)
        if Tc == 25.0:
            f2 = {k: v["diff_pct"] for k, v in rows.items()}
            out["F2_vs_L2_detail"] = rows
        f2T[f"{Tc:g}"] = rows
    out["F2_vs_L2_pct"] = f2
    out["F2_vs_L2_pct_T"] = {T: {k: v["diff_pct"] for k, v in rows.items()} for T, rows in f2T.items()}
    out["T3"] = dict(max_abs_pct=max(abs(v) for v in f2.values()), tol_pct=0.05, proto_max_pct=0.001,
                     pass_=max(abs(v) for v in f2.values()) < 0.05)
    # L1 vs L0 (battery losses at 25 °C, F3 check)
    a = analytic.a_coefficient(RP63, 25.0)
    out["L1_vs_L0_pct"] = {f"{P}": dict(exact=100 * (f(analytic.tte_L1(P, RP63)) / f(analytic.tte_L0(P, RP63)) - 1),
                                        first_order=-100 * f(a["total"]) * P,
                                        proto=PROTO["L1_vs_L0"].get(P)) for P in P_T3}

    # ------------------------------------------------------------------------------- F3 a coefficient
    out["a_coeff"] = dict(ohm=f(a["ohm"]), cut=f(a["cut"]), total=f(a["total"]), dUdz_ze=f(a["dUdz_ze"]),
                          unit="1/W", proto=dict(ohm=PROTO["a_ohm"], cut=PROTO["a_cut"], total=PROTO["a_total"],
                                                 dUdz_ze=PROTO["dUdz_ze"]),
                          dev_total_pct=100 * (f(a["total"]) / PROTO["a_total"] - 1))
    aT = {}
    for soh in (1.0, 0.8):
        pp = aged(RP63, 1 - soh)
        for Tc in (25, 10, 0, -10, -20):
            ac = analytic.a_coefficient(pp, float(Tc))
            aT[f"{Tc}|{soh}"] = dict(ohm=f(ac["ohm"]), cut=f(ac["cut"]), total=f(ac["total"]),
                                     ratio_to_25C_new=f(ac["total"]) / f(a["total"]))
    out["a_coeff_T_soh"] = aT

    # ------------------------------------------------------------------------------- F4 brown-out
    r = scenario("S4", Tamb_C=-20.0).run()
    b = analytic.brownout_root(1.06, 3.0 * f(RP63.m_b), r.T_end_C + K0, RP63)
    out["F4_check"] = dict(case="S4 social, fresh, T_amb = -20 C, battery T from the L2 run",
                           T_end_C=r.T_end_C, z_star=f(b["z_star"]), s_off_F4_pct=100 * f(b["s_off"]),
                           z_end_L2=r.z_end, s_off_L2_pct=100 * r.s_off, event_L2=r.event,
                           proto=dict(z_star=PROTO["F4_z_star"], s_off_F4_pct=PROTO["F4_s_off"],
                                      s_off_L2_pct=PROTO["F4_ode_s_off"], T_end_C=PROTO["F4_T_end"]),
                           dev_z_star_pct=100 * (f(b["z_star"]) / PROTO["F4_z_star"] - 1))
    smap = {}
    for (Ta, soh), (so_p, tte_p) in PROTO["stranded"].items():
        pp = aged(RP63, 1 - soh)
        rr = scenario("S4", Tamb_C=float(Ta), p=pp).run()
        bb = analytic.brownout_root(rr.P_end_W, 3.0, rr.T_end_C + K0, pp)
        smap[f"{Ta}|{soh}"] = dict(event=rr.event, tte_h=rr.tte_h, s_off_L2_pct=100 * rr.s_off,
                                   s_off_F4_pct=100 * f(bb["s_off"]), T_end_C=rr.T_end_C,
                                   F4_minus_L2_pp=100 * (f(bb["s_off"]) - rr.s_off),
                                   proto_s_off_pct=so_p, proto_tte_h=tte_p,
                                   dev_s_off_pp=100 * rr.s_off - so_p, dev_tte_pct=100 * (rr.tte_h / tte_p - 1))
    br = [v for v in smap.values() if v["event"] == "brownout"]
    out["F4_map_check"] = dict(cells=smap, n_brownout=len(br),
                               max_abs_F4_minus_L2_pp_brownout=max(abs(v["F4_minus_L2_pp"]) for v in br),
                               max_abs_dev_s_off_vs_proto_pp=max(abs(v["dev_s_off_pp"]) for v in smap.values()),
                               max_abs_dev_tte_vs_proto_pct=max(abs(v["dev_tte_pct"]) for v in smap.values()))
    zs = np.array([0.5, 0.2, 0.1])
    sb = {}
    for Tc in (25, 0, -10, -20):
        for soh in (1.0, 0.8):
            pp = aged(RP63, 1 - soh)
            sb[f"{Tc}|{soh}"] = dict(z=zs, dPb_safe_W=analytic.safe_burst(zs, Tc + K0, 1.06, pp))
    out["safe_burst_W"] = dict(load="S4 social mean load 1.06 W, v1 = R1 P/U", table=sb)

    # ------------------------------------------------------------------------------- F5 thermal cap
    f5 = {}
    for Ta in (5, 15, 25, 30, 35, 40):
        d = analytic.thermal_cap("GM", RP63, float(Ta))
        rr = scenario("S7", Tamb_C=float(Ta)).run()
        e = dict(P_star_W=f(d["P_star"]), Ts_star_C=f(d["Ts_star_C"]), rho_perf=f(d["rho_perf"]),
                 tte_cap_h=f(d["tte_cap_h"]), tte_cap_L1_h=f(d["tte_cap_L1_h"]), sharp_cap_W=f(d["sharp_cap_W"]),
                 dTTE_dTamb_h_per_K=f(d["dTTE_dTamb"]), P_nonsoc_W=f(d["P_nonsoc"]), P_demand_W=f(d["P_demand"]),
                 ode_tte_h=rr.tte_h, ode_event=rr.event, ode_rho_perf=rr.rho_perf, ode_P_mean_W=rr.P_mean_W,
                 tte_cap_vs_ode_pct=100 * (f(d["tte_cap_h"]) / rr.tte_h - 1))
        if Ta in PROTO["F5"]:
            e["proto"] = PROTO["F5"][Ta]
            e["dev_P_star_mW"] = 1e3 * (e["P_star_W"] - PROTO["F5"][Ta]["P_star"])
        f5[f"{Ta}"] = e
    out["F5_cap"] = f5
    out["F5_dTTE_dTamb_ode_25_35"] = (f5["35"]["ode_tte_h"] - f5["25"]["ode_tte_h"]) / 10.0

    # ------------------------------------------------------------------------------- P_cut table
    pc = {}
    for key in PROTO["Pcut"]:
        Ta, soh = key.split("|")
        t = analytic.pcut_table(aged(RP63, 1 - float(soh)), float(Ta))
        pp_, ps_ = PROTO["Pcut"][key]
        pc[key] = dict(z=t["z"], P_cut_W=t["P_cut"], P_cut_sys_W=t["P_cut_sys"], proto_P_cut_W=pp_,
                       proto_P_cut_sys_W=ps_,
                       max_abs_dev_W=float(max(np.max(np.abs(t["P_cut"] - pp_)), np.max(np.abs(t["P_cut_sys"] - ps_)))))
    out["Pcut_table"] = pc

    # ------------------------------------------------------------------------------- F2 steady vs L2 (25 °C)
    st = {}
    for sc, prim in (("S3", "RD"), ("S4", "SO"), ("S5a", "VW"), ("S5b", "VL"), ("S5c", "V5"), ("S6", "NV"),
                     ("S7", "GM"), ("S8a", "VC"), ("S8b", "VC-W"), ("S16", "AU"), ("CM", "CM")):
        u = PRIMITIVES[prim]
        Pst = f(power.steady_power(u, RP63, 25.0))
        Tb = 25.0 + f(RP63.k_env * RP63.R_th) * Pst
        rr = scenario(sc).run()
        tl1 = f(analytic.tte_L1(Pst, RP63, Tb))
        st[sc] = dict(P_W=Pst, T_batt_C=Tb, F2_h=tl1, L0_h=Eu / Pst, L2_h=rr.tte_h, L2_event=rr.event,
                      F2_vs_L2_pct=100 * (tl1 / rr.tte_h - 1), L0_vs_L2_pct=100 * (Eu / Pst / rr.tte_h - 1))
    out["F2_steady_vs_L2_pct"] = st
    write_json("analytic.json", out)
    write_pcut_tex(pc)

    # ------------------------------------------------------------------------------- report
    print("F2 vs L2 (T3, 25 C):", {k: round(v, 5) for k, v in f2.items()}, "% ; 0/-20 C:",
          {T: {k: round(v, 4) for k, v in rows.items()} for T, rows in out["F2_vs_L2_pct_T"].items() if T != "25"})
    print("L1 vs L0 %:", {k: (round(v["exact"], 3), round(v["first_order"], 3)) for k, v in out["L1_vs_L0_pct"].items()})
    print(f"a = {out['a_coeff']['ohm']:.5f} + {out['a_coeff']['cut']:.5f} = {out['a_coeff']['total']:.5f} 1/W "
          f"(PROTO 0.00818); U'(z_e) = {out['a_coeff']['dUdz_ze']:.3f}")
    print("a ratios:", {k: round(v["ratio_to_25C_new"], 2) for k, v in aT.items()})
    fc = out["F4_check"]
    print(f"F4: z*={fc['z_star']:.4f} s_off F4 {fc['s_off_F4_pct']:.2f} % vs L2 {fc['s_off_L2_pct']:.2f} % "
          f"(T_end {fc['T_end_C']:.1f} C)")
    fm = out["F4_map_check"]
    print(f"F4 map: {fm['n_brownout']} brown-out cells, max |F4 - L2| = {fm['max_abs_F4_minus_L2_pp_brownout']:.3f} pp;"
          f" vs PROTO max |ds_off| {fm['max_abs_dev_s_off_vs_proto_pp']:.3f} pp, max |dTTE| "
          f"{fm['max_abs_dev_tte_vs_proto_pct']:.3f} %")
    for Ta, e in f5.items():
        print(f"F5 {Ta} C: P*={e['P_star_W']:.4f} rho={e['rho_perf']:.3f} cap {e['tte_cap_h']:.2f} h (L1 "
              f"{e['tte_cap_L1_h']:.2f}) vs ODE {e['ode_tte_h']:.2f} h {e['ode_event']} (rho {e['ode_rho_perf']:.3f}) "
              f"sharp {e['sharp_cap_W']:.2f}")
    for k, v in pc.items():
        print(f"Pcut {k}: {np.round(v['P_cut_W'], 2)} / {np.round(v['P_cut_sys_W'], 2)}  max dev {v['max_abs_dev_W']:.3f} W")
    for k, v in st.items():
        print(f"steady {k:4s}: P {v['P_W']:.3f} F2 {v['F2_h']:.2f} L0 {v['L0_h']:.2f} L2 {v['L2_h']:.2f} "
              f"({v['L2_event']}) F2-L2 {v['F2_vs_L2_pct']:+.2f} % L0-L2 {v['L0_vs_L2_pct']:+.2f} %")
    print(f"done in {time.time() - t0:.1f} s")


def write_pcut_tex(pc: dict) -> None:
    label = {"25|1.0": "25\\,\\textdegree C, new", "-10|1.0": "$-$10\\,\\textdegree C, new",
             "-10|0.8": "$-$10\\,\\textdegree C, SOH 0.8", "-20|0.8": "$-$20\\,\\textdegree C, SOH 0.8"}
    rows = []
    for k, v in pc.items():
        a = " / ".join(f"{x:.1f}" for x in v["P_cut_W"])
        b = " / ".join(f"{x:.1f}" for x in v["P_cut_sys_W"])
        rows.append(f"{label[k]} & {a} & {b} \\\\")
    txt = "\n".join([
        "% F4 power limits (generated by analysis/run_analytic.py). P_cut = V_cut(U - V_cut)/(R0 + R1): loads >= 3 tau1;",
        "% P_cut,sys = V_sys(U - V_sys)/R0: sub-second peaks. Battery at the stated temperature, v1 = 0.",
        "\\small",
        "\\begin{tabular}{@{}lcc@{}}",
        "\\toprule",
        "Condition & $P_{\\mathrm{cut}}$ (W) at $z=0.5/0.2/0.1$ & $P_{\\mathrm{cut,sys}}$ (W) at $z=0.5/0.2/0.1$ \\\\",
        "\\midrule",
        *rows,
        "\\bottomrule",
        "\\end{tabular}",
        "",
    ])
    write_tex("tab_pcut.tex", txt)


if __name__ == "__main__":
    main()
