# Group "global" (AN9 Sobol + AN10 structural variants): digest for the paper writer

Scripts: `analysis/run_sobol.py` (AN9, 625 s with 2 workers), `analysis/run_structural.py` (AN10, 594 s) with `analysis/structural_variants.py` (union model) and this digest (`analysis/summary_global.py`). Outputs: `results/sobol.json`, `results/structural.json`, `fig_sobol`, `fig_structural`, `paper/tables/tab_structural.tex` (T4). Core untouched; no core bug found. Seeds fixed (2026); RP-63, s0 = 100 % unless stated.

## AN9 Sobol (§10.3): what drives TTE
- Design: 24 factor groups, SALib Saltelli design (second order), S_i by Saltelli 2010, S_Ti by Jansen; 95 % bootstrap CIs (1000 resamples). Surrogate scenarios N = 1024 (51,200 evals each); full L2 (rk4, h = 2 s) N = 256 (12,800 runs); S10a-G N = 512 (k = 26, 27,648 runs).
- Surrogate (F2 + F4 + F5 incl. q_batt and the thermal event) vs full L2 on 200 points each: max |err| 0.77 % over the scenarios that use it (S3, S4, S5a, S5b, S5c, S6, S8a, S8b, S11a; spec: ≤ 1 %). S11c fails (6.9 %, 3/200 points, all 7–10 W calls at the 68 °C thermal-shutdown boundary where TTE jumps), so S11c is reported from full L2 (N = 256); surrogate kept as `S11c_surrogate`.

| Scenario | method | mean TTE (h) | CV | top S_Ti (S_i) | ΣS_i |
|---|---|---|---|---|---|
| S3 | surr. | 30.84 | 0.19 | display 0.68 (0.66), floor 0.29 (0.26), SoC 0.03 (0.03) | 0.98 |
| S4 | surr. | 14.85 | 0.14 | SoC 0.50 (0.49), display 0.24 (0.23), floor 0.12 (0.12) | 0.99 |
| S5a | surr. | 25.04 | 0.16 | display 0.31 (0.30), floor 0.26 (0.25), Wi-Fi 0.17 (0.16) | 0.97 |
| S5b | surr. | 20.69 | 0.14 | LTE 0.29 (0.28), display 0.26 (0.25), floor 0.23 (0.22) | 0.98 |
| S5c | surr. | 15.23 | 0.13 | NR 0.54 (0.53), display 0.16 (0.16), floor 0.14 (0.14) | 0.97 |
| S6 | surr. | 10.45 | 0.15 | display 0.61 (0.59), SoC 0.22 (0.21), peripherals 0.07 (0.06) | 0.96 |
| S8a | surr. | 5.74 | 0.21 | LTE 0.87 (0.85), SoC 0.05 (0.02), peripherals 0.03 (0.02) | 0.96 |
| S8b | surr. | 8.82 | 0.15 | Wi-Fi 0.47 (0.45), SoC 0.25 (0.20), peripherals 0.12 (0.11) | 0.95 |
| S11a | surr. | 10.67 | 0.18 | γ_sig 0.37 (0.35), LTE 0.32 (0.30), SoC 0.18 (0.16) | 0.96 |
| S11c | L2 | 4.03 | 0.34 | LTE 0.62 (0.64), γ_sig 0.32 (0.33), SoC 0.01 (-0.00) | 1.00 |
| S7 | L2 | 6.69 | 0.19 | k_env 0.40 (0.38), R_skin 0.38 (0.35), throttle 0.18 (0.14) | 0.93 |
| S9 | L2 | 9.28 | 0.56 | R0,ref 0.47 (0.45), γ_R 0.14 (0.09), Ea,R0 0.08 (0.07) | 0.93 |
| S10a | L2 | 7.39 | 0.47 | R0,ref 0.47 (0.42), m_b 0.25 (0.22), γ_R 0.13 (0.08) | 0.95 |
| S10a-G | L2 | 9.98 | 0.39 | R0,ref 0.30 (0.24), SOH 0.23 (0.19), T_amb 0.19 (0.13) | 0.93 |

- S10a-G (social; T_amb ~ U(−20, 10) °C, SOH ~ U(0.75, 1), m_b explicit): TTE S_i/S_Ti T_amb 0.13/0.19, SOH 0.19/0.23, m_b 0.09/0.14, R0,ref 0.24/0.30; interaction share 1 − ΣS_i = 0.07 [-0.07, 0.21]; S_ij T_amb×SOH -0.006 ± 0.078, T_amb×m_b +0.019 ± 0.078, SOH×m_b +0.043 ± 0.086.
- Same for s_off (stranded charge, pp): T_amb 0.19/0.28, SOH 0.07/0.14, m_b 0.13/0.22, R0,ref 0.35/0.43; interaction share 0.12 [-0.04, 0.26]; S_ij T_amb×SOH -0.006, T_amb×m_b +0.034, SOH×m_b +0.033.
- S9 s_off (mean 53.2, median 54.3 pp): top S_Ti (S_i) R0,ref 0.55 (0.53), γ_R 0.16 (0.10), Ea,R0 0.09 (0.08).
- S10a s_off (mean 39.2, median 34.6 pp): top S_Ti (S_i) R0,ref 0.50 (0.46), m_b 0.27 (0.24), γ_R 0.14 (0.09).
- Reading. 25 °C constant loads: variance is additive (ΣS_i ≥ 0.95) and carried by the load parameters of the active components; every cell/thermal group has S_Ti ≤ 0.02 (an energy-budget problem; Q_n is small only because its range is ±2 %). Gaming: thermal design (k_env, R_skin, throttle). Cold + aged: R0,ref, γ_R, E_a,R0 and the burst size m_b, with TTE CV 0.47–0.56 vs 0.13–0.21 at 25 °C (a voltage problem). Same three regimes as §10.1. Interactions in S10a-G are modest and their CIs include 0.

## AN10 structural variants (§10.4): ΔTTE (%) vs X0, same solver and h
- X0 of the union model vs core LSODA on all 23 deterministic scenarios: max |ΔTTE| 0.0001 %, events identical: True; vs core run_batch (same h): |Δ| = 0.0e+00 h, 0.0e+00 h. Variant effects at h = 1 s vs 2 s differ by ≤ 0.0000 pp.

| X | variant | S4 | S7 | S7 35°C | S12 | S13 | S14 | S9 −20/0.8 | S10a −20/0.8 | Δs_off S10a −20/0.8 (pp) |
|---|---|---|---|---|---|---|---|---|---|---|
| X1 | L0 energy bucket | +0.9 | +0.5 | +0.1 | +0.3 | +3.9 | +4.7 | +79.0 | +31.3 | -23.3 |
| X2 | L1 quasi-static ECM, isothermal | -0.0 | -22.7 | -57.0 | -1.0 | -5.0 | -10.6 | -36.5 | +23.8 | -20.3 |
| X3 | L2 isothermal, no throttling | -0.0 | -22.9 | -57.0 | -1.0 | -5.1 | -10.7 | -10.8 | -11.5 | +9.1 |
| X4 | constant current I = P/3.87 V | +0.1 | +0.1 | +0.1 | +0.0 | +0.2 | +0.1 | +16.2 | +11.4 | -8.5 |
| X5 | 2RC (R·Q-scaled fit) | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +2.7 | +1.4 | -1.0 |
| X5n | 2RC, R_a + R_b = R1,ref | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | -0.0 | +0.0 | -0.0 |
| X6 | R0 split: 35 mΩ PCM + 25 mΩ cell | -0.0 | -0.3 | -0.2 | -0.0 | -0.1 | -0.1 | +52.7 | +22.9 | -18.1 |
| X7 | no burst envelope | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +24.9 | -20.8 |
| X8 | single threshold V_sys = V_cut | -1.1 | -1.7 | -1.2 | -0.2 | -1.1 | -1.3 | +0.0 | -20.9 | +17.0 |
| X9 | KiBaM instead of RC | -0.2 | -0.3 | -0.1 | -0.1 | -0.2 | -0.2 | -29.8 | -1.6 | +0.9 |
| X10 | RRC tail as ODE state | +0.0 | +0.0 | +0.0 | +0.0 | -0.4 | -0.5 | -1.3 | +0.0 | +0.0 |
| X11a | brightness law γ_b = 1 | -2.9 | +0.0 | +0.0 | +0.1 | -0.2 | -0.3 | +13.5 | -3.4 | +0.4 |
| X11b | brightness law γ_b = 2.9 | +2.2 | +0.0 | +0.0 | -0.3 | -0.2 | -0.1 | -18.9 | +2.6 | -0.3 |
| X12 | SOC output D2 (s0, s_off) | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | +0.0 | -1.4 |
| X13 | SoC leakage, k_lk = 0.04 /K | -2.6 | -0.8 | -2.4 | -0.4 | -1.3 | -1.3 | +16.2 | +9.7 | -1.1 |
| X14 | audio on DSP (a = 0 in AU) | +0.0 | +0.0 | +0.0 | +0.0 | +2.7 | +0.5 | +0.0 | +0.0 | +0.0 |
| X15 | lognormal sessions, CV 1.5 | – | – | – | – | +1.6 | +2.4 | – | – | – |
| X16 | diurnal CTMC | – | – | – | – | +8.1 | +6.0 | – | – | – |
| X17 | calendar law t^0.75 | +0.1 | – | – | +0.0 | – | – | – | – | – |
| X18a | OCV stretch 4.40 V / 3.85 V | -0.5 | -0.5 | -0.5 | -0.1 | -0.5 | -0.6 | -4.1 | -0.9 | +0.2 |
| X18b | OCV stretch 4.48 V / 3.89 V | +0.5 | +0.5 | +0.5 | +0.1 | +0.4 | +0.5 | +5.4 | +1.4 | -0.6 |

- X17 (A1, 25 °C): SOH at 1/2/3 y 0.894, 0.834, 0.784 vs X0 0.883, 0.834, 0.795; years to 0.8: 2.66 vs 2.86 (k_H, k_L recalibrated to SOH 0.8 at 1000 FCE; q_cal(1000 d, C0) 0.0958 vs 0.0782).
- Reading. ≥ 25 °C, fresh, from 100 %: the ECM-form variants (X2–X6, X9, X18) move TTE by ≤ 1.0 % (X2/X3 excluded where they remove throttling); L0 is ≤ 1.9 % long on constant loads and schedules (CTMC days +3.9 / +4.7 %); V_sys = V_cut costs up to 1.7 %; the load laws matter more (X10/X11/X13 up to 8.5 %; X14 on S16 +63 %). Throttling is indispensable for gaming (X2/X3 -23 % at 25 °C, -57 % at 35 °C). In the cold on an aged pack the form dominates: L0 +79 % (misses 44 pp stranded charge), R0 split +53 %, KiBaM -30 %, no envelope +25 %, V_sys = V_cut -21 % (S9/S10a at −20 °C). 2RC with equal total R1 (X5n) changes nothing; isothermal S10b (0 °C, aged) -53 %: self-heating matters in the cold.
- Usage laws: lognormal sessions leave the mean unchanged within MC error (+1.6 ± 1.3 %, +2.4 ± 1.6 %) but raise the SD +15% / +11%; the diurnal chain adds +8.1 ± 1.4 % / +6.0 ± 2.0 % (idle nights). X4 (constant current) is harmless from 100 % but +5.1 % from 50 %; S12 from 20 % is non-smooth (X4, X12 and X8 move the shutdown across a segment boundary).
- CTMC days: X0 mean 39.75 / 21.36 h (S13/S14, 400 paths, h = 10 s, paths from stochastic_common.UsageSim; [PROTO] MC L2 39.53 / 21.22 h with another RNG stream); variant deltas use common random numbers (paired SE in the JSON).

## Deviations from the spec and modelling choices (all documented in the JSON meta / docstrings)
- S1 estimator: Saltelli 2010 (SALib's), not Jansen first order: Jansen-S1 (kept as `S1_jansen`) had CI half-widths ≈ 0.12 per group at N = 256 and ΣS_i = 0.07 on S10a-G. Outputs centred before estimating.
- Surrogate includes q_batt in the F5 fixed point (as T(z), P(z) profiles) and the thermal event; the plain spec form (q_batt = 0) had 1.3 % error on S11c in a development check. S11c uses full L2 (above).
- Conditions chosen where §10.3 names none: S9 and S10a at −10 °C, SOH 0.8 (the §10.1 cold case); S11a/S11c at −110 dBm; S10a-G N = 512.
- X1 uses the F5 quasi-steady powers (the paper's analytic L0); X2 is literal L1 at T_amb, so no throttling and no burst envelope. X5 also run normalised (X5n). X6 keeps the z-shape and r_growth on both parts. X9 keeps R1 as series resistance (as AN6). X10's ODE implies an exponential tail timer (promotions λ(1 − r)E_pro), which dominates its effect. X13 is relative to 25 °C, so it lowers power in the cold (−0.08 W at T_s = −15 °C): this is why it lengthens the cold TTEs. X16 uses 5.5 h⁻¹ from 07:00 to 24:00. X17 recalibrates k_H, k_L; its TTE column is at the 2-year SOH.
- h = 10 s for S1, S2a, S2b, S16 and the CTMC days (§6.3); CTMC paths come from stochastic_common.UsageSim (exp law), not phonebatt CTMC.sample_paths, so X15/X16 share random numbers.
- [PROTO] checks: X0 reproduces §5.4 (e.g. S9 −20 °C/0.8 12.33 h, s_off 44.2 % vs 12.33 h, 44.2 %; S10a −20/0.8 10.23 h, 23.3 % vs 10.23 h, 23.3 %); A1 SOH 0.883, 0.834, 0.795 = [PROTO]. No [PROTO] exists for Sobol or X-variants.
