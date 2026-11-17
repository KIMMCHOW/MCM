# Group "scenarios" (AN3 + AN12): digest for the paper writer

Scripts `analysis/run_scenarios.py` (~100 s) and `analysis/run_stranded.py` (~30 s), both re-run end to end from
deleted outputs (`cd /home/user/MCM && PYTHONPATH=src:analysis python3 analysis/<script>`). Outputs: `results/scenarios.json`,
`results/stranded.json`; `paper/figures/fig_{trajectories,tte_s0,throttle,stranded_map}.pdf` (+ PNG previews);
`paper/tables/tab_baseline.tex`, `tab_stranded.tex`. Solver: LSODA (§6.2); S13/S14 via `run_batch` (RK4, h = 10 s,
400 paths, seed 7). Core: HEAD 911187e + uncommitted review fixes B1/B2. No core bug found (no core_issues file).

## Regression vs [PROTO] (§5.4, §12.2): 265 comparisons, all within §15.3 tolerances
- Deterministic: max TTE deviation +0.09 % (S7 15 °C, 4.404 vs printed 4.40 h), all others within ±0.04 %; s_off
  ±0.05 pp; T_end/Ts_end ±0.05 K (PROTO prints 1 decimal); every event identical; 14 demand powers within ±0.5 mW.
- S13/S14 MC L2: mean 39.70 vs 39.53 h (+0.43 %), 21.42 vs 21.22 h (+0.93 %); median 40.04 vs 39.38, 21.31 vs
  20.87 h; SD 7.44 vs 7.05, 4.80 vs 4.77 h; brown-out share 10.5 vs 7 %, 8.0 vs 11 %. Means within the 3 % tolerance;
  the rest is path-draw noise (core: 200/200 identical events on identical paths). One-off check, 2 x 1000 extra
  paths (seeds 11/12, not stored): S13 39.82/39.62 h, brown-outs 9.2/11.7 %; S14 21.38/21.35 h, 8.9/9.3 %.
- §5.4 "Mean P" of S7 (2.94/1.65/3.89 W) is the throttled end power P_end (ours 2.937/1.652/3.885 W); time-mean
  delivered power is 3.06/1.71/3.90 W (tab_baseline shows demand -> time mean).

## Baseline (25 °C, new, s0 = 100 %; tab_baseline.tex)
- S1/S2a/S2b 442.7/378.3/173.8 h; S3 38.88; S4 16.64; S5a/b/c 32.32/25.49/21.10; S6 11.83; S8a/b 6.93/10.83;
  CM 10.40 (brown-out); S16 79.20; S12 34.92 h. All others cut-off; s_off <= 1.1 %.
- s0 = 50/20 % runs of every constant scenario are in the JSON (e.g. S4 7.90/3.06 h).

## S7 gaming / throttling (fig_throttle, log time axis)
- 15/25/35 °C: TTE 4.40/5.68/10.29 h (brown-out/cut-off/cut-off); P_end 3.88/2.94/1.65 W; time-mean rho_perf
  1.00/0.853/0.523; T_s,end 38.2/42.4/44.7 °C. Throttle engages within ~10 min (tau_th = 7.5 min).
- rho_perf is below F5 (0.864/0.531) because F5 neglects battery heat, which raises T_s late in the discharge.
- No throttling (and no burst envelope, as in PROTO): 4.40/4.43/4.45 h; T_s,end 48.2 °C (25), 58.1 °C (35).

## S9 cold hike (T0 = 20 °C; all cut-off during a 30 s photo burst): s_off / TTE
- SOH 1.0: 25 °C 1.7 %/26.99 h; 10 °C 2.7/26.66; 0 °C 3.8/26.33; -10 °C 7.2/25.33; -20 °C 16.9/22.67.
- SOH 0.9: 2.1/24.16; 3.9/23.66; 6.4/23.00; 13.1/21.33; 30.4/17.17.
- SOH 0.8: 2.7/21.33; 5.4/20.66; 9.9/19.66; 22.2/17.00; 44.2/12.33. All six [PROTO] cells match.
- Cold-soaked start (T0 = T_amb) gives the same 17.00 h (-10 °C) and 12.33 h (-20 °C) at SOH 0.8.
- O1 has no effect on S9 (shutdown is the explicit 4.1 W burst vs V_cut, not the envelope): O7 is the remedy.

## S10a stranded-charge map and O1 governor (fig_stranded_map, tab_stranded.tex)
- Map matches [PROTO] to ±0.05 pp / ±0.04 %: 0.4 % (25 °C, SOH 1.0) to 23.3 % / 10.23 h (-20 °C, SOH 0.8); all
  brown-out except 25 °C / SOH 1.0.
- Online O1 cap dPb_eff = min(dPb, max(0, 0.9 dPb_safe)) (F4 on the state) is exactly the envelope-off trajectory
  (proof in the script docstring; checked V_pk,gov - V_sys >= 0.010 V on all 15 runs). s_off drops to 0.4-2.4 %
  everywhere: -20 °C/0.8 23.3 -> 2.4 % (TTE +25 %, 10.23 -> 12.78 h); -10 °C/0.8 8.5 -> 1.5 % (+7.0 %, 12.15 -> 13.00 h).
- Cost: the cap binds only in the last 8.6 % of the run at -10 °C/0.8 (from displayed 10.7 %) and 25 % at -20 °C/0.8
  (from 28.8 %); time-mean delivered burst amplitude 97.7 % / 92.7 % (29 % / 21 % at shutdown).
- Fixed derating at -10 °C/0.8, m_b = 1/0.75/0.5/0.25/0: s_off 8.5/4.8/2.6/1.5/1.5 %, TTE 12.15/12.61/12.88/13.00/13.00 h
  ([PROTO] points match). Key O1 message (panel c): the online cap reaches the m_b = 0 result while delivering 98 %
  of the burst amplitude on average.
- S10b (S12 day, aged): 25 °C SOH 1.0/0.9/0.8 34.92/34.05/29.19 h (cut-off); 0 °C 34.50/29.70/24.78 h (brown-out,
  s_off 1.9/2.4/4.7 %); matches [PROTO].
- Navigation at -10 °C (new, cold-soaked): 11.43 h ([PROTO] 11.43); brown-out at s_off 3.3 % (§5.4 gives no event);
  T_end -5.2 °C.

## S11 weak signal ([PROTO] -12 % / -23 % at -110 dBm reproduced)
- S11a SO-L at -85/-95/-105/-110/-120 dBm: 14.46/13.72/13.05/12.73/12.15 h (P 1.218 -> 1.447 W).
- S11c VC: 6.93/6.17/5.58/5.33/4.87 h (P 2.518 -> 3.610 W).
- S11b (new, no [PROTO]): Wi-Fi cliff at RSSI -85/-90 dBm: 15.17/9.80 h (-8.8 % / -41 % vs S4).

## S12 / S15 initial-SOC sweep (fig_tte_s0; "S15" block of scenarios.json)
- S12 from 100/80/50/20/10 %: 34.92/28.68/11.79/8.70/2.71 h ([PROTO] 100/50/20 % reproduced). TTE(50)/TTE(100) =
  0.34 with steps: a time-of-day effect (every run starts at 07:00 into heavy use), not a battery effect.
- S13 (CTMC mean): 39.70/31.16/19.36/8.01/4.29 h.
- 25 °C new: TTE(50)/TTE(100) = 0.475/0.476/0.471 (S4/S5a/S7), 0.49 (S13). Max deviation from proportionality
  -2.6 to -2.9 pp, almost exactly the D1-vs-D2 gap (energy fraction at s = 0.5 is 0.477): "nearly linear" (§10.5)
  confirmed.
- Cold/aged S4: linear with an offset. -10 °C/0.8: the 50-100 % chord crosses zero at s0 = 12 %. -20 °C/0.8: the
  chord crosses zero at 25 %, and a cold-soaked start at s0 <= 25 % browns out at t = 0. That is above the 23.3 %
  reached from 100 %, where the pack has self-heated by 3.5 K.
- S9 -20 °C/0.8: threshold-like. From s0 <= 45 % it dies at the first or second photo burst (0.16-0.5 h); above
  that, ~0.25 h per displayed percent; chord crosses zero at 44.8 %.

## Caveats
- Model outputs only: no observations enter these scenarios and nothing was tuned. S13 at each s0 reuses the
  same 400 paths (common random numbers).
