# Group "uq_sens" (AN7 + AN8): digest for the paper writer

Scripts `analysis/run_uq.py` (AN7) and `analysis/run_sensitivity.py` (AN8). Outputs: `results/uq.json`,
`results/sensitivity.json`, `fig_uq_violins`, `fig_tornado`, `tab_tte_uq.tex` (T3), `tab_elasticities.tex`.
Both scripts re-ran from deleted outputs and cache with identical numbers. Cost: UQ ≈ 63 CPU-min (clean run 22 min on
3 workers, as `--only 'det:S(1|2a|2b)\|' --no-aggregate` then a plain run; cache `results/cache/uq/`); sensitivity
50 s. No core bug; see `core_issues_uq_sens.md` (hazard H1).

## Set-up and deviations from §9 (also in `uq.json` → `deviations`)
- θ: one joint scrambled-Sobol' sample of all 70 `UQ_SPEC` factors (N = 2048, seed 2026), run with L2 `run_batch`.
  The default-θ "nominal" (LSODA) matches 29 §5.4 [PROTO] TTEs within 0.08 %.
- Constant loads use h = 10 s, not 2 s. Against h = 2 s (256 θ, 5 cells): ≤ 2e-6 relative (≤ 0.03 s), except one
  θ that browns out after 34 s (0.3 s). Labels agree on ≥ 98.8 %. S9/S12/S10b: h = 2 s; standby h = 10 s.
- Cold standby N = 512. s0 = 50/20 % only for S4, S5a, S7, S12. S11a/c only at −110 dBm, 25 °C. S9 at
  25/0/−10 °C (T0 = 20 °C), plus −10/−20 °C at SOH 0.8. t* = 16 h for S12/S10b from 100 %, else 0.9 × nominal.

## AN7: deterministic scenarios (T3 = median [5 %, 95 %], h)
- **25 °C, from 100 %:** S4 14.6 [11.7, 18.8] (nominal 16.64); S7 6.6 [4.7, 9.0] (5.68); S12 33.9 [25.2, 35.4]
  (34.92; Pr(< 16 h) 0.002); S5c 15.0 [12.3, 19.0] (21.10); S16 59 [41, 90] (79.2); S1 373 [271, 584] (442.7).
  CV 0.10 (S12) to 0.33 (S11c).
- **The default θ is optimistic: median/nominal 0.62–0.97** (S7 1.16).
  - This is energy, not voltage: L0 shows the same ratios. Several defaults sit at the low-power end of their §4
    range: NR P_tail 0.255 in U(0.249, 0.593); T_tail 10.3 in U(10.2, 21.4); P_d0 0.12 in U(0.08, 0.25);
    ΔP_awake 0.09 in LU(0.05, 0.25) (median 0.112); Wi-Fi P_act 0.13 in U(0.10, 0.30).
  - Nominal inside the 90 % band: 16/19 cells (25 °C), 14/17 (0 °C), 13/18 (−10 °C). It is above q95 for S2b,
    S5a, S5c (and S5b, S8b at −10 °C). S7 runs longer because low-T_lim / high-R_skin draws throttle more.
- **Cold, new, 100 %:** S4 13.8 [9.8, 18.1] (0 °C), 13.0 [6.7, 17.6] (−10 °C); brown-out ends 89–95 % of θ.
  s_off median 4.2/9.5 %, q95 30/54 % (nominal 1.4/2.8 %). S16 and standby are unchanged.
- **Aged + cold:** S10a (−10 °C, SOH 0.8) 7.8 [1.3, 12.8] h (nominal 12.15); s_off median 34 %, q95 89 %;
  3.2 % of θ dead at the first burst. S9 at −20 °C, SOH 0.8: 4.8 [0.3, 14.9] h, s_off median 76 % (nominal 44 %).
  Pr(day dies before 23:00), SOH 0.8 (S10b) at 25/0/−10 °C: 0.31/0.84/0.92; new pack (S12): 0.002/0.18/0.37.
- **Low start SOC, cold:** TTE = 0 (brown-out at the first burst) for 34 % of θ (S4, −10 °C, 20 %), S12 27 %,
  S7 0/−10 °C 57/77 %, S7 25 °C 12 %; that is what a "0.0" lower bound in T3 means.
- Cold-standby s_off ≥ −0.3 % (P1, spec-consistent). S11c: 1 % of θ hit the 68 °C thermal shutdown.

## Delta method (§9.2: CV² ≈ CV_Eu² + Σ s_i² CV_Pi²)
- 25 °C: CV_L2/CV_delta = 0.87–1.32 over cells (median 1.09); exact L0 sample 0.92–1.13. The delta method is low
  for standby (1/P convexity; S1 0.217 vs 0.254) and high for throttled gaming (0.227 vs 0.198).
- Cold: 1.3–2.1 (−10 °C, new), 3.4–6.6 (aged + cold). The voltage limit carries the uncertainty; L0 fails there.

## S13/S14 variance split (§9.3; F6–F7 per θ, N_θ = 2048, 25 °C)
- **Statement (hardware θ only):** typical day, usage randomness 71 % of Var(TTE), hardware 29 % (SD 6.97 h; L0
  mixture 26.9/37.8/49.8 h). Heavy day 77/23 % (SD 4.96 h; 14.3/22.3/30.6 h).
- Full L2, 50 θ × 200 paths (MC-noise corrected): usage 66 % [55, 78] (typical), 74 % [65, 83] (heavy). F6–F7
  on the same θ gives 71/76 %, so the L2 hardware share is ≈ 5 pp higher.
- Adding the §4.6 usage-model factors (pick-freeze, Jansen totals; path / hardware / usage model / interaction):
  typical 38/19/44/−1 % (CI of interaction [−4, 3]); heavy 52/16/30/2 %. Mixtures 25.6/39.7/57.8 and
  14.1/23.2/34.7 h.
- **Pre-registered L2 check: FAILS.**
  - Mean within 3 % of r̄·μ_Dynkin (r̄ = pooled L2/L0 = 0.927/0.914): 48/44 % of θ (max 19.6/17.8 %); SD
    within 10 %: 64/70 % of θ (MC SE 1.2 % mean, 5 % SD).
  - Cause: L2/L0 depends on θ (0.745–0.977). Brown-out ends 58 % of paths across θ (9.7 % at default θ);
    high R0·A_R·m_b and small Δ_sys give the lowest ratios.
  - The ratio correlates 0.96/0.93 with the S4 L2/L0 ratio of the same θ (using it: max error 9.9/8.0 %). L0 MC
    matches F6 (≤ 2.4/3.5 %): the misses are the voltage layer. L2-equivalent (× r̄): 24.9/35.0/46.2, 13.1/20.4/28.0 h.

## ABC (§9.4; 14 Pro black-screen idle within ±5 % of 55.03 h)
- Prior predictive 28.2/40.9/65.7 h (default θ 54.48 h). The Bonn device range 54.4–57.3 h is inside the 90 %
  prior. 155 of 2048 accepted (7.6 %).
- Posterior ΔP_awake 0.053/0.082/0.122 W (prior 0.054/0.112/0.231; width ×0.39); P_d0 0.086/0.126/0.161 W
  (×0.50); others ×0.98–1.04. 25 °C intervals: median ×0.95 (S12 ×0.64); medians +1 to +13 %. S13 mixture
  28.5/40.0/51.7 h, hardware share 23 %.

## AN8: elasticities (§10.1; LSODA, ±5 %; voltages in % per +10 mV, T_lim in % per +1 K)
- All 36 §10.1 [PROTO] entries reproduced within 0.0005.
- Extra: S4 25 °C v0 −0.464, κ_RR −0.160, δU +0.38, V_cut −0.12. S4 −10 °C aged δU +1.0, V_cut −0.66, Δ_sys
  +0.66 (brown-out-limited: only V_sys matters), γ_R −0.097. S7 T_lim −4.46 %/K, R_skin +0.767.
- Smooth metrics: S12 E_u/P̄ (35.91 h) v0 −0.252, P_SoC,max −0.247, R_skin +0.198, T_lim −1.16 %/K; S13 (41.25 h) Q_n +0.98.

## AN8: OAT tornado (§10.2; ΔTTE % at the low / high end)
- **S4 25 °C (16.64 h):** RAT 5G −23.5 (LTE −13.1); P_SoC,max +19.4/−14.0; 60 Hz +19.0; brightness ∓0.2
  +9.8/−13.0; ΔP_awake +4.0/−13.2; dark theme +9.8. Resistances < 1 %.
- **S10a −10 °C, SOH 0.8 (12.15 h):** R0,ref +7.1/−52.4; m_b +6.0/−32.7; γ_R +3.3/−28.2; 5G −25.0; E_a,R0
  +5.1/−22.6; then the load factors; 35 °C +9.2 (cut-off instead of brown-out).
- **S7 (5.68 h):** ambient −10/35 °C −29.4/+81.1; R_skin −18.2/+50.2; k_env −14.2/+38.9; T_lim +22.3/−11.3.
- **S12 / S13 (L0 metrics):** RAT +3.9/−17.3, +6.4/−17.0; ambient −6.2/+13.5, −7.0/+13.1 (via throttling only).
- **Reading:** at 25 °C, user choices (radio, refresh, brightness, theme) span as much as the largest hardware
  uncertainties. Cold and aged, resistance and bursts dominate.
- **Captions:** tornado blue/red = low/high end, italic = behaviour/environment, upright = hardware (§4 ends).
  Violins: (a) log TTE at 25 °C, black bar = default θ; (b) TTE / 25 °C default, darker = colder; (c) shares.
