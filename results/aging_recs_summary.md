# Group "aging_recs" (AN13 + AN14 + AN15): digest for the paper writer

`analysis/run_aging.py` (~300 s, 2 workers), `run_recommendations.py` (~40 s), `run_generalization.py` (~2 s), re-run end to end
in that order from deleted outputs. Outputs: `results/{aging,recommendations,generalization}.json`, `fig_{aging,levers,
generalization}.pdf` (+ PNG), `tab_aging.tex`. Core: HEAD 911187e plus the uncommitted review fixes. No core bug found.

## AN13 aging (§11): every [PROTO] row reproduced
- **C0 calibration:** k_c 1.4132e-3 d^-1/2, k_H 2.6970e-3, k_L 1.2883e-3 FCE^-1/2, q_cal(1000 d) 0.0782, mean g 1.749,
  SOH(1000 FCE) 0.8000 (= [PROTO]). With N_80 = 500 / 800, A1 reaches SOH 0.8 at 1.55 / 2.35 yr (2.86 at 1000).
- **A1–A7:** SOH at 1/2/3 yr within ±0.0005 and years to 0.8 within ±0.005 of §11.3 (tab_aging.tex).
  2-yr change vs A1 (pp): A2 +1.89, A3 +0.76, A4 −12.44, A5 −5.98, A6a −3.73, A6b −5.04, A6c +1.91, A7 −1.40;
  A6c − A6b +6.96. Loss split at 2 yr, A1: calendar 0.080 / SEI 0.074 / plating 0.013; A6b plating 0.094.
- **Joint UQ** (N = 128 Sobol': k_c, E_cal, g_hi, κ_HV, w_L, E_H, E_L, β_L, γ_R, ΔT_ch; each sample recalibrated):
  - SOH at 2 yr: A1 0.816–0.851, A4 0.683–0.791. Robust in sign (5–95 %): A2 − A1 +0.4…+2.9 pp and A3 − A1 +0.05…+1.3
    (positive in 100 % of samples); A4, A5, A6a, A7 negative in 100 %; A6c − A6b +1.2…+11.4. Not robust: A6c − A1 −6.1…+3.7 (58 %).
  - P(SOH 0.8 within 3 yr): A1 62 %, A2 16 %, A3 36 %, A4 100 %, A5 97 %.
- **β_L ends** (0.5 vs 2.64), SOH at 2 yr: A6a 0.830 vs 0.796, A6b 0.780 vs 0.783, A6c 0.799 vs 0.853, A1 0.833 vs 0.834.
  The fast-charge penalty and the JEITA advantage over A1 rest on the upper-end (LFP) β_L; the 0 °C-charging penalty does not.
- **TTE vs age:** S13 Dynkin (capacity only) 41.25 → 34.52 h after 2 yr of A1 (A4 29.49 h); S12 L2 25 °C 34.92 → 29.72 h.
- **SOH grid (fig b):** S4 TTE relative to new at SOH 0.8 is 0.794 / 0.779 / 0.751 / 0.658 at 25 / 0 / −10 / −20 °C
  (0.589 / 0.555 / 0.488 / 0.340 at SOH 0.6).
  Channel split at −10 °C, SOH 0.8: capacity alone 0.800, resistance alone 0.939 (0.993 at 25 °C). The stranded charge is all
  resistance channel: s_off 2.76 % with capacity loss only vs 8.54 % with both.
  F11 closed form vs L2 (25 °C, SOH 0.8): 0.7965 vs 0.7944 (max |diff| 0.0054).
- **Winter shutdown (fig c; the s_off ≥ 5 % displayed threshold, S4 cold-soaked, is our illustrative choice):**
  - Reached at SOH 0.898 (−10 °C) / 0.733 (0 °C); a new pack already strands 6.4 % at −20 °C. At −10 °C it is reached after
    A1 0.77 yr, A2 0.98, A4 0.26, A5 0.42.
  - s_off at −10 °C after 1/2/3 yr: A1 5.4 / 7.1 / 8.8 %; A4 8.7 / 13.4 / 17.9 %.

## AN14 recommendations (§12)
- **Levers:** all 13 §12.1 values reproduced within 0.05 pp (vectorised engine = core `lever_effects`, 0.0 pp difference).
  Base: S12 P̄ 0.4954 W (E_u/P̄ 35.91 h); S13 Dynkin mean 41.25 h.
  S13 ranking: skip gaming +36.9 %, 60 Hz +4.7, dark theme +4.6, voice call +4.0, prefer Wi-Fi +3.7, brightness −0.1 +3.1,
  OS frame cap +2.7, video call on Wi-Fi +1.6, OS coalescing +0.7, halve background +0.4. Avoid: 5G standby −4.5,
  weak signal −5.7, brightness +0.2 −7.5.
- **90 % parameter intervals** (N = 1024, §4 distributions):
  - Skip gaming S13 +19.7…+36.0 (median +26.8). The nominal value sits at the top because the R_skin and k_env priors are
    skewed towards more throttling (median steady game power 2.59 W vs 3.12 W nominal). Same effect for 60 Hz (median +3.4)
    and the OS frame cap (+1.3); fig_levers marks the medians.
  - Skip gaming ranks first in 100 % of samples. Every gain is positive in ≥ 98 % of samples; every "avoid" item is negative in 100 %.
  - With the §4.6 usage factors added, skip gaming spans +5.6…+71 %.
- **Non-smooth L2 S12 TTE** (not used for ranking): skip game +37.7 %, dark +1.4, 60 Hz +0.67, weak signal −0.79 (= [PROTO]).
- **S12 energy shares at the F5 fixed point:** SoC 40.0 %, display 27.2, radio 14.5, awake 7.2, suspend 4.0, GNSS 3.4, audio 1.7,
  camera 1.3, sensors 0.6. These differ from [PROTO] (see deviation 1).
- **Scenario levers (L2):**
  - Video on Wi-Fi vs 5G 32.32 vs 21.10 h (+53 %, 90 % +34…+102); video call Wi-Fi vs LTE 10.83 vs 6.93 h (+56 %); voice +77 %.
  - S9 kept warm: s_off 7.2 / 16.9 % (new, −10 / −20 °C) and 22.2 / 44.2 % (SOH 0.8) vs 1.7 / 2.7 % at 25 °C.
  - Hot car: 10.29 vs 5.68 h at ρ_perf 0.52. Weak signal: social −12.0 %, video call −23.2 %.
- **OS strategies:**
  - O1 (−10 °C, SOH 0.8, m_b 1 / 0.5 / 0): s_off 8.5 / 2.6 / 1.5 %, TTE 12.15 / 12.88 / 13.00 h. At −20 °C: 23.3 → 2.4 %,
    10.23 → 12.78 h. The AN12 online governor reaches the m_b = 0 result while delivering 97.7 % of the burst amplitude.
  - O2: game power 3.116 → 2.722 W; S12 +3.4 % (0.4–5.2), S13 +2.7 %. L2 S7 gains +13.7 %, but mean performance falls to
    0.78 of requested vs 0.85 under reactive throttling, so the gain is not free.
  - O3: standby +2.4 / +25.4 / +66.7 % (Wi-Fi / LTE / 5G), daily +0.1 / +0.7 %. O4: avoided −4.5 % (5G standby), −34.7 %
    (5G streaming); all cellular on 5G −8.6 %. O5: −5.7 % (S13), −12.0 % (social). O6: +4.1 / +4.7 %.
  - O7: F4 at the measured battery temperature predicts the L2 stranded charge within 0.2 pp (−20 °C, SOH 0.8: 23.1 vs 23.3 %).
  - O8: IG band from full charge 30–54 h (median 41). Coverage of the nominal 90 %: 91.1 % (L0), 91.9 % (half charge), 88.3 % (400 L2 paths).
  - O9: 80 % limit +1.89 pp, optimised charging +0.76, hot climate −12.44, car −5.98, JEITA +6.96.

## AN15 generalisation (§13). Device classes are illustrative orders of magnitude [unverified]
- **RP-63:** Π1 0.0047 at 1 W (0.036 for S4 cold-aged); Π2 1.19 gaming (0.32 social); Π3 0.022; Π4 0.110 → 1.18 at −20 °C,
  SOH 0.8, z = 0.2 (1.35 with RC polarisation); Π5 0.175 (0.230 heavy).
- **Size invariance:** per cell-equivalent, Π1 = c·RQ_Σ/Ū and Π4 depend only on C-rates and R·Q, not on pack size or N_s.
- **Classes (Π1 / Π2 / Π4 / Π4 at −20 °C, SOH 0.8 / Π5):**
  - Tablet 0.0085 / 0.89 / 0.09 / 1.4 / 0.24. Laptop 0.0075 / 1.8 / 0.15 / 2.0 / 0.18.
  - Watch 0.0078 / 1.2 / 0.31 / 4.3 / 0.16. Earbuds 0.011 / 0.13 / 0.11 / 2.1 / 0.04.
  - E-reader 0.0017 / 0.17 / 0.13 / 1.1 / 0.08. Drone 0.032 / 0.30 / 0.32 / ∞ (hover alone breaches V_sys) / 0.03.
  - Input-range corners move each Π by about one decade (JSON). E-reader example: about 48 days on 5 Wh.
- **Map:** drones have the largest Π1 and Π3 (0.32), watches the largest Π4; cold + aged + low charge pushes every class to
  Π4 ≳ 1. Laptops and watches throttle (Π2 > 1); earbuds and e-readers are the most predictable.

## Deviations from [PROTO] and caveats
1. **S12 energy shares:** the §12.1 [PROTO] shares (display 32.1, SoC 29.4 %) are reproduced exactly by final_proto's formula.
   That formula takes the game's T_s from its unthrottled demand, giving a SoC power of 0.58 W instead of 2.37 W and
   P̄ = 0.421 W, which contradicts its own 0.4954 W lever base. Use the fixed-point shares.
2. **L2 S12 "5G standby":** −0.88 % vs [PROTO] −1.5 %. final_proto 'recs' used SB-5's own λ (1/60 s⁻¹); the package lever (and
   the matching smooth metric) keeps each state's λ.
3. **Π4:** the [PROTO] values 0.11 / 1.2 omit RC polarisation (§7 table convention); including it gives 1.35.
   Π1 0.0047 vs printed 0.005, Π3 0.022 vs 0.02 and A3 − A1 +0.76 vs 0.7 pp are rounding differences.
4. **A6b/A6c** charge at 0 °C + ΔT_ch but hold overnight at 25 °C (as [PROTO]). Usage-factor intervals: `ctmc_ci90_usage_pct`.
