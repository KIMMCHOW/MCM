# Group "stochastic" (AN5 + AN11): digest for the paper writer

Scripts: `analysis/run_ctmc.py` (AN5, 253 s), `analysis/run_usage_fluct.py` (AN11, 220 s), shared
`analysis/stochastic_common.py`. Both were rerun end to end from deleted outputs (exit 0, numbers identical to the
first run; all seeds fixed). Outputs: `results/ctmc.json`, `results/usage_fluct.json`, `fig_ctmc_tte`,
`fig_cv_scaling`, `fig_usage_fluct` (PDF + PNG), `paper/tables/tab_ctmc.tex`. Core: HEAD 911187e + uncommitted
review fixes (meta `git_src_dirty: true`). No core bug found, so no core_issues file. L2 = `run_batch` (rk4_batch,
h = 10 s); L0 = exact event-driven energy depletion; RP-63, 25 °C, new, s0 = 100 %.

## AN5 closed forms (F6–F8): every §7 [PROTO] number reproduced to its printed rounding
- Typical / heavy: P̄ 0.4400 / 0.8224 W; σ² 0.2492 / 0.8388 W²h; μ_fluid 40.43 / 21.63 h; μ_Dynkin 41.25 /
  22.50 h; CLT SD 7.21 / 5.18 h (CV 0.175 / 0.230); τ_c 0.234 / 0.362 h; CV_P 1.659 / 1.310; screen-on 6.50 /
  11.53 h/day; IG 5/50/95 % 30.25/40.59/54.49 and 14.72/21.85/32.50 h; IG Pr(TTE < 24 h) 0.0015 / 0.650.
- Energy shares by activity and component: within 0.0005 of [PROTO] (plus wake-ups 0.001, sensors 0.003).
  Hub variance shares: game 0.807, nav 0.102, social 0.032, vcall 0.032, audio 0.021 (heavy: game 0.940).
  Deviation matrix = closed form to 8e-16 (σ²_hub 0.236495 / 0.820905); two-state formula equal to 1e-16.
- μ_fluid elasticities (F8 hub) d_game/d_social/d_audio/λ_on: −0.257/−0.169/+0.065/−0.531 (= [PROTO]). Full-chain
  finite differences: μ_fluid −0.242/−0.171/+0.063/−0.530; μ_Dynkin −0.221/−0.175/+0.061/−0.531. The SD
  elasticity to d_game is +0.354: longer games shorten the mean but widen the spread.

## AN5 closed forms vs Monte Carlo
- L0 MC (20 000 paths, seed 1), typical: 41.30 ± 0.05 h ([PROTO] 41.34), SD 7.06 (7.04), 5/50/95 % 29.52/41.41/
  52.76 (29.58/41.46/52.69), Pr(< 24 h) 0.0073 (0.009). Heavy: 22.45 ± 0.04 h (22.52), SD 5.11 (5.11),
  Pr(< 24 h) 0.626 (0.623). The differences from [PROTO] are RNG streams, all ≤ 2 SE.
- Dynkin minus MC: −0.055 h (z = −1.1) and +0.055 h (z = +1.5). The pathwise Dynkin identity residual is
  +0.025 ± 0.051 h, so the identity holds; the depletion-state approximation Σπ_x P_x h_x/P̄ costs ≤ 0.006 h.
- CLT SD is 2.2 % / 1.4 % too high (MC/CLT 0.978 / 0.986 ± 0.005): a finite-horizon effect, not noise.
- MC skewness −0.06 / +0.16 vs IG +0.54 / +0.73. IG quantiles miss by +0.7/−0.8/+1.7 h (KS D 0.047 / 0.035).
  N(μ_Dynkin, σ_CLT) fits better: errors ≤ 0.35 h, D 0.011 / 0.017.
- L2 MC (2000 paths, seed 2026), typical / heavy: mean 39.60 ± 0.16 / 21.49 ± 0.11 h, SD 7.08 / 4.93 h, 5/50/95 %
  27.80/39.73/50.92 and 13.50/21.50/29.75 h, Pr(< 24 h) 0.015 / 0.704, brown-out 9.7 / 9.2 %, s_off 0.76 / 0.97 %.
- L2/L0 on the same paths: 0.964 [0.963, 0.965] / 0.956 [0.954, 0.957]. Jensen (F10) gives 0.986 / 0.982, about
  40 % of the loss; the rest is unthrottled game starts and the end-of-discharge brown-out or cut-off.
- [PROTO] configuration (400 paths, seed 7): 39.70 h (+0.43 % vs 39.53) and 21.42 h (+0.93 % vs 21.22), within
  the 3 % tolerance; SD 7.45 vs 7.05 and 4.80 vs 4.77 h; brown-out 10.5 vs 7 % and 8.0 vs 11 %. Causes: the B1
  fix plus path noise (binomial SE 1.4 pp; the 2000-path value is 9.7 / 9.2 %).
- Rate scaling κ: Dynkin 40.63/41.25/43.71 h, CLT SD 3.61/7.21/14.43 h at κ = 0.25/1/4 (= [PROTO]). MC/CLT SD
  over κ = 0.25…4: 0.998/0.984/0.978/0.957/0.912 (typical), 0.995/1.005/0.986/0.962/0.923 (heavy): the √κ law
  holds within 2 % for κ ≤ 1 and over-predicts by 9 % at κ = 4 (fitted slope 0.47).
  - L2/L0 means: 0.920/0.940/0.961/0.978/0.992. Fast switching costs energy (short games never throttle), and
    brown-outs reach 23 % at κ = 0.25 (4 % at κ = 4). Step check κ = 0.25, h 10 vs 5 s: +0.008 % in the mean.

## AN11 usage fluctuations (§10.5)
- Burstiness, 1.2 W constant vs 0.4 ↔ 2.0 W (L0 unchanged):
  - 25 °C, new: −0.47 % up to 5 min. At 30–60 min the phase mean is −0.1 to −0.8 % (energy-equivalent −0.5 %);
    a single run spreads +0.6/−2.2 % by phase. L1 brackets −0.47/−0.53 % and F10 −0.44 % agree.
  - −10 °C, SOH 0.8: −2.6 % (10 s), −3.0 % (10 min), −1.4 % (1 h); L1 −2.65 % (short), −3.41 % (long),
    −2.63 % (long, phase temperatures); F10 −1.89 %.
  - Adding 3 W bursts (no effect at 25 °C): the constant load browns out at 10.63 h (s_off 9.2 %). The
    alternation costs −7.3 % at 10 s (L1 −7.48 %), −8.4 % at 10 min, −3.6 % at 1 h (thermal L1 −5.6 % plus
    phase snapping), and s_off rises to 15.6 %. F10 misses this: cold burstiness acts through the brown-out.
- Order effect, S12 reversed. 25 °C: +5.4 % in L2 and +5.5 % in L0, so it is pure scheduling.
  - −10 °C, SOH 0.8: forward 11.96 h (brown-out at 32 % in the 18:00 game), reversed 23.48 h (10.7 %).
    That is +96 % in L2 vs +26 % in L0, so +70 pp comes from the voltage limit.
  - Start-time scan: 25 °C 26.3–43.3 h; cold aged 12.0–28.5 h, s_off up to 39 %.
- Initial SOC (S4 with 3 W bursts):
  - TTE(50 %)/TTE(100 %): 0.475/0.470/0.463/0.444 at 25/0/−10/−20 °C new; 0.431 (−10 °C) and 0.330 (−20 °C)
    at SOH 0.8. At −20 °C, SOH 0.8 the phone dies at once for s0 ≤ 25 %.
  - TTE is linear for s0 ≥ 0.2 within 1.3–1.5 % (4.7 % at −20 °C, SOH 0.8). Upper-half fit intercepts
    5.2/6.1/7.3/10.5/12.5/25.8 % vs the F4/F2 empty point 0.4/1.4/2.8/6.5/8.7/23.6 % (≈ +4 pp OCV-shape offset).
- Heavy-tailed sessions, L0 MC SD vs F9 (typical): CV 0.5 5.66 vs 5.72; CV 1 lognormal 6.89 vs 7.21; CV 1.5
  8.22 vs 9.17; CV 3 10.73 vs 16.06 h.
  - F9 is exact asymptotically (hub at 50 E_u: MC/F9 0.94–1.01), but at the real E_u it is 0.87 (CV 1.5) and
    0.63–0.66 (CV 3): a battery cannot outlast a long session. L2 at CV 1.5 on coupled paths: SD ×1.10 / ×1.14.
- Per-visit jitter ±30 % on b, APL, λ, d_dl: state-power CV ≤ 10 %, σ² +0.2 %, SD ratio 0.998–0.999 (closed
  form, L0 and L2), mean shift ≤ 0.10 h.
- Usage-model parameters (4096 Sobol' samples, F6–F7):
  - μ_Dynkin 5/50/95 % 32.0/44.1/58.5 h (typical), 18.0/23.8/31.5 h (heavy); share of the TTE variance 54 % /
    37 %; mixture 5/50/95 % 27.8/43.6/64.1 h vs nominal IG 30.3/40.6/54.5 h.
  - First-order indices p_game 0.45, λ_on 0.28 (typical; heavy 0.54, 0.24, d_game 0.12). OAT: λ_on ×0.7/1.3 →
    +22.2/−12.4 %; p_game 5/95 % → +20.4/−20.1 %; F8 first order within 1–2 pp.
  - Caveat: LU(0.7, 1.3) has a geometric mean of 0.954, not the listed 1.0, so the sample median is above the
    nominal 41.25 h (a property of the §4.6 ranges, not of the code).

## Take-aways
1. F6 is exact to MC noise and F7 good to about 2 %; √κ holds for κ ≤ 1; the TTE law is closer to normal than IG.
2. Path randomness gives a CV of 0.18 / 0.23. The usage-model ranges double the variance; game frequency and the
   session rate dominate. Per-visit jitter does not matter; heavy-tailed sessions add 10–16 % to the SD.
3. At 25 °C only the energy matters. In the cold on an aged pack, when the heavy use happens decides the day:
   order +96 %, bursts −7 %, stranded charge up to 39 %.
