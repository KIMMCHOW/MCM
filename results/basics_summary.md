# basics: AN0, AN1, AN2, AN4, AN16 (digest for the paper writer)

Scripts in `analysis/` (run with `PYTHONPATH=src:analysis`); all deterministic (no seed); total ≈ 15 s. Model code:
HEAD `911187e` plus the uncommitted review fixes in `src/` (`meta.git_src_dirty = true`). [PROTO] comparisons
below are regression checks, not validation.

## AN0 parameters (`params.json`, `paper/tables/tab_params.tex`)
* 93 `Params` fields with value, unit, distribution, source keys and confidence. An independent transcription of
  the §4 tables (84 values and UQ ranges) agrees with `phonebatt.params` exactly: **0 mismatches**.
* T1: 38 rows in 5 groups (symbol, meaning, value, unit, UQ, source). Every `\cite` key resolves in `refs.bib`;
  fits the 467 pt text width (checked with the paper preamble). Caption should explain: U uniform, LU
  log-uniform, T triangular (mode = value); fit = our fit to open data/PyBaMM, der. = derived, assum. =
  assumption, **cal.** = the single [CAL] constant ΔP_awake (V1.1).
* Derived: z_e 0.0263, U(1) 4.4496 V, mean OCV 3.8699 V, E_n 18.189 Wh, E_u 17.788 Wh (all = [PROTO]).

## AN1 primitive powers (`power.json`, `fig_power_breakdown.pdf`, `tab_scenarios.tex`)
* All 19 primitive powers match §5.1 (max 0.47 mW off the 3-decimal spec column = rounding; tolerance 0.5 mW);
  every component equals `final_proto.py` exactly.
* Game at 25 °C: demand 3.906 W, F5 cap P* 3.116 W, ρ_perf 0.864, non-SoC 0.743 W; GM-L capped 3.21 W.
* Fig (a): stacked components before throttling, F5 cap tick, published max-load band 2.41–5.84 W (NBC22–24).
  SoC dominates game/burst, radio dominates VC (1.2 W real-time LTE), display dominates RD/SO/NV, floor dominates
  audio. SB-W 0.040 W = 0.226 %/h, inside the 0.15–0.5 %/h standby band (BG9).
* Fig (b), V6: model components vs the power implied by the published runtime (F2 inverted on the device's own
  battery, 25 °C; the rated-energy value of `device_benchmarks.json` is stored too). **9 of 14 model rows pass
  ±25 %** (the spec's "11 of 16" adds the GM-band and standby rows); V6 TTEs = [PROTO] within 0.03 %. Failures:
  Android idle (implied 0.38 W S24, 0.69 W Pixel 8 vs model 0.245 W), Android local video (0.74–0.76 vs 0.48 W),
  iPhone audio (implied 0.167 vs 0.224 W, −25.5 %, marginal). This is the vendor floor spread; nothing re-tuned.
* T2: 19 primitives (inputs, ΔP_b, P with the capped value in parentheses for GM/GM-L, scenarios that use each;
  C = usage CTMC). Footnotesize; fits the width.

## AN2 cell maps (`cell_maps.json`, `fig_cell_maps.pdf`)
* OCV hardly depends on T (entropic shift ≤ 30 mV over −20…45 °C); cold acts through resistance: R0(0.5) =
  218 / 117 / 60 / 42 mΩ at −20 / 0 / 25 / 40 °C, ×1.76 higher at z = 0 than at z = 0.5.
* (a) OCV and the loaded quasi-static voltage at 1.06 W: the loaded curves fan out with T, the OCV does not.
  (b) R0(z,T), dashed = −20 °C, SOH 0.8. (c) P_cut,sys(z, T_batt) map (log z) with the brown-out boundaries
  P + ΔP_b = 4.1 W (S4), 4.5 W (S6), 7.9 W (S7 unthrottled) and the L2 shutdown points (open = brown-out).
* All 9 L2 brown-outs (S4, S6 at 0/−10/−20 °C; S7 at 15/0/−10 °C) lie on the F4 boundary: |z_end − z*| ≤ 0.0012.
  S7 browns out at s_off 1.9 / 4.3 / 9.4 % (15 / 0 / −10 °C); S6 at −20 °C at 8.0 %. The runs reproduce
  [PROTO] within 0.09 % TTE and 0.05 pp s_off.

## AN4 analytic checks (`analytic.json`, optional `tab_pcut.tex`)
* **T3, F2 vs isothermal L2 (25 °C, 0.5–3.5 W): max |Δ| = 0.0009 %** ([PROTO] ≤ 0.001 %, tolerance 0.05 %).
  Extra rows: 0 °C ≤ 0.006 %; −20 °C ≤ 0.053 % at 3.5 W (RC with τ1 = 135 s not yet equilibrated at the start;
  T3 itself is defined at 25 °C).
* **F3: a = 0.00495 (ohmic) + 0.00322 (cut-off) = 0.00818 W⁻¹** (−0.03 % vs [PROTO]); U′(z_e) = 8.820 V. Exact
  L1/L0 losses −0.41 / −0.89 / −3.15 % at 0.5 / 1.06 / 3.5 W ([PROTO] −0.4 / −0.9 / −3.1); first order −aP gives
  −0.41 / −0.87 / −2.86 %. a grows ×2.1 (0 °C), ×2.9 (−10 °C), ×4.2 (−20 °C), ×1.5 (SOH 0.8) and **×6.3 (−20 °C
  with SOH 0.8)**: the spec's "×3–5" holds for one stressor, not for the two combined.
* **F4: z* = 0.0881, s_off 6.35 % vs L2 6.38 %** (T_batt −16.6 °C), all = [PROTO]. On the 15-cell S10a map F4 is
  within 0.17 pp of L2 in all 14 brown-out cells (largest: 23.09 vs 23.26 % at −20 °C, SOH 0.8); the L2 map
  reproduces [PROTO] within 0.045 pp and 0.037 % TTE. Safe burst on the social load at z = 0.5 / 0.2 / 0.1:
  27.1 / 17.4 / 12.7 W at 25 °C new, but 3.9 / 2.2 / 1.4 W at −20 °C, SOH 0.8.
* **F5 thermal cap (GM)** (all = [PROTO]; sharp-cap bound 3.27 / 1.45 W at 25 / 35 °C):

  | T_amb | P* (W) | ρ_perf | E_u/P* (h) | L2 TTE (h), event | L2 ρ_perf |
  |---|---|---|---|---|---|
  | 15 °C | 3.902 | 0.999 | 4.56 | 4.40, brown-out | 0.999 |
  | 25 °C | 3.116 | 0.864 | 5.71 | 5.68, cut-off | 0.853 |
  | 35 °C | 1.728 | 0.531 | 10.30 | 10.29, cut-off | 0.523 |

  dTTE/dT_amb (L0 formula) is 0.33 h/K at 25 °C and 1.08 h/K at 35 °C; the L2 secant over 25→35 °C is 0.46 h/K.
  Caveat for the text: E_u/P* ≈ L2 at 25 °C comes from two errors that cancel. P* neglects battery heat (the L2
  mean power is 3.06 W, not 3.12) and L0 neglects resistive loss; F2 at P* gives 5.58 h (−1.8 %).
* **P_cut table**: every entry within 0.05 W of [PROTO] (rounding). 25 °C new: 20.2 / 12.5 / 9.0 W (cut-off),
  28.3 / 18.6 / 13.9 W (burst); −20 °C, SOH 0.8: 3.2 / 2.1 / 1.6 and 5.3 / 3.6 / 2.7 W.
* F2 at the thermal steady state vs L2 at 25 °C: within 0.04 % for every non-throttled constant scenario (S3–S6,
  S8a/b, S16, CM); S7 −1.8 % (the q_batt effect above); L0 overestimates TTE by 0.2–1.9 %. This supports the
  "energy problem at 25 °C" reading and the F2 Sobol surrogate.

## AN16 architecture (`fig_architecture.pdf`)
usage u(t) → component load model → constant-power ECM (z, v1; OCV, R0, R1‖C1, E2 root, burst peak) ↔ thermal node;
thermal → skin T_s and throttle → SoC power; slow aging states → Q_max↓, R↑ (dashed: I, T, z drive aging over
days); outputs: SOC (D1/D2), TTE = first of five events, s_off, T, T_s, ρ_perf, P_k. Orange = electro-thermal
feedback, violet = aging coupling, ink = signal flow. Matches the existing caption in `sec_intro.tex`.

## Deviations from [PROTO] and issues
* No [PROTO] target missed. Largest deviations: primitive P 0.47 mW (rounding), S7 at 15 °C TTE +0.09 %, F4 z*
  +0.02 %, a −0.03 %, P_cut 0.05 W (rounding).
* No core bug found (hence no `core_issues_basics.md`). Notes: (i) the qualitative "a × 3–5" of §7 F3 is exceeded
  when cold and aging combine (×6.3); (ii) F5's E_u/P* agreement with L2 is partly a coincidence (above).
