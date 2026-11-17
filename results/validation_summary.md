# Validation (AN6) summary

`analysis/run_validation.py` -> `results/validation.json` (FINAL_SPEC §8, pre-registered). Only the declared [CAL] steps were fitted (V1.1 ΔP_awake; V2 HG2 OCV/thermal at 25 °C; V4 Oxford OCV). [PROTO] numbers are regression targets, not validation. **Main model (L2): 31 of 59 criteria pass**; structural variants (X9, L1, X5): 22 of 34.

## V1 Bonn iPhones (ΔP_awake 0.090 W; exact CAL 0.0877 W, [PROTO] 0.0877)
| case | model h | obs h (range) | ε | ε blind 0.10 W | 90 % UQ, CAL fixed | (a) | (b) |
|---|---|---|---|---|---|---|---|
| 14 Pro idle | 54.48 | 55.03 (54.4–57.3) | -1.0 % | -5.1 % | 35.3–63.0 h | -- | pass |
| 14 Pro camera | 8.96 | 8.43 (8.3–8.6) | +6.3 % | +5.5 % | 5.4–10.0 h | pass | pass |
| 14 Pro GPS | 24.98 | 48.49 (45.1–51.3) | -48.5 % | -49.5 % | 17.7–29.3 h | fail | fail |
| 14 Pro camera+GPS | 7.48 | 8.38 (8.1–8.5) | -10.8 % | -11.3 % | 4.8–8.2 h | pass | fail |
| 13 Pro idle | 51.85 | 50.82 (50.6–51.4) | +2.0 % | -2.2 % | 33.6–59.9 h | pass | pass |
| 13 Pro camera | 8.52 | 9.53 (9.2–9.7) | -10.5 % | -11.2 % | 5.1–9.5 h | pass | fail |
| 13 Pro GPS | 23.77 | 29.76 (28.9–31.2) | -20.1 % | -21.7 % | 16.8–27.9 h | fail | fail |
| 13 Pro camera+GPS | 7.11 | 8.33 (8.0–8.4) | -14.7 % | -15.2 % | 4.5–7.8 h | pass | fail |
* Regression: model TTEs within 0.06 % of [PROTO]; additivity (c) model +0.41 / +0.43 % ([PROTO] +0.41), data +0.15 / -0.10 % (14 / 13 Pro, camera 60 fps as in the combined run).
* Reading: idle and camera pass (a) on both phones; GPS fails on both (pre-registered for the 14 Pro: the data imply ≈ +0.03 W, the KW19 prior is 0.27 W; the 13 Pro also fails at −20 %). Camera+GPS passes (a) on both phones, narrowly on the 13 Pro (−14.7 %); on the 14 Pro the GPS excess is partly offset by the camera's low power.
* (b) is strict: the §4 ranges of P_d0 and ΔP_cam lie mostly above their defaults, so the CAL-fixed 90 % band sits below the point prediction. It misses the 13 Pro camera range and both camera+GPS ranges, which lie above the point prediction, besides GPS. The ABC posterior (330/4096 prior samples within ±5 % of the 14 Pro idle) covers 4/8.
* V1.4 SOC-shape test, 28 traces: RMS window-slope error D1 4.4 %, D2 2.0 %; D2 better on 27/28. The data favour the energy-fraction SOC (D2), as expected. 14 Pro idle observed 31.1/30.8/32.4 min per % (spec 31.1/30.8/32.4). D1 window OCVs match [PROTO] to 0.3 mV.

## V2 LG HG2 (CAL at 25 °C only; f_Q = Q(T)/Q(25 °C); criteria |Δf_Q| ≤ 0.05 (T ≥ 0), 0.10 (−20))
| T, rate | f_Q obs | L2 | X9 | RMSE_V L2 / X9 (mV) |
|---|---|---|---|---|
| 40 °C 0.5C | 0.972 | 1.001 | 1.007 | 34 / 27 |
| 10 °C 0.5C | 0.932 | 0.998 | 0.985 | 74 / 47 |
| 0 °C 0.5C | 0.909 | 0.995 | 0.962 | 103 / 53 |
| -10 °C 0.5C | 0.862 | 0.990 | 0.922 | 168 / 90 |
| -20 °C 0.5C | 0.730 | 0.980 | 0.846 | 312 / 196 |
| 40 °C 2C | 0.978 | 1.009 | 1.030 | 64 / 33 |
| 10 °C 2C | 0.929 | 0.984 | 0.945 | 85 / 68 |
| 0 °C 2C | 0.882 | 0.956 | 0.867 | 109 / 119 |
| -10 °C 2C | 0.729 | 0.917 | 0.749 | 128 / 196 |
| -20 °C 2C | 0.492 | 0.843 | 0.522 | 138 / 321 |
* As pre-registered, L2 under-predicts the cold loss: R_Σ(−20 °C) = 0.165 Ω (spec 0.165) vs apparent 0.250 Ω; effective E_a 20.5 kJ/mol (spec 20.5) vs data 26.4. X9 (KiBaM: c 0.166, κ_K 40 h⁻¹, E_a,K 35 kJ/mol, R0+R1 kept) closes most of the 2C gap. At 0.5C it is still 0.04–0.12 too high (−20 °C: 0.85 vs 0.73).
* The data also contain a non-kinetic loss: the C/20 ratio is 0.97/0.95/0.93/0.92/0.90 at 40/10/0/−10/−20 °C (sequential-test aging and thermodynamic capacity). The model's C/20 ratio is 1.00. Dividing it out (fQ_obs_over_C20) leaves L2 +0.022 at 0 °C 0.5C (diagnostic only).
* 25 °C RMSE_V (CAL temperature): 16 mV (0.5C), 34 mV (2C). UDDS dynamic check: not run (no drive-cycle file in lg_hg2_raw).

## V3 NASA B0005 (aging -> TTE), V4 Oxford (pulsed load)
* V3: observed 2 A TTE 55.8 -> 39.7 min (-28.8 %); model 55.5 -> 39.6 min (-28.7 %; Q only -28.7 %, spec −28.5 %; R only -0.03 %). Both runs end on the `charge` event (z = 0 before 2.7 V; core_issues C2), so V3 confirms Q↓ only, not R↑ through the cut-off. RMSE_V 86/148 mV. 4 °C check (B0055): model Q ratio 1.00 vs 0.74 (no criterion; clearly off, same cause).
* V4: RMSE_V 82 mV for the R·Q-scaled cell ECM (R0·Q 89 mΩ·Ah). L1 gives 113 mV, so the RC branch helps. But OCV-only gives 33 mV: the Kokam power cell's measured 1-s step resistance is 21 mΩ vs the model's 84 mΩ. R·Q scaling (A12) does not transfer from energy cells to a power cell.

## V5 PyBaMM DFN (Ai2020, isothermal; ΔTTE % vs DFN; criterion ≤ 3 % at ≥ 0 °C)
| case | DFN h | L1 | L2 | X5 | X9 |
|---|---|---|---|---|---|
| 0.1C_25C | 10.736 | +0.4 | +0.4 | +0.4 | +0.1 |
| 0.5C_25C | 2.088 | +1.3 | +1.3 | +1.3 | -0.2 |
| 1C_25C | 1.007 | +2.4 | +2.4 | +2.3 | -0.7 |
| 0.1C_0C | 10.720 | +0.5 | +0.5 | +0.5 | -0.6 |
| 0.5C_0C | 2.008 | +4.2 | +4.2 | +4.2 | -1.3 |
| 1C_0C | 0.694 | +44.6 | +44.8 | +44.6 | +29.5 |
| 0.1C_-10C | 10.626 | +1.3 | +1.3 | +1.3 | -0.7 |
| 0.5C_-10C | 1.034 | +100.4 | +100.6 | +100.4 | +81.6 |
| 1C_-10C | 0.212 | +362.0 | +363.5 | +362.3 | +275.8 |
* At 25 °C every structure is within 2.4 % (X9 within 0.7 %). At 0 °C, 0.5C only X9 passes. At 1C (0 °C) and ≥ 0.5C (−10 °C) the DFN delivers only 1.67 / 1.21 of 2.47 Ah, a transport-limited collapse (mechanism not isolated here) that no ECM represents: structural error +30 % to +360 %. Sustained phone loads (≤ 0.21C, gaming) lie between the passing 0.1C and failing 0.5C cases at 0 °C.
* 10 s pulses (9 cases): RMSE_V L2 11–15 mV, X5 marginally lower, L1/X9 14–37 mV. R1 seen instantly overshoots the pulse drop; the 25/0 °C pulses are near in-sample for the [RUN] ECM fits. Local L2 = core LSODA to 2e-08.

## V6 published band (±25 %), V7 cold phenomenology
* V6: 11/16 pass ([PROTO] 11/16; model hours within 0.03 % of [PROTO]). Fails: 14P_audio_BT -25 %, Pixel8_local_video +57 %, Pixel8_idle_minbright +183 %, S24_local_video +53 %, S24_idle_minbright +54 %. Every iPhone case passes except BT audio; X14 (DSP audio) gives +21.6 %. Android screen-on floors are 1.6–2.9× ours (not tuned).
* V7 (pre-registered config: Q 2.7 Ah, R_Σ(20 °C) 0.25 Ω, ΔP_b 3 W, P̄ = E_u/159 min = 3.86 W): 0 °C TTE loss -51 % (LEE20 −15 %), s_off 59/68/78/100 % at 0/−5/−10/−15 °C (want 15–35 %), recovered 1429 mAh (LEE20 330): all three fail. Even at 25 °C the phone browns out at 14 %.
* V7 sensitivity (not fitted): with LEE20's own shut-off resistance 0.13 Ω (C5) the 0 °C loss is -11 %, s_off(−5 °C) 18 %, recovery 414 mAh, all inside the bands. The outcome hinges on the read-off phone resistance; the mechanisms (brown-out, recovery on warming) are reproduced.
