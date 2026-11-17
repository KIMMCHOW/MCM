# Smartphone Li-ion cell: electrical parameters and battery model families

Research note for 2026 MCM Problem A (smartphone battery drain). Scope: the electrical and thermal-electrical
behaviour of the phone cell and pack, the open-circuit-voltage (OCV) curve, equivalent-circuit (Thevenin)
parameters with their SOC and temperature dependence, rate and temperature effects on usable capacity, and the
continuous-time model families we could use. Load (screen, SoC, radio) power is out of scope here.

Compiled 2026-11-17. Every number below either comes from a source in the reference list (section e), with a
source key in brackets, or was **computed in this session** (key `[RUN]`) from open-licensed PyBaMM parameter
sets with the scripts in `docs/research/scripts/cell_ecm_pybamm_*.py`. The raw outputs are in
`docs/research/data/cell_ecm_pybamm/*.json`. Values that rest only on background knowledge are marked
**unverified** and given as ranges.

Conventions: SOC `z` in [0, 1]. Current `I > 0` on discharge. `T` is in kelvin inside the equations. `Q` is the
capacity in A·h (or C when multiplied by 3600).

Confidence levels:
- **H**: read directly in a primary or official document during this session.
- **M**: from a reputable secondary source, a search-result excerpt of a primary source, or computed from an open
  physics model.
- **L**: our own inference or synthesis, or unverified.

---

## (a) Parameter table

### a.1 Chemistry, voltages, capacity

| Name | Symbol | Value | Range | Unit | Conditions / notes | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Cathode | – | LiCoO2 (high-voltage, doped and coated) | – | – | LCO "continues to dominate" portable electronics | [LIN26] | H |
| LCO theoretical capacity | q_th | 274 | – | mAh/g | – | [LIN26] | H |
| LCO practical capacity at 4.2 V vs Li | q_4.2 | ~140 | – | mAh/g | about 51% of theoretical; higher cut-off voltages extract more Li | [LIN26] | H |
| LCO cathode energy density | – | 740 (at 4.45 V); 840 (at 4.55 V) | – | Wh/kg (cathode) | material level, not cell level | [LIN26] | H |
| Anode | – | graphite; Si/C–graphite blends in recent Android flagships | Si/C 800–1200 mAh/g | – | the full cell stays limited by LCO capacity | [LIN26]; press for products [GIZ25] | H / L |
| Max charge voltage, older phones | V_max | 4.35 | 4.2–4.35 | V | Galaxy S2/S3/S4 era; PinePhone pack 4.35 ± 0.025 V | [HOQ16], [PINE] | H |
| Max charge voltage, 2017+ | V_max | 4.40 | 4.35–4.40 | V | Galaxy S8 cell rated 4.4 V but derated to 4.35 V by Samsung | [QNOVO17] | M |
| Max charge voltage, current flagships | V_max | 4.45 | 4.40–4.48 | V | Consistent with the 4.45–4.48 V LCO cut-offs reported as commercial in reviews. **Not verified for any specific 2024–26 phone model.** | [LIN26] (material), unverified (product) | L |
| Nominal voltage vs V_max | V_nom | 3.80 at 4.35 V; 3.85 at 4.40 V | – | V | rule quoted by a battery-firmware company | [QNOVO17], [PINE] (3.80 V at 4.35 V) | M / H |
| Nominal voltage, current flagships | V_nom | 3.87 | 3.85–3.89 | V | iPhone 16 Pro Max reported at 3.875 V (press citing filings); Galaxy S24 Ultra listed at 3.88 V (aftermarket listings) | press [MR24], [PA24] | L–M |
| Cell discharge cut-off (spec) | V_cut,cell | 3.0 | 2.75–3.1 | V | rated-capacity test at 0.2C to 3.0 V ± 0.05 V | [PINE]; PyBaMM sets 3.0 / 3.105 V | H |
| Phone software cut-off ("SOC = 0") | V_sys | 3.4 | 3.2–3.5 | V | Qualcomm fuel-gauge device-tree **examples**: 3400 mV ("FG will steer SOC to zero"); 3300 mV ("minimum allowed system voltage. FG will report SOC = 0 at this voltage") | [QC-FG], [QC-SMB] | M (examples, not defaults) |
| Pack over-discharge protection | V_OD | 2.5 | 2.4–2.6 | V | release at 3.00 V; protection IC R5478N with 2 × 8205 MOSFETs | [PINE] | H |
| Pack over-charge protection | V_OC | 4.375 | 4.325–4.425 | V | for a 4.35 V pack | [PINE] | H |
| Pack over-current trip | I_OC | 6 | ≤8 | A | 2.75 Ah pack | [PINE] | H |
| Capacity: iPhone 16 / 16 Pro / 16 Plus / 16 Pro Max | Q | 3.561 / 3.582 / 4.674 / 4.685 | – | Ah | Brazil Anatel filings reported by press; Apple does not publish mAh | [MR24] | M |
| Capacity: Pixel 9 / 9 Pro / 9 Pro XL | Q | 4.700 / 4.700 / 5.060 typical; 4.558 / 4.558 / 4.942 minimum | – | Ah | official spec sheet | [PIX9] | H |
| Capacity: Galaxy S24 Ultra | Q | 5.000 typical; 4.855 rated | – | Ah | certification filing reported by press | [PA24] | M |
| Capacity: Si/C flagships | Q | 6.0–7.15 | – | Ah | OnePlus 13 6.0 Ah (805 Wh/L reported); Honor Magic V6 7.15 Ah | press [GIZ25] | L |
| Baseline capacity for our model | Q_n | 4.7 | 3.5–6.0 | Ah | synthesis of the rows above | – | L |
| Baseline energy | E_n | ≈18.2 (4.7 Ah × 3.87 V) | 13.5–23 | Wh | derived | – | L |
| Rated vs typical capacity gap | – | 2.3–3.0 | – | % | Pixel 9 Pro XL 2.3%, Pixel 9 3.0%, S24U 2.9% | derived from [PIX9], [PA24] | M |

### a.2 OCV curve

| Name | Symbol | Value | Range | Unit | Conditions / notes | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Equilibrium OCV table, LCO/graphite pouch (Enertech, 4.2 V) | U(z) | Appendix A1 (21 points) | 3.000–4.200 | V | 25 °C. Computed from the electrode OCPs (data from Ai et al. 2020 / Rieger 2016) and the PyBaMM electrode-SOH stoichiometry window (x: 0.0051→0.848; y: 0.966→0.430) | [RUN] on [PYB-AI20] | M |
| Parametric OCV in Chen–Rincón-Mora form, fitted to the above | b0, b1, a0..a3 | 0.60117, 36.426, 3.6193, 0.44311, −0.48605, 0.63127 | RMSE 7.0 mV; max 18 mV | V, – | `U = −b0·exp(−b1 z) + a0 + a1 z + a2 z² + a3 z³` | [RUN] | M |
| Plett "combined" model fitted to the same OCV, z ∈ [0.02, 0.98] | K0..K4 | 3.45064, −0.0122655, 0.506218, −0.110213, −0.0723812 | RMSE 14.4 mV; max 58.7 mV | V | worse than the Chen–RM form near both ends | [RUN] | M |
| Original Chen–Rincón-Mora polymer Li-ion OCV | – | U = −1.031·exp(−35 z) + 3.685 + 0.2156 z − 0.1178 z² + 0.3201 z³ | 2.65–4.10 | V | Coefficients cross-checked against a GSEIM implementation that cites [CRM06]. The full paper could not be retrieved (TLS chain failure), so the exact cell model is unverified. | [CRM06] via [GSEIM] | M |
| Mean OCV over z (Ai2020 fit) | Ū | 3.820 | – | V | 4.2 V-class LCO/graphite | [RUN] | M |
| High-voltage stretch (our construction) | ΔV, p | 0.242, 3.86 (4.45 V / 3.87 V) | (0.192, 5.44) for 4.40/3.85; (0.272, 2.90) for 4.48/3.89 | V, – | `U_HV(z) = U_4.2(z) + ΔV·z^p`. Matches the target V_max and the target mean OCV, and leaves z < 0.3 almost unchanged. **Assumption; no measured 4.45 V OCV curve was found.** | [RUN] | L |
| OCV slope dU/dz | – | 1.23 (z = 0.1); 0.51 (0.2); 0.36 (0.5); 0.89 (0.8); 1.12 (0.95) | – | V per unit SOC | needed for the sensitivity of time-to-empty (TTE) to the cut-off voltage | [RUN] | M |
| OCV hysteresis, LCO/graphite | M_h | **unverified**, probably ≤ ~20 mV mid-SOC | 0–40 | mV | LCO is a single-phase layered oxide. Hysteresis in graphite cells is "primarily related to the graphite electrode" (staging peaks); that paper measured NMC(-LCO)/graphite cells. | [OVEJ19] (mechanism) | L |

### a.3 Equivalent-circuit (Thevenin) parameters

`[RUN]` values are fitted to PyBaMM DFN (Doyle–Fuller–Newman, isothermal) responses to a 1C, 60 s pulse followed by 1 h of rest. A
coulomb-counted OCV drift correction is applied. PyBaMM sets have **zero contact resistance**, so these are
electrode-only values.

| Name | Symbol | Value | Range | Unit | Conditions / notes | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Ohmic+fast resistance × capacity | R0·Q | 94.7 (Ai2020, 41.5 mΩ × 2.28 Ah); 84.1 (Marquis2019, 123.6 mΩ × 0.68 Ah) | 80–100 | mΩ·Ah | 25 °C, z = 0.5, 1RC fit | [RUN] | M |
| Cell R0 scaled to a 4.7 Ah phone cell | R0,cell | 19 | 15–35 | mΩ | `R0·Q / Q`. Lower bound for real cells: thick, energy-optimised electrodes, tabs and contacts are not in the model. | [RUN] + scaling (L) | L–M |
| Pack AC impedance spec (1 kHz) | Z_1k | ≤ 100 | – | mΩ | 2.75 Ah phone pack incl. protection board, 50% SOC, 23 °C | [PINE] | H |
| Protection-board (PCM) conduction resistance | R_PCM | 40 | 20–60 | mΩ | 2 × 8205 dual MOSFET (SOT-23-6), min/typ/max | [PINE] | H |
| Modern dual-FET protection switch | R_SS(on) | 3.2 (EFC8811R, 12 V); 7 (ECH8693R, 24 V, smartphone/tablet target) | 3–45 | mΩ | product-page titles; gate-drive dependent | [ONSEMI] | M |
| Legacy phone pack breakdown | – | cell 50; weld 1; PTC 25 (18–30); protection PCB 50; total ≈ 130 | – | mΩ | Siemens AG (2015) table quoted by Battery University | [BU802a] | L–M |
| **Recommended pack R0 (cell + PCM + connector), new, 25 °C, mid-SOC** | R0,pack | **60** | 35–130 | mΩ | synthesis: 15–35 (cell) + 5–60 (PCM) + a few (flex, connector, PTC if present) | – | L |
| Polarization R1 (1RC) | R1 | Ai2020 14.4 mΩ (R1·Q = 33 mΩ·Ah); Marquis2019 18.9 mΩ (12.8 mΩ·Ah) | for 4.7 Ah: 3–10 (model); real cells likely 5–25 | mΩ | 25 °C, z = 0.5 | [RUN] | M / L |
| Time constant (1RC) | τ1 = R1·C1 | 35 (Ai2020); 29 (Marquis2019) | 25–40 at 25 °C | s | z = 0.5 | [RUN] | M |
| 2RC split (Ai2020, 25 °C, z = 0.5) | R1, τ1; R2, τ2 | 2.5 mΩ, 6.7 s; 14.1 mΩ, 47 s | – | mΩ, s | RMSE 0.07 mV (2RC) vs 0.32 mV (1RC) | [RUN] | M |
| Original Chen–RM resistances, high-SOC limits | R_s,∞; R_TS,∞; R_TL,∞ | 74.46; 46.69; 49.84 | – | mΩ | `R_s(z) = 0.07446 + 0.1562·exp(−k z)` Ω; the high-SOC limits assume the same `a + b·exp(−k z)` form for all three. a0 and b0 verified via [GSEIM]; the exponent k ≈ 24.37 and the C_TS/C_TL coefficients are **unverified**. | [CRM06] via [GSEIM] | M / L |
| R0 rise at low SOC (25 °C, 0.1 s pulse) | R0(z)/R0(0.5) | Ai2020: 1.007 (z = 0.8), 1.19 (0.2), 1.39 (0.1). Marquis2019: 0.96, 1.15, 1.27 | – | – | 1C pulse | [RUN] | M |
| SOC-shape fit | A, z0 | A = 0.80, z0 = 0.139 | – | –, – | `R0(z) = R0,∞·[1 + A·exp(−z/z0)]`, fitted to the Ai2020 points at z = 0.1 and 0.2 | [RUN] | M |
| Activation energy, R0 (ohmic + fast) | Ea,R0 | **18** | 15–25 | kJ/mol | DFN Marquis2019 25 → −20 °C: 18.0. Barcellona LCO 10 Ah, ac resistance, 20–46 °C: 15.6 (derived). Low-T review: R_ohm(−20 °C) = 278.6% of R_ohm(20 °C), giving 15.8 (derived; number from a search excerpt). | [RUN], [BARC22], [SUN23] | M |
| Activation energy, R1 (charge transfer + diffusion) | Ea,R1 | **30** | 25–55 | kJ/mol | DFN Marquis2019: 31.4. Ai2020 R1 (25/0 °C): 27.7. Barcellona low-frequency R at 25 °C: 29.2 (effective), with an exponential part equal to 55 kJ/mol at 25 °C. PyBaMM exchange-current Ea: 35 (Chen2020 graphite), 37.5 / 39.6 (Marquis graphite / LCO), 53.4 / 43.6 (Ecker2015). | [RUN], [BARC22], [PYB] | M |
| R0 ratio, −20 °C vs 25 °C | – | 3.6 (Marquis2019, z = 0.5) | 2.8–4 | – | R1 ratio 9.5; τ1 ratio 4.5 (29 s → 131 s) | [RUN] | M |
| Barcellona LCO resistance law | – | `R(T) = a·exp(−0.075·T[°C]) + c`; new cell, 50% SOC, low frequency: a = 15.2 mΩ, c = 2.1 mΩ. Measured 5.40 mΩ at 20 °C and 2.59 mΩ at 46 °C. | 20–46 °C only | mΩ | 10 Ah LCO power cell (General Electronics Battery 8773160K, 4.2/2.75 V) | [BARC22] | H |
| Low-T dominant resistance | – | R_ct "most significantly increased and becomes predominant" below −10 °C | – | – | symmetric-cell study | [ZXJ03] | H (qualitative) |
| Resistance growth with ageing | – | +59% (20 °C), +33% (46 °C) after 20,000 Ah throughput at C/C0 = 0.972 | – | – | 10 Ah LCO; end of life is commonly defined as 80% capacity **or** 200% resistance | [BARC22] | H |
| Coulombic efficiency | η | ≈ 0.99 on charge; 1 on discharge | – | – | energy efficiency ≈ 95% | [PLETT-N] | M |

### a.4 Capacity versus C-rate and temperature

| Name | Symbol | Value | Range | Unit | Conditions / notes | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Rate effect at 25 °C, Q(rate)/Q(0.1C) | – | 0.2C: 0.998 / 0.997. 0.5C: 0.990 / 0.991. 1C: 0.978 / 0.981. 2C: 0.943 / 0.962. | – | – | Ai2020 / Marquis2019 DFN, isothermal, to the cell cut-off | [RUN] | M |
| Effective Peukert exponent at 25 °C | k | 1.01 (0.1→1C); 1.04 (0.5→2C) | 1.0–1.05 | – | derived from the row above | [RUN] | M |
| Peukert fit, Dualfoil Li-ion (2001) | a, b | a = 1550.2, b = 1.09 | – | (A/m², min) | 1–12 h discharges; Voc 4.31 V, cut-off 3.19 V. Peukert max/mean error 21% / 6% vs 5% / 2% for the diffusion model. | [RV01] | H |
| Peukert exponent, Itsy Li-ion (Dualfoil) | b | ≈ 1.04 | – | – | derived by us from the constant-load rows T1 and T10 of the Rakhmatov et al. table reproduced in [JH09] | derived from [JH09] | M |
| Datasheet minimum capacity at 0 °C | Q(0 °C)/Q_n | ≥ 0.85 | – | – | 0.2C to 3.0 V after 4 h soak | [PINE] | H |
| Datasheet minimum at −20 °C | – | ≥ 1 h at 1.5 A (≥ 0.55·Q_n) | – | – | 2.75 Ah pack, 0.55C to 3.0 V; "can boot normally" | [PINE] | H |
| Datasheet minimum at 60 °C | – | ≥ 0.95 | – | – | 0.2C | [PINE] | H |
| DFN capacity vs T (Marquis2019, isothermal), Q/Q(25 °C, same rate) | – | 0.2C: 0.997 (10 °C), 0.994 (0), 0.989 (−10), 0.980 (−20). 0.5C: 0.994, 0.987, 0.974, 0.753. 1C at −20 °C: 0.27. | – | – | no self-heating. The steep drop is a **cut-off effect**: polarization grows until V hits the cut-off. | [RUN] | M |
| DFN capacity vs T (Ai2020, isothermal) | – | 0.5C at 0 °C: 0.993; 1C at 0 °C: 0.772; 0.5C at −10 °C: 0.545 | – | – | below −10 °C the solver failed or gave non-physical results, so not used | [RUN] | M–L |
| Charging temperature limits | – | 0–15 °C: ≤ 0.2C; 15–50 °C: ≤ 0.5C. Discharge allowed −10 to 55 °C. | – | – | pack spec | [PINE] | H |

### a.5 Thermal, entropic, self-discharge, ageing, standby

| Name | Symbol | Value | Range | Unit | Conditions / notes | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Entropic coefficient (cell) | ∂U/∂T | −0.25 at z = 0.5; table in Appendix A2 | −0.05 (z = 1) to −0.71 (z = 0.05) | mV/K | From Ai2020 electrode entropic fits (Enertech LCO/graphite). The graphite fit is documented as accurate only for 0.43 < sto < 0.9936, so the z < 0.45 values are extrapolated. | [RUN] on [PYB-AI20] | M–L |
| LCO vs other cathodes | – | LiCoO2 has a "much larger entropy change" than NMC and LFP | – | – | reversible heat is a significant fraction of total heat | [VISW10] | H (qualitative) |
| Cell-average density, specific heat | ρ, c_p | 2470; 1080.2 | – | kg/m³; J/(kg·K) | Ai2020 set (2.28 Ah, 15.3 cm³ → m ≈ 38 g, m·c_p ≈ 41 J/K) | [PYB-AI20] | M |
| Cell heat capacity per Ah (derived) | C_th/Q | ≈ 18 | 15–22 | J/(K·Ah) | gives ≈ 85 J/K for a 4.7 Ah cell alone. The phone body adds much more (outside this topic). | derived | L |
| Self-discharge (spec bound) | – | ≥ 85% retained after 28 d at 23 °C (≥ 90% recoverable) | – | – | upper bound, about 15%/month. Real cells are typically far lower (**unverified**: ~1–3%/month). | [PINE] | H (bound) |
| Protection-IC quiescent current | I_q | 4 | ≤ 7 | µA | at Vcc = 3.9 V | [PINE] | H |
| Phone standby current (airplane mode, idle) | – | within ~10 mA | – | mA | Galaxy S-series measurements | [HOQ16] | M |
| Cycle-life design targets | – | 80% capacity at 500 cycles (iPhone 14 and earlier); 80% at 1000 cycles (iPhone 15 and later) | – | – | "ideal conditions" | [APPLE] | H |
| Cycle-life spec, generic phone pack | – | ≥ 80% after 1000 cycles at 0.2C/0.2C, 23 °C | – | – | – | [PINE] | H |
| KiBaM, Li-ion (Itsy pocket computer) | c, k′ | c = 0.166; k′ = 0.122 min⁻¹ (k = k′·c(1−c) = 0.0169 min⁻¹) | – | –, min⁻¹ | KiBaM least-squares fitted to the RV diffusion model with β = 0.273 min^−½ at I = 1 A | [JH09] | H (as reported) |
| Rakhmatov–Vrudhula diffusion model, Itsy | α, β | α ≈ 40,375 mA·min (≈ 673 mAh); β = 0.273 min^−½ | – | – | α derived by us as I × L_ideal from the table in [JH09] (222.7 mA × 181.3 min) | [JH09] | M |
| Rakhmatov–Vrudhula diffusion model, Dualfoil Li-ion (2001) | α, β | 271.47; 10.39 | – | per-area units (A·min/m²; β = w/√D) | different β convention from [JH09] | [RV01] | H |

---

## (b) Key equations and model forms

**1. SOC (coulomb counting) with a constant-power load.** Phones draw roughly constant *power*, so current
depends on the terminal voltage:

```
dz/dt = −η(I) · I(t) / (3600 · Q_eff(T, age))       η = 1 on discharge, ≈ 0.99 on charge   [PLETT-N]
I(t)  = P_load(t) / V(t)
```

**2. Thevenin ECM (1RC or 2RC), with a hysteresis state as an option.** [PLETT-N], [CRM06]:

```
dv_k/dt = −v_k / (R_k C_k) + I / C_k                 k = 1 (τ ≈ 5–10 s), 2 (τ ≈ 30–50 s at 25 °C)
V = U(z, T) + h − I·R0(z, T) − Σ_k v_k
dh/dt = −|η I γ / (3600 Q)| · (h + M·sgn(I))         (Plett one-state hysteresis; M ≈ 0 for LCO baseline)
```

Combining (1) and (2) with V = P/I gives a quadratic for I:

```
R0·I² − (U − Σv_k)·I + P = 0
I = [ (U − Σv) − sqrt( (U − Σv)² − 4·R0·P ) ] / (2·R0)
```

A real solution exists only if `P ≤ P_max = (U − Σv)² / (4·R0)`. The phone browns out when either
`V < V_sys` (about 3.4 V) or `P > P_max`. P_max shrinks at low z (U falls and R0 rises), at low T (Arrhenius),
and with age. This is the "unexpected shutdown" mechanism Apple describes [APPLE].

**3. OCV(z, T).**

```
Chen–Rincón-Mora form [CRM06]: U(z) = −b0·exp(−b1 z) + a0 + a1 z + a2 z² + a3 z³
Plett "combined" [PLETT04]:     U(z) = K0 − K1/z − K2 z + K3 ln z + K4 ln(1−z)
     (sign convention as commonly quoted; the original was not read here. Our fitted K's in a.2 use
      U = K0 + K1/z + K2 z + K3 ln z + K4 ln(1−z), so only the signs of K1 and K2 differ)
"Combined+3" [PILLAI23]:        U(s) = k0 + k1/s + k2/s² + k3/s³ + k4/s⁴ + k5 s + k6 ln s + k7 ln(1−s)
     with scaled SOC s = (1 − 2ε) z + ε, ε = 0.175 [AHMED20]
Temperature [PLETT-N]:          U(z, T) = U0(z) + T·U_rel(z)  ≡  U(z, Tref) + (T − Tref)·∂U/∂T(z)
```

**4. Arrhenius temperature dependence and SOC dependence of the resistances.**

```
R_j(z, T) = R_j,ref(z) · exp[ (Ea,j / R_g)·(1/T − 1/T_ref) ]      Ea,R0 ≈ 18 kJ/mol, Ea,R1 ≈ 30 kJ/mol
τ_j(T) similarly (Ea,τ1 ≈ 21 kJ/mol from [RUN])
R0,ref(z) = R0,∞ · [1 + A·exp(−z/z0)],   A ≈ 0.8, z0 ≈ 0.14            [RUN]
Alternative [BARC22]: R(T) = a·exp(−b·T[°C]) + c, b = 0.075 /°C (valid 20–46 °C)
Usable capacity: Q_eff(T) = Q_n·f_Q(T), with f_Q from datasheet points (0 °C ≥ 0.85, −20 °C ≥ 0.55 at 0.55C)
  or emergent from the cut-off when R(T) is modeled (preferred: physics, not curve fitting)
```

**5. Lumped thermal model with Bernardi heat generation** [BERN85] (the ECM form used by PyBaMM's Thevenin model,
which has an "Entropic change" term):

```
C_th dT/dt = I²·R0 + Σ_k v_k²/R_k − I·T·∂U/∂T(z) − (T − T_amb)/R_th
             (irreversible Joule + polarization)  (reversible; > 0 on discharge when ∂U/∂T < 0)
```

At I = 1 A, T = 298 K and ∂U/∂T = −0.25 mV/K, the reversible heat is about 0.075 W. Joule heat at 60 mΩ is
0.06 W. **At phone currents the two terms are of the same order** [RUN-derived].

**6. Peukert** [DS06], [RV01]: `t_life = a / I^b` (b ≈ 1.01–1.09 for Li-ion at room temperature). It is only valid
at constant current and constant temperature.

**7. Kinetic Battery Model (KiBaM)** [MM93], [JH09]. Two wells: available charge y1 (fraction c of the capacity)
and bound charge y2.

```
dy1/dt = −I + k(h2 − h1),   dy2/dt = −k(h2 − h1),   h1 = y1/c,  h2 = y2/(1−c)
equivalently  dγ/dt = −I,  dδ/dt = I/c − k′δ,  k′ = k/(c(1−c));  empty when γ = (1−c)·δ
constant I:  δ(t) = (I/(c k′))·(1 − e^{−k′t})
```

**8. Rakhmatov–Vrudhula diffusion model** [RV01], [JH09]. Constant-current form with β = π√D/w ([JH09]
convention):

```
α = I·[ L + 2 Σ_{m≥1} (1 − e^{−β² m² L}) / (β² m²) ]          (battery empty at t = L)
```

KiBaM is an approximation of this model. Matching the first-order term gives c = 1/3 and k′ = β², which is a
poor approximation; the least-squares fit at I = 1 A gives c = 0.166 and k′ = 0.122 min⁻¹ for β = 0.273 min^−½,
with errors up to 80% for t < 10 min [JH09].

**9. Single-particle model (SPM / SPMe)** [MQ19]. Solid diffusion in a representative sphere for each electrode:

```
∂c_s/∂t = (D_s/r²) ∂/∂r (r² ∂c_s/∂r),   −D_s ∂c_s/∂r |_{r=R} = j = ± I / (a_s F A L)
V = U_p(c_p,surf/c_p,max) − U_n(c_n,surf/c_n,max) − (2RT/F)[asinh(j_p/2j0,p) + asinh(j_n/2j0,n)] − I·R_contact
j0 ∝ k_ref·exp[−Ea/R (1/T − 1/Tref)]·sqrt(c_e c_s,surf (c_max − c_s,surf))
```

The SPMe adds one PDE for the electrolyte and an algebraic voltage correction [MQ19]. DFN is the full
pseudo-two-dimensional (P2D) model.

### Which family suits a phone, and why

- **Ideal coulomb counting** gives no voltage, so it cannot express shutdown, cold, or ageing effects on power.
  It is only a baseline.
- **Peukert** gives a single algebraic lifetime. It fails under variable load and temperature [DS06], [RV01], and
  phone C-rates are small (our estimate: about 0.05–0.5C) where k ≈ 1.01. Not adequate.
- **KiBaM and the RV diffusion model** capture rate-capacity and recovery effects with 2–3 parameters, and are
  continuous-time and physically motivated (two charge wells, or a diffusion PDE). They still have no OCV, no
  resistance, and no temperature dependence, so they do not predict the 3.4 V brown-out. Useful as an
  interpretable "available charge" sub-model.
- **Thevenin 1RC/2RC ECM with Arrhenius R(T), SOC-dependent R0(z) and OCV(z, T), coupled to a lumped thermal ODE
  and a constant-power load** is our choice. It is a 3–4-state ODE system. Each element has a physical meaning
  (ohmic, charge-transfer, diffusion polarization; Butler–Volmer and Arrhenius kinetics). It is cheap enough for
  Monte Carlo uncertainty quantification and sensitivity analysis. It captures the cut-off-limited TTE, cold
  shutdown and ageing (Q↓, R↑). The RC states play the role of KiBaM's bound charge.
- **SPMe / DFN (PyBaMM)** are too heavy and too parameter-hungry to serve as the main model. They are the
  **high-fidelity physics reference** for checking the reduced ECM (rate, temperature and pulse behaviour), as we
  did in `[RUN]`.

---

## (c) Key findings

1. **The phone cell is a high-voltage LCO/graphite pouch.** Typical values are V_max ≈ 4.40–4.48 V,
   V_nom ≈ 3.85–3.89 V and Q ≈ 3.5–5.1 Ah (iPhone 16 family 3.56–4.69 Ah; Pixel 9 4.7–5.06 Ah;
   S24U 5.0 Ah). Si/C anodes push recent Android flagships to 6–7 Ah. Specs come as "typical" vs
   "minimum/rated" capacity, about 2.5–3% apart.
2. **The phone shuts down well above the cell's 3.0 V cut-off.** Fuel gauges map SOC = 0 to a *loaded* system
   voltage of about 3.3–3.4 V (Qualcomm DT examples). Pack protection only trips at 2.5 V. Because the OCV curve is
   steep at low SOC (dU/dz ≈ 1.2 V per unit at z = 0.1), the usable charge is set by U(z) − I·R ≥ V_sys. A
   physical ECM must model the voltage, not just charge.
3. **The resistance budget is dominated by the pack, not the electrodes.** Open DFN models give a cell R0·Q of
   about 85–95 mΩ·Ah (≈ 19 mΩ for a 4.7 Ah cell). The protection board alone is 20–60 mΩ in a real phone pack, and
   the pack spec is ≤ 100 mΩ at 1 kHz. Legacy breakdowns total ≈ 130 mΩ. Baseline R0,pack = 60 mΩ
   (35–130 mΩ) should be a primary sensitivity parameter.
4. **Rate-capacity (Peukert) effects are negligible at room temperature for phone loads.** At 0.2–0.5C, DFN loses
   only 0.2–1% of capacity versus 0.1C (effective k ≈ 1.01). This is one of the "surprisingly small" factors.
   KiBaM-style recovery matters only for high-current bursts near empty.
5. **Cold is a voltage problem, not lost lithium.** In isothermal DFN at 0.2C, capacity to the cell cut-off falls
   only 2% at −20 °C. At 0.5C it falls 25%, and at 1C 73%. R0 rises ×3.6 and R1 ×9.5 from 25 to −20 °C
   (Ea ≈ 18 and 31 kJ/mol). With a 3.4 V system cut-off and pulsed loads, P_max collapses, which is the
   unexpected-shutdown mechanism described by Apple. The datasheet guarantees ≥ 85% at 0 °C (0.2C) and
   ≥ 55% at −20 °C (0.55C).
6. **SOC dependence of R0 matters only at the end.** It is about 0–5% between z = 0.3 and 1, +15–20% at z = 0.2,
   and +27–40% at z = 0.1. This compounds finding 2.
7. **The Chen–Rincón-Mora OCV form fits LCO/graphite well.** RMSE is 7 mV over the full range, versus 14 mV
   (max 59 mV) for Plett's combined model on the same data. Use the Chen–RM form (or a table), with a z^p stretch
   for 4.45 V cells.
8. **Reversible (entropic) heat is not negligible for LCO.** ∂U/∂T is −0.05 to −0.7 mV/K, and LCO has the
   largest entropy change among common cathodes [VISW10]. At 1 A it is comparable to Joule heat, though the total
   is only about 0.1–0.2 W (our estimate at ~1 A), so phone self-heating from the battery itself is small next to
   the SoC and display.
9. **Hysteresis is small and comes mainly from graphite staging for LCO.** The baseline can ignore it and treat it
   as an OCV uncertainty of ±10–20 mV (unverified magnitude).
10. **Ageing enters as Q↓ and R↑.** The design targets are 80% capacity at 500 cycles (iPhone ≤ 14) and 1000
    cycles (iPhone 15+). Resistance grew 33–59% in a 10 Ah LCO cell at only 2.8% capacity loss, so for phones the
    power fade (R) can matter before the capacity fade does. The end-of-life convention is 80% Q or 200% R.
11. **PyBaMM is a usable physics reference.** It is pip-installable (`pip install pybamm`; 26.10.1.1 is installed
    here and is the latest on PyPI as of 2026-11-17) and BSD-3-licensed. Its LCO/graphite pouch sets are
    **Ai2020** (Enertech 2.28 Ah pouch, measured OCP and entropic data; best match to a phone cell, but its kinetic
    Ea values are only 5 kJ/mol), **Marquis2019** (Kokam SLPB78205130H geometry, Dualfoil LCO/graphite kinetics,
    realistic Ea 18.5–42.8 kJ/mol) and **Ramadass2004** (LCO/graphite, 1 Ah, flagged "Frankenstein" by PyBaMM).
    PyBaMM also ships a Thevenin ECM with a lumped thermal model and an entropic term (`ECM_Example`, 100 Ah, as a
    template).

---

## (d) Caveats

- **No open dataset of a 2024–26 flagship phone cell was found.** The 4.45 V OCV curve is constructed (z^p
  stretch), not measured. Absolute R0 for a modern phone pack is a synthesis from a 2.75 Ah PinePhone datasheet
  (4.35 V generation), a legacy breakdown, and DFN scaling. Treat it as uncertain by ×2.
- **`[RUN]` values come from PyBaMM's isothermal DFN** with literature parameter sets that PyBaMM itself
  says "do not claim to be representative". Dualfoil-era kinetics in Marquis2019 are idealised, there is no
  contact resistance, and there is no self-heating. Ai2020 solver runs failed or were non-physical below −10 °C.
  Ai2020 uses Ea = 5 kJ/mol for kinetics and diffusion, so its temperature trend is unrealistically weak; its R1
  trend comes mostly from the electrolyte.
- **The PyBaMM OCV windows are 4.2 V (Ai2020) and 4.1 V (Marquis2019).** Marquis2019's OCV is a Dualfoil fit
  with an unusual plateau near 3.69 V, so we prefer Ai2020 for OCV.
- **The entropic table uses the Ai2020 graphite fit outside its stated validity (sto < 0.43)** for z < ~0.45.
- **The Chen & Rincón-Mora (2006) paper could not be downloaded.** The rincon-mora.gatech.edu certificate chain
  failed under the proxy, and the archive mirrors were blocked. Only the OCV coefficients and the asymptotic
  resistances were cross-checked through the GSEIM implementation. The remaining coefficients (exponents,
  capacitances) and the exact cell identity are unverified, and none of those are used.
- **Some figures come from press or aftermarket sources** (iPhone mAh via Anatel reports, S24U rated capacity,
  nominal voltages of 3.875 or 3.88 V, Si/C phones). The 278.6% R_ohm ratio [SUN23] comes from a search-engine
  excerpt; MDPI full texts were blocked (HTTP 403).
- **The Qualcomm cut-off values are device-tree *examples*, not defaults.** OEMs tune them. Apple's thresholds are
  undisclosed.
- **Hysteresis magnitude for LCO/graphite and typical real self-discharge rates** are unverified.
- **The Rakhmatov–Vrudhula and KiBaM parameters belong to a 2000s Itsy / Dualfoil Li-ion cell.** Their absolute
  values do not transfer to modern phones; only the structure does.
- **Licensing.** PyBaMM is BSD-3; Ai2020 [AI20] and PyBaMM [SULZER21] are CC BY; [LIN26], [OVEJ19], [SUN23] and
  [BARC22] (post-print) are open access; the Oxford degradation dataset is ODbL. The PinePhone datasheet, Google
  spec sheet, onsemi pages and Qualcomm/AOSP device-tree docs are freely accessible but carry no explicit open
  data licence (AOSP kernel docs are GPL-2.0). Cite them as published specifications.

---

## (e) References

1. **[LIN26]** Z. Lin, Y. Ying, H. Li, Y. Ren, T. Liu, P. Hou, H. Huang, "Decoding High-voltage LiCoO2: From
   Degradation to Stabilization Toward Durable Li-ion Batteries," *Advanced Materials* 38 (2026).
   doi:10.1002/adma.202523570. Open access, CC BY 4.0. Read in full (PolyU repository PDF).
2. **[PIX9]** Google, "Pixel 9, Pixel 9 Pro, Pixel 9 Pro XL Spec Sheet" (2024), footnote 6.
   https://cdn.uc.assets.prezly.com/05fc1248-3c32-4f54-b0ce-f8752f25e6f4/-/inline/no/Pixel%209,%20Pixel%209%20Pro,%20Pixel%209%20Pro%20XL%20Spec%20Sheet.pdf.
   Official, freely accessible.
3. **[MR24]** MacRumors, "iPhone 16 Battery Capacities Revealed," 19 Sep 2024 (citing Anatel filings).
   https://www.macrumors.com/2024/09/19/iphone-16-battery-capacities-revealed/. Values seen in search excerpts
   from several outlets. Press.
4. **[PA24]** PhoneArena / SamMobile reports on Galaxy S24 Ultra EB-BS928ABY certification (4,855 mAh rated,
   5,000 mAh typical). https://www.phonearena.com/news/bad-samsung-galaxy-s25-ultra-battery-capacity-news-gets-certified_id161550
   and https://sammobile.com/news/galaxy-s24-plus-ultra-batteries-pictured. Press, via search excerpts.
5. **[GIZ25]** Gizmochina, "Top 10 latest smartphones with silicon carbon batteries in 2025" (7 Feb 2025).
   https://www.gizmochina.com/2025/02/07/top-10-latest-smartphones-with-silicon-carbon-batteries-in-2025/.
   Press, via search excerpts; low confidence.
6. **[QNOVO17]** Qnovo, "What the label on a lithium-ion battery means" (4 May 2017).
   https://qnovo.com/113-label-lithium-ion-battery-means. Company technical blog, freely accessible.
7. **[PINE]** ShenZhen SaiJiaoYang Energy & Science Technology Co., "Lithium battery specification, model
   QZ01-396172, 2750 mAh" (PinePhone battery), rev. A0, incl. PCM parameter table, p. 12.
   https://xnux.eu/devices/feature/bat-pp/PinePhone_Battery_model_QZ01-396172-2750.pdf. Manufacturer datasheet,
   freely accessible; no explicit licence.
8. **[HOQ16]** M. A. Hoque, M. Siekkinen, J. Koo, S. Tarkoma, "Accurate Online Full Charge Capacity Modeling of
   Smartphone Batteries," arXiv:1604.05689v2 (2016). https://arxiv.org/abs/1604.05689. arXiv licence.
9. **[QC-FG]** Qualcomm / AOSP kernel msm, "qpnp-fg.txt" device-tree binding (property
   `qcom,fg-cutoff-voltage-mv`, example 3400).
   https://android.googlesource.com/kernel/msm/+/973f4134d9deb396415846f902848f0a32cb4cfa/Documentation/devicetree/bindings/power/qpnp-fg.txt.
   GPL-2.0 kernel documentation.
10. **[QC-SMB]** Qualcomm / AOSP kernel msm, "smb1360-charger-fg.txt" (cut-off example 3300 mV, "minimum allowed
    system voltage"). https://android.googlesource.com/kernel/msm/+/973f4134d9deb396415846f902848f0a32cb4cfa/Documentation/devicetree/bindings/power/smb1360-charger-fg.txt.
    GPL-2.0.
11. **[APPLE]** Apple Support, "iPhone battery and performance" (article 101575 / HT208387).
    https://support.apple.com/en-us/101575. Official, freely accessible.
12. **[CRM06]** M. Chen, G. A. Rincón-Mora, "Accurate electrical battery model capable of predicting runtime and
    I–V performance," *IEEE Trans. Energy Conversion* 21(2), 504–511 (2006). doi:10.1109/TEC.2006.874229.
    Author PDF https://rincon-mora.gatech.edu/publicat/jrnls/tec05_batt_mdl.pdf (could not be retrieved:
    TLS chain error). Metadata verified via Crossref.
13. **[GSEIM]** GSEIM circuit simulator documentation, battery subcircuit `s_battery_1` (lists a0_vsrc = 3.685,
    a1 = 0.2156, a2 = −0.1178, a3 = 0.3201, b0_vsrc = −1.031, b1 = 35, a0_rs = 0.07446, b0_rs = 0.1562,
    a0_rts = 0.04669, a0_rtl = 0.04984; cites [CRM06]).
    https://gseim.github.io/build/html/_downloads/a8f47b624b4e8ffabf95c823df4b8cb7/sub_s_battery_1.pdf.
    Seen in search excerpts only (direct fetch returned 404).
14. **[PLETT04]** G. L. Plett, "Extended Kalman filtering for battery management systems of LiPB-based HEV battery
    packs: Part 2. Modeling and identification," *J. Power Sources* 134, 262–276 (2004).
    doi:10.1016/j.jpowsour.2004.02.032. Metadata via Crossref; the combined-model form is attributed to Plett by
    [PATT14] and [PILLAI23].
15. **[PLETT-N]** G. L. Plett, *ECE4710/5710 Lecture Notes, Ch. 2 "Equivalent-Circuit Cell Models"*, Univ. of
    Colorado Colorado Springs (© 2011–2018). http://mocha-java.uccs.edu/ECE5710/ECE5710-Notes02.pdf. Freely
    accessible; copyright author.
16. **[PATT14]** B. Pattipati, B. Balasingam, G. V. Avvari, K. R. Pattipati, Y. Bar-Shalom, "Open circuit voltage
    characterization of lithium-ion batteries," *J. Power Sources* 269, 317–333 (2014).
    doi:10.1016/j.jpowsour.2014.06.152.
17. **[PILLAI23]** P. Pillai, J. Nguyen, B. Balasingam, "Performance Analysis of Empirical Open-Circuit Voltage
    Modeling in Lithium Ion Batteries, Part-2: Data Collection Procedure," arXiv:2306.16547 (2023).
    https://arxiv.org/abs/2306.16547. arXiv licence; 4 Ah test cell, C/64 OCV test.
18. **[AHMED20]** M. S. Ahmed, S. A. Raihan, B. Balasingam, "A scaling approach for improved state of charge
    representation in rechargeable batteries," *Applied Energy* 267, 114880 (2020).
    doi:10.1016/j.apenergy.2020.114880.
19. **[AI20] / [PYB-AI20]** W. Ai, L. Kraft, J. Sturm, A. Jossen, B. Wu, "Electrochemical Thermal-Mechanical
    Modelling of Stress Inhomogeneity in Lithium-Ion Pouch Cells," *J. Electrochem. Soc.* 167, 013512 (2020).
    doi:10.1149/2.0122001JES. Open access (CC BY). Implemented as PyBaMM parameter set "Ai2020"
    (LCO data from Rieger et al. 2016).
20. **[MQ19] / [PYB-MQ19]** S. G. Marquis, V. Sulzer, R. Timms, C. P. Please, S. J. Chapman, "An Asymptotic
    Derivation of a Single Particle Model with Electrolyte," *J. Electrochem. Soc.* 166, A3693–A3706 (2019).
    doi:10.1149/2.0341915jes. arXiv:1905.12553. PyBaMM set "Marquis2019" (Kokam SLPB78205130H).
21. **[PYB] / [SULZER21]** V. Sulzer, S. G. Marquis, R. Timms, M. Robinson, S. J. Chapman, "Python Battery
    Mathematical Modelling (PyBaMM)," *J. Open Research Software* 9(1), 14 (2021). doi:10.5334/jors.309 (CC BY).
    Software: https://pypi.org/project/pybamm/, version 26.10.1.1, BSD-3-Clause licence (verified in the installed
    `LICENSE.txt`). Ea values cited from the parameter-set source files (Chen2020, Ecker2015, Marquis2019).
22. **[RUN]** This session: PyBaMM 26.10.1.1 DFN/SPM simulations. Scripts:
    `docs/research/scripts/cell_ecm_pybamm_{probe,ocv_fit,rc_fit2,capT}.py`. Outputs:
    `docs/research/data/cell_ecm_pybamm/*.json`.
23. **[BARC22]** S. Barcellona, S. Colnago, G. Dotelli, S. Latorrata, L. Piegari, "Aging effect on the variation of
    Li-ion battery resistance as function of temperature and state of charge," *J. Energy Storage* 50, 104658
    (2022). doi:10.1016/j.est.2022.104658. Post-print (open):
    https://re.public.polimi.it/bitstream/11311/1227399/5/ (read in full).
24. **[SUN23]** B. Sun, X. Qi, D. Song, H. Ruan, "Review of Low-Temperature Performance, Modeling and Heating for
    Lithium-Ion Batteries," *Energies* 16(20), 7142 (2023). doi:10.3390/en16207142. CC BY. Full text blocked
    (403); the 278.6% figure comes from a search excerpt.
25. **[BU802a]** Battery University, "BU-802a: How does Rising Internal Resistance affect Performance?" (table
    sourced to Siemens AG 2015).
    https://www.batteryuniversity.com/article/bu-802a-how-does-rising-internal-resistance-affect-performance/.
    Secondary; freely accessible.
26. **[ONSEMI]** onsemi product pages: ECH8693R "Dual N-Channel Power MOSFET for 1-2 Cells Lithium-ion Battery
    Protection 24V, 14A, 7mΩ"
    (https://www.onsemi.com/products/discrete-power-modules/mosfets/low-medium-voltage-mosfets/ech8693r) and
    EFC8811R "…1-Cell…, 12 V, 27 A, 3.2 mΩ"
    (https://www.onsemi.com/products/discrete-power-modules/mosfets/small-signal-mosfets/efc8811r). Freely
    accessible.
27. **[ZXJ03]** S. S. Zhang, K. Xu, T. R. Jow, "The low temperature performance of Li-ion batteries," *J. Power
    Sources* 115, 137–140 (2003). doi:10.1016/S0378-7753(02)00618-3. Abstract only.
28. **[JOW18]** T. R. Jow, S. A. Delp, J. L. Allen, J.-P. Jones, M. C. Smart, "Factors Limiting Li+ Charge Transfer
    Kinetics in Li-Ion Batteries," *J. Electrochem. Soc.* 165, A361–A367 (2018). doi:10.1149/2.1221802jes. Open
    access. Background only; no numbers used.
29. **[VISW10]** V. V. Viswanathan, D. Choi, D. Wang, W. Xu, S. Towne, R. E. Williford, et al., "Effect of entropy change of lithium intercalation in cathodes and anodes on Li-ion battery thermal
    management," *J. Power Sources* 195, 3720–3729 (2010). doi:10.1016/j.jpowsour.2009.11.103. OSTI record
    https://www.osti.gov/biblio/978956 (abstract-level).
30. **[BERN85]** D. Bernardi, E. Pawlikowski, J. Newman, "A General Energy Balance for Battery Systems,"
    *J. Electrochem. Soc.* 132, 5–12 (1985). doi:10.1149/1.2113792.
31. **[OVEJ19]** V. J. Ovejas, A. Cuadras, "Effects of cycling on lithium-ion battery hysteresis and
    overvoltage," *Scientific Reports* 9, 14875 (2019). doi:10.1038/s41598-019-51474-5. CC BY 4.0.
32. **[MM93]** J. F. Manwell, J. G. McGowan, "Lead acid battery storage model for hybrid energy systems," *Solar
    Energy* 50, 399–405 (1993). doi:10.1016/0038-092X(93)90060-2.
33. **[JH09]** M. R. Jongerden, B. R. Haverkort, "Which battery model to use?" *IET Software* 3, 445–457
    (2009). doi:10.1049/iet-sen.2009.0001. Author version: https://ris.utwente.nl/ws/files/6728809/IET-final-beta.pdf
    (read in full).
34. **[RV01]** D. Rakhmatov, S. Vrudhula, "An analytical high-level battery model for use in energy management of
    portable electronic systems," *Proc. IEEE/ACM ICCAD 2001*, 488–493. doi:10.1109/ICCAD.2001.968687. PDF:
    https://cecs.uci.edu/~papers/compendium94-03/papers/2001/iccad01/pdffiles/08c_2.pdf (read in full).
35. **[DS06]** D. Doerffel, S. A. Sharkh, "A critical review of using the Peukert equation for determining the
    remaining capacity of lead-acid and lithium-ion batteries," *J. Power Sources* 155, 395–400 (2006).
    doi:10.1016/j.jpowsour.2005.04.030. Abstract and secondary summaries only.
36. **[ECK15]** M. Ecker, S. Käbitz, I. Laresgoiti, D. U. Sauer, "Parameterization of a Physico-Chemical Model of a
    Lithium-Ion Battery: II. Model Validation," *J. Electrochem. Soc.* 162, A1849–A1857 (2015).
    doi:10.1149/2.0541509jes (Part I is in the same issue, not checked here). Used only via the PyBaMM
    "Ecker2015" Ea values (Kokam SLPB 75106100).
37. **[CHEN20]** C.-H. Chen, F. Brosa Planella, K. O'Regan, D. Gastol, W. D. Widanage, E. Kendrick, "Development
    of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,"
    *J. Electrochem. Soc.* 167, 080534 (2020). doi:10.1149/1945-7111/ab9050. PyBaMM "Chen2020" (LG M50,
    NMC811 / graphite-SiOx, 5 Ah).
38. **[OXF-BDD]** C. Birkl, D. Howey, "Oxford Battery Degradation Dataset 1" (2017), 8 × Kokam SLPB533459H4 740 mAh
    pouch cells at 40 °C. https://ora.ox.ac.uk/objects/uuid:03ba4b01-cfed-46d3-9b1a-7d4a7bdf6fac. ODbL 1.0 /
    DbCL 1.0. Cathode chemistry not verified here. A candidate small-pouch validation dataset.

---

## Appendix A1: Equilibrium OCV, Ai2020 LCO/graphite pouch, 25 °C ([RUN])

| z | 0.00 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| U (V) | 3.0000 | 3.5421 | 3.6396 | 3.6654 | 3.6954 | 3.7168 | 3.7331 | 3.7481 | 3.7622 | 3.7763 | 3.7925 |

| z | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 | 0.95 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|
| U (V) | 3.8122 | 3.8365 | 3.8675 | 3.9089 | 3.9519 | 3.9945 | 4.0409 | 4.0879 | 4.1384 | 4.2000 |

For comparison, the Marquis2019 (Dualfoil LCO) equilibrium OCV at z = 0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 is
3.105, 3.682, 3.697, 3.745, 3.799, 3.904 and 4.100 V. It is flatter in mid-SOC, which is less representative.

## Appendix A2: Cell entropic coefficient ∂U/∂T (mV/K), Ai2020 fits ([RUN])

| z | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 | 0.70 | 0.80 | 0.90 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ∂U/∂T | −0.706 | −0.534 | −0.479 | −0.372 | −0.307 | −0.248 | −0.248 | −0.205 | −0.153 | −0.097 | −0.046 |

## Appendix A3: 1RC fits to DFN pulse responses ([RUN]; 1C, 60 s pulse, 1 h rest; R in mΩ, τ in s)

| Set (Q) | T (°C) | z | R0 | R1 | τ1 | RMSE (mV) |
|---|---|---|---|---|---|---|
| Ai2020 (2.28 Ah) | 25 | 0.8 / 0.5 / 0.2 | 42.6 / 41.5 / 49.5 | 17.9 / 14.4 / 16.8 | 32 / 35 / 36 | ≤ 0.6 |
| Ai2020 | 0 | 0.8 / 0.5 / 0.2 | 50.4 / 49.0 / 57.1 | 44.5 / 40.0 / 45.0 | 89 / 99 / 96 | ≤ 1.0 |
| Marquis2019 (0.68 Ah) | 25 | 0.8 / 0.5 / 0.2 | 119.5 / 123.6 / 142.4 | 20.7 / 18.9 / 14.8 | 23 / 29 / 26 | ≤ 0.25 |
| Marquis2019 | 0 | 0.8 / 0.5 / 0.2 | 258.6 / 255.7 / 276.5 | 50.8 / 64.8 / 40.9 | 64 / 66 / 61 | ≤ 0.7 |
| Marquis2019 | −10 | 0.8 / 0.5 / 0.2 | 345.1 / 336.5 / 357.3 | 81.4 / 106.1 / 66.4 | 111 / 89 / 84 | ≤ 1.0 |
| Marquis2019 | −20 | 0.8 / 0.5 / 0.2 | 460.0 / 447.6 / 463.5 | 130.8 / 178.3 / 112.3 | 166 / 131 / 89 | ≤ 2.0 |

R0 at 0.1 s from separate 10 s pulses at 25 °C, z = 0.1: Ai2020 56.8 mΩ and Marquis2019 156.2 mΩ.

## Appendix A4: Recommended baseline for our ECM (synthesis; replace with fitted values when data allow)

| Quantity | Baseline | Sensitivity range |
|---|---|---|
| Q_n (new) | 4.7 Ah | 3.5–6.0 Ah |
| V_max / V_nom / V_sys cut-off | 4.45 / 3.87 / 3.40 V | 4.40–4.48 / 3.85–3.89 / 3.2–3.5 V |
| OCV | Chen–RM fit (A1) + stretch ΔV = 0.242 V, p = 3.86 | ± 20 mV offset (hysteresis and cell-to-cell) |
| R0,pack(25 °C, z ≥ 0.3) | 60 mΩ | 35–130 mΩ |
| R0(z) shape | 1 + 0.80·exp(−z/0.139) (normalised) | A 0.5–2.0 |
| R1, τ1 (25 °C) | 10 mΩ, 30 s | 3–25 mΩ, 10–60 s |
| Ea,R0 / Ea,R1 | 18 / 30 kJ/mol | 15–25 / 25–55 kJ/mol |
| ∂U/∂T | Appendix A2 (≈ −0.25 mV/K mid-SOC) | × (0.5–1.5) |
| η (charge) | 0.99 | 0.98–1.00 |
| Ageing | Q(N) to 80% at 500 or 1000 equivalent full cycles; R grows 1.3–2× at end of life | – |
