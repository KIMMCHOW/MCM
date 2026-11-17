# Proposal P1: Physics-first coupled electro-thermal-aging ECM (model specification)

2026 MCM Problem A, *Modeling Smartphone Battery Drain*. Design proposal written 2026-11-17.

**Status.** This is an implementation-ready specification. Every parameter has a value, a unit, a source key and an
uncertainty range. The source keys resolve to the research notes in `docs/research/` (summary list in §16).
The tags have fixed meanings:

* **[DERIVED]**: our own arithmetic on cited numbers; the arithmetic is shown.
* **[ASSUMPTION]**: a modelling choice with no measured source. It always carries a UQ range.
* **[unverified]**: background knowledge we could not check in this session. It is given as a range.
* **[PROTO]**: a number computed in this session by a scratch implementation of exactly the equations below
  (SciPy LSODA, rtol 1e-7, script in `docs/design/scripts/`). PROTO numbers are regression targets for
  `src/phonebatt`. They are **not** validation results.

**Design philosophy.** The phone draws **power**, not current, because the PMIC regulates every rail. The battery is
a Thevenin cell whose resistance depends on temperature and state of charge. That cell sits inside a lumped thermal
body that is heated by everything the phone does. Time-to-empty (TTE) is the first time the *voltage* fails. That
can happen at the fuel-gauge cut-off, in an instantaneous brown-out during a load burst, or by power
infeasibility. Running out of charge comes after all three. Battery aging enters as a slow state that shrinks
capacity and raises resistance. The model contains no curve fit to TTE data. Every term is a physical law or a
cited measurement, and the open datasets are held out for validation.

---

## 0. One-page summary

| Item | Choice |
|---|---|
| Fast states | `x = (z, v1, T)`: coulombic SOC, RC polarisation voltage, lumped device/battery temperature |
| Slow states | `(q_cal, q_cyc)`: calendar and cycle capacity loss. These give `q_loss = q_cal + q_cyc` and `r_growth = γ_R·q_loss`. They are frozen during one discharge |
| Load | battery-terminal power `P(t) = Σ_k P_k(u(t), T_s)`: sum of component sub-models (display, SoC with DVFS and thermal throttling, Wi-Fi/LTE/5G with RRC tails and signal penalty, GNSS, camera, audio, floor) |
| Current | `I = 2P / (Ũ + √(Ũ² − 4R0P))` with `Ũ = U(z,T) − v1`. This is the stable root of `R0 I² − Ũ I + P = 0` |
| Cell | Chen–Rincón-Mora OCV fitted to LCO/graphite (Ai2020) and stretched to 4.45 V; Arrhenius `R0(z,T)`, `R1(T)`, `τ1(T)`; entropic heat |
| Thermal | `C_th dT/dt = P + I²R0 + v1²/R1 − I T ∂U/∂T − (T − T_amb)/R_th`; skin hot spot `T_s`; sigmoid throttling of SoC power near 43 °C |
| TTE | first of: `V ≤ V_cut = 3.40 V` (fuel-gauge empty); `V_pk ≤ V_sys = 3.30 V` under a load burst `P + ΔP_b` (brown-out); `Ũ² < 4R0P` (infeasible power); `z ≤ 0` |
| Reference phone | **RP-24**: generic 2024 6.3-inch OLED flagship, 4.70 Ah / 18.2 Wh (usable OCV energy 17.79 Wh), pack R0 60 mΩ |
| Scenarios | 12 scenarios (S1–S12) built from 11 activity primitives, plus a stochastic CTMC day (S12-R) |
| Analytics | energy bucket; exact quasi-steady quadrature (within 0.4 % of the ODE [PROTO]); brown-out threshold formula; thermal-cap fixed point; CTMC renewal-CLT mean and variance of TTE |
| Validation | Bonn iPhone runs (blind, pre-registered predictions in §7), LG HG2 temperature data, NASA aging data, Oxford dynamic load, PyBaMM DFN, published-spec band |

**Headline [PROTO] behaviour of RP-24 (fresh, 25 °C, from 100 %).**

| Scenario | TTE |
|---|---|
| standby | 17.8 days |
| reading | 35.1 h |
| social over Wi-Fi | 16.2 h |
| video over Wi-Fi / LTE / 5G | 29.7 / 23.2 / 20.4 h |
| navigation | 11.2 h |
| gaming (thermally capped) | 5.6 h |
| video call over LTE | 6.8 h |
| typical day | 34.7 h |

The same social load on an aged pack (SOH 0.8) at −20 °C browns out while still showing **23.6 %**. Cold, aging and
peak power interact much more than additively.

---

## 1. Notation, states, inputs, outputs

### 1.1 Units and conventions

Internally the model uses SI units, with three exceptions: `T` is in kelvin (°C only in tables), capacity `Q` is in
A·h, and the factor 3600 appears explicitly. `I > 0` means discharge, and every power is measured at the **battery
terminals**. Constants: `R_g = 8.314 J mol⁻¹ K⁻¹`, `T_ref = 298.15 K`. The Arrhenius factor is

```
A(T; E) = exp[ (E/R_g)·(1/T − 1/T_ref) ]        (> 1 when T < T_ref)
```

### 1.2 Time scales

| Process | Time scale | Source |
|---|---|---|
| Load bursts (CPU boost, radio promotion, camera start) | ms to 1 s | HU12, LEE20 |
| RC polarisation `τ1` | 18 s (45 °C), 20 s (40 °C), to 135 s (−20 °C) | RUN (cell_ecm) |
| Thermal `τ_th = R_th·C_th` | 7.5 min (5–15) | HE19 |
| Usage switching | minutes | ASSUMPTION |
| Discharge (TTE) | 3 h to 18 days | this model |
| Aging (SOH) | weeks to years | SCH14, AP-BATT |

The spread of time scales fixes the design:

* Bursts are handled **algebraically**: they see `R0` only, while `v1` stays frozen.
* RC and thermal dynamics are **ODE states**.
* SOH is a **slow state**. It is held constant within one TTE run and integrated in the aging layer (§10).

### 1.3 State vector

| State | Symbol | Unit | Range | Initial condition |
|---|---|---|---|---|
| Coulombic SOC, relative to the current maximum capacity `Q_max` | `z` | – | [0, 1] | `z0 = z_e + s0(1 − z_e)` from the displayed SOC `s0` (§2.11) |
| RC polarisation voltage | `v1` | V | ≥ 0 on discharge | 0 (rested) |
| Lumped device and battery temperature | `T` | K | – | `T0 = T_amb` (soaked), unless the scenario sets it |
| Calendar capacity loss (slow) | `q_cal` | – | [0, 0.3] | 1e-4 (new) or the scenario value |
| Cycle capacity loss (slow) | `q_cyc` | – | [0, 0.3] | 1e-4 (new) or the scenario value |

Derived quantities: `q_loss = q_cal + q_cyc`, `SOH = 1 − q_loss`, `r_growth = γ_R·q_loss` and
`Q_max = Q_n(1 − q_loss)`.

### 1.4 Usage input vector u(t) and environment e(t)

`u(t)` is piecewise constant. It switches between activity primitives (§5.1), which are either scheduled or driven
by a CTMC.

| Field | Symbol | Type / range | Meaning |
|---|---|---|---|
| Screen on | `s_scr` | {0, 1} | display powered |
| Awake | `a` | {0, 1} | application processor held awake. Defined as `a = max(s_scr, audio, g>0, c_r, c_f, κ=rt)` |
| Brightness | `b` | [0, 1] | OS brightness-slider fraction |
| Average picture level | `APL` | [0, 1] | linear-light mean pixel level. Light UI ≈ 0.7; dark UI ≈ 0.15 |
| Effective refresh rate | `f_r` | Hz, {60, 90, 120} | rate the panel actually runs at (LTPO static content = 60) |
| SoC demanded throughput | `x_c` | [0, 1] | fraction of maximum CPU+GPU throughput the workload asks for |
| Wi-Fi associated | `w` | {0, 1} | |
| Cellular camping RAT | `c` | {off, LTE, NR} | NR = 5G NSA low-band |
| Data path | `ρ` | {wifi, cell, none} | RAT carrying the traffic |
| Traffic class | `κ` | {burst, rt} | `burst`: Poisson bursts with tail timers (browse, stream, sync). `rt`: continuous real-time media (voice or video call) |
| Burst rate | `λ` | s⁻¹ | arrival rate of traffic bursts or chunks |
| Mean DL / UL rate | `d_dl`, `d_ul` | Mbit/s | time-averaged throughput |
| Signal strength | `RSRP` | dBm | cellular reference signal received power |
| GNSS duty | `g` | [0, 1] | fraction of time the GNSS engine tracks |
| Rear / front camera | `c_r`, `c_f` | {0, 1} | camera pipeline active |
| Audio playback / Bluetooth audio | `au`, `bt` | {0, 1} | |
| Sensor power | `P_sens` | W | sum of active sensors |
| Burst amplitude | `ΔP_b` | W | characteristic short-peak power above the mean in this activity (§3.8) |
| Ambient temperature (environment) | `T_amb(t)` | K | |
| Enclosure factor (environment) | `k_env` | – | multiplies `R_th` (case, hand, pocket) |

### 1.5 Outputs

* `z(t)`, and the displayed SOC `s(t)` (§2.11).
* Terminal voltage `V(t)`, current `I(t)` and power `P(t)` with its component breakdown `P_k(t)`.
* Temperatures `T(t)` and `T_s(t)`, and the throttle factor `σ(t)`.
* `TTE`, the event type (`cutoff` / `brownout` / `infeasible` / `charge`), and the **stranded charge**
  `s_strand = s(TTE)`: the charge displayed at shutdown, which comes back once the phone warms (LEE20 C12).
* For UQ: the TTE distribution, prediction intervals, and `Pr(TTE < t_target)`.

---

## 2. Governing equations

### 2.1 Power demand at the battery terminals

```
P(t) = P_floor + P_disp + P_SoC(T_s) + P_net + P_GNSS + P_cam + P_aud + P_sens                         (E1)
P_floor = P_susp + a·ΔP_awake
```

* **Battery-side values.** All component parameters are battery-side values, so they already contain the PMIC loss.
  `η_PMIC = 0.83` (CH13) is used only to convert rail-level sources (LX CPU coefficients, GU24 ODPM rails) to battery
  side. It is never applied to battery-side measurements a second time (device_power §2.1).
* **Additivity.** This is assumption A10, and the Bonn data support it: `1/TTE` for camera + GPS is predicted from
  the single-sensor runs to within 0.15 % (validation_data a.2).

### 2.2 Current under constant power, and feasibility

Let `Ũ = U(z,T) − v1`. Power balance `P = V·I` and Ohm's law `V = Ũ − I·R0` give `R0 I² − Ũ I + P = 0`. The
physical (high-voltage) root, written in a cancellation-free form, is

```
D = Ũ² − 4·R0·P                                                                                  (E2)
I = 2P / (Ũ + √D),      V = (Ũ + √D)/2        (requires D ≥ 0)
```

* **Power ceiling.** `P_max = Ũ²/(4R0)` is the largest power the cell can deliver, at `V = Ũ/2`. A `D < 0` means
  voltage collapse.
* **The cut-off binds first.** RP-24 always has `V_cut = 3.40 V > Ũ/2 ≈ 1.8–2.2 V`, so the cut-off event is reached
  before `D < 0` whenever the demand changes smoothly.
* **The binding limit.** The physically binding limit is therefore the cut-off-limited power
  `P_cut(V*) = V*(Ũ − V*)/R0` (thermal_aging T3).
* **When D < 0 is reachable.** It can be reached only if a demand step lands beyond `P_max`. The RHS then uses the
  clamped current `I = Ũ/(2R0)`, and the step is flagged as event E4.

### 2.3 Charge conservation

```
dz/dt = − I / (3600 · Q_max),        Q_max = Q_n·(1 − q_loss)                                       (E3)
```

Coulombic efficiency is 1 on discharge. It is 0.99 on charge, used only in the aging layer (PLETT-N).

### 2.4 Polarisation (1RC Thevenin)

```
dv1/dt = (R1(T)·I − v1) / τ1(T),        C1 = τ1/R1                                                   (E4)
```

The 2RC form (τ ≈ 7 s and 47 s, RUN) is the structural variant M4.

### 2.5 Open-circuit voltage U(z, T)

```
U(z,T) = −b0·exp(−b1 z) + a0 + a1 z + a2 z² + a3 z³  +  ΔV_hv·z^p  +  δU  +  (T − T_ref)·∂U/∂T(z)     (E5)
```

* **Fit.** Chen–Rincón-Mora form (CRM06) fitted to the Ai2020 LCO/graphite equilibrium OCV (RUN, RMSE 7 mV):
  b0 = 0.60117, b1 = 36.426, a0 = 3.6193, a1 = 0.44311, a2 = −0.48605, a3 = 0.63127.
* **High-voltage stretch.** `ΔV_hv = 0.242 V`, `p = 3.86` [ASSUMPTION; RUN]. These give `U(1) = 4.450 V` and a mean
  OCV of 3.870 V. The stretch leaves `z < 0.3` almost unchanged.
* **Offset.** `δU` is a ±20 mV offset that stands in for hysteresis and cell-to-cell spread (OVEJ19 mechanism;
  magnitude unverified). It is 0 at baseline.
* **Entropic coefficient.** `∂U/∂T(z)` is linear interpolation in the table of cell_ecm Appendix A2 (mV/K):

| z | 0 | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.50 | 0.60 | 0.70 | 0.80 | 0.90 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ∂U/∂T (mV/K) | −0.449 | −0.706 | −0.534 | −0.479 | −0.372 | −0.307 | −0.248 | −0.248 | −0.205 | −0.153 | −0.097 | −0.046 |

Values on the RP-24 OCV [PROTO]:

| z | 0.10 | 0.20 | 0.30 | 0.50 | 0.80 | 1.00 |
|---|---|---|---|---|---|---|
| U (V) | 3.644 | 3.694 | 3.728 | 3.815 | 4.088 | 4.450 |

### 2.6 Resistances

```
R0(z,T) = R0∞·[1 + A_R·exp(−z/z0)] · A(T; E_a,R0) · (1 + r_growth)                                 (E6)
R0∞ = R0,ref / (1 + A_R·exp(−0.5/z0))       so that R0(0.5, T_ref) = R0,ref = 60 mΩ
R1(T) = R1,ref·A(T; E_a,R1)·(1 + r_growth)
τ1(T) = τ1,ref·A(T; E_a,τ)
```

Parameter values: `A_R = 0.80`, `z0 = 0.139`, `E_a,R0 = 18 kJ/mol`, `R1,ref = 10 mΩ`, `E_a,R1 = 30 kJ/mol`,
`τ1,ref = 30 s`, `E_a,τ = 21 kJ/mol` (RUN, BARC22, SUN23).

RP-24 values [PROTO]:

| T (°C) | 40 | 25 | 10 | 0 | −10 | −20 |
|---|---|---|---|---|---|---|
| R0 at z = 0.5 (mΩ) | 42.4 | 60.0 | 88.1 | 116.6 | 157.6 | 218.1 |
| R1 (mΩ) | 5.6 | 10.0 | 19.0 | 30.3 | 50.0 | 86.0 |
| τ1 (s) | 20.0 | 30.0 | 47.0 | 65.1 | 92.6 | 135.3 |

At 25 °C, `R0(z)` is 81.6 mΩ at z = 0.1 and 69.9 mΩ at z = 0.2.

* **Consistency check against phone-level data.** The total `(R0 + R1)(−5 °C)/(R0 + R1)(25 °C)` is 2.5. LEE20 report
  2.0–3.3 on four phones.
* **Why a single Arrhenius term.** `R0` includes the protection board and contacts (PINE: PCM 20–60 mΩ), whose
  resistance hardly varies with temperature. Putting one Arrhenius factor on the whole of `R0` matches the
  phone-level ratios better. The split alternative is structural variant M5 (§9.4).

### 2.7 Lumped thermal balance (Bernardi heat, Newton cooling)

```
C_th·dT/dt = P + I²·R0 + v1²/R1 − I·T·∂U/∂T(z) − (T − T_amb)/(k_env·R_th)                         (E7)
T_s = T_amb + (R_skin/R_th)·(T − T_amb)                    (quasi-steady skin hot spot)
```

* **`P` (= I·V)**: everything delivered to the phone becomes heat inside it. Emitted light and RF are neglected
  (thermal_aging §2.1).
* **`I²R0 + v1²/R1`**: irreversible battery heat. The capacitor term `C1·v1·dv1/dt` is stored, not dissipated
  (BERN85).
* **`−I·T·∂U/∂T`**: reversible entropic heat. It is positive on discharge, and comparable to the Joule heat at ≈ 1 A
  for LCO (cell_ecm finding 8).
* **Parameters**: `C_th = 150 J/K`, `R_th = 3.0 K/W` (`τ_th = 7.5 min`) and `R_skin = 5.5 K/W` (HE19, KANG19, PU23,
  MAL99).

### 2.8 Thermal throttling (feedback from T_s to load)

The SoC demand (§3.3) is `P_SoC,dem`. Above the floor `P_fl` it is scaled continuously:

```
P_SoC = P_SoC,dem                                              if P_SoC,dem ≤ P_fl
P_SoC = P_fl + (P_SoC,dem − P_fl)·σ((T_thr − T_s)/w_thr)       otherwise,   σ(y) = 1/(1 + e^(−y))         (E8)
```

* **Parameters**: `T_thr = T_lim − 2 K = 43 °C` (`T_lim = 45 °C`, CHI15), `w_thr = 1 K` and `P_fl = 0.5 W`
  [ASSUMPTION].
* **The chain is explicit.** `T → T_s → P → I` has no algebraic loop.
* **Equilibrium.** The steady state is the fixed point `T_s* = T_amb + R_skin·(P(T_s*) + q_b)` (§6.3, F5). With these
  parameters the capped total power is ≈ 3.1–3.3 W at 25 °C, inside the measured 2.2–3.8 W (CHI15).

### 2.9 Load bursts and brown-out

Within each activity, short bursts (≪ τ1) raise the power to `P_pk = P + ΔP_b`. During a burst `v1` stays frozen, so

```
D_pk = Ũ² − 4R0·P_pk,      V_pk = (Ũ + √D_pk)/2       (V_pk := 0 if D_pk < 0)                     (E9)
brown-out  ⇔  V_pk ≤ V_sys  ⇔  P + ΔP_b ≥ P_cut,sys = V_sys·(Ũ − V_sys)/R0
```

* **`V_sys` vs `V_cut`.** `V_sys = 3.30 V` is the "minimum allowed system voltage" (QC-SMB). It differs from the
  fuel-gauge cut-off `V_cut = 3.40 V` (QC-FG; LEE20 "≈ 3.4 V").
* **The mechanism.** E9 is what produces shutdowns at 15–35 % displayed SOC in the cold and on aged packs (LEE20,
  HE19, AP-BATT). A mean-power model cannot produce them.

### 2.10 Aging coupling

Within a TTE run, `q_loss` and `r_growth` are constants. They scale `Q_max` (E3) and multiply `R0` and `R1` (E6). The
aging ODEs that evolve them over days to years are in §10.

### 2.11 Displayed SOC (output map) and stranded charge

`z_e` is the coulombic SOC at which the rested OCV equals the cut-off, `U(z_e, T_ref) = V_cut`. For RP-24,
`z_e = 0.0263` [PROTO]. Two candidate display maps are carried:

```
D1 (coulombic, primary):  s = (z − z_e)/(1 − z_e)
D2 (state of energy):     ς = ∫_{z_e}^{z} U(ζ,T_ref) dζ / ∫_{z_e}^{1} U(ζ,T_ref) dζ
```

* **Usage.** D1 sets initial conditions and is the reported "SOC(t)". D2 is tested against the iOS traces (§7, V1.4).
* **iOS top-hold map.** Used for Bonn comparisons only (validation_data b.4):
  `s_iOS = min(1, s/(1 − h))`, with `h = 0.065`, floored to 1 % steps.
* **Stranded charge.** `s_strand = s(t_TTE)`.

### 2.12 Complete fast system

For `x = (z, v1, T)` and given `(u(t), e(t), θ, q_loss)`:

```
dz/dt  = − I(x,u)/(3600 Q_max)
dv1/dt = (R1(T) I(x,u) − v1)/τ1(T)
dT/dt  = [P(u,T_s) + I²R0 + v1²/R1 − I T ∂U/∂T − (T − T_amb)/(k_env R_th)] / C_th
with I(x,u) from (E2) using P(u, T_s(T)).
```

* **Well-posedness.** On the feasible set `F = {D > 0, 0 ≤ z ≤ 1}` the RHS is locally Lipschitz in `x`: the sigmoid,
  the exponentials and the square root of a positive argument are smooth, and the `∂U/∂T` interpolant is piecewise
  linear. It is piecewise continuous in `t`, because `u` is piecewise constant. A unique solution therefore exists
  up to the first exit from `F` or the first event.
* **Monotonicity.** `z` is strictly decreasing whenever `P > 0`.
* **Bounded temperature.** `T` stays within `T_amb + k_env R_th·(P_max + q_b,max)`.

---

## 3. Component power sub-models (battery-side)

Source keys resolve in §16. Each range listed here is the UQ range of §4.

### 3.1 Floor

| Term | Form | Default | Range | Source |
|---|---|---|---|---|
| Suspend floor: screen off, airplane-equivalent (SoC suspend, PMIC, RAM self-refresh) | constant | 20 mW | 10–30 | CH13 24–27 mW on 2008–13 phones; modern value unverified |
| Awake overhead: AP awake, no display | `a·ΔP_awake` | 100 mW | 50–250 (log-uniform) | CH10 Freerunner 269 vs 69 mW; AOSP cpu.awake 50 mA (illustrative); L |

### 3.2 Display (OLED)

```
P_disp = s_scr·[ P_d0 + P_c,ref·(b/0.5)^γ_b·(APL/APL_ref) + κ_RR·max(f_r − 60, 0) ]                    (E10)
```

**Physics.** OLED pixel power is linear in linear-light RGB (DZ11, R² > 0.95), so the content term scales with APL.
Black pixels cost only the fixed overhead `P_d0` (CH10, CH13).

**Parameters.**

| Symbol | Value | Range | Source |
|---|---|---|---|
| `P_d0` | 0.15 W | 0.10–0.25 | CH13 223 mW constant term (rail, S3); AOSP screen.on |
| `P_c,ref` (content power at `b = 0.5`, `APL_ref = 0.7`) | 0.25 W | 0.12–0.50 | GU24 total display ≈ 0.30 W (Pixel 7a, 6.1") / 0.45 W (S23+, 6.6"); CH13 300–600 mW at 50 % (4.8") |
| `APL_ref` | 0.7 | – | DZ11: 65 % white pixels in mobile web |
| `γ_b` | 2.0 | 1.0–2.9 | see below |
| `κ_RR` | 2.8 mW/Hz | 1.7–3.3 | AT20 [DERIVED]: 0.16–0.18 W from 60 to 120 Hz on a fixed-rate panel (M/L) |

**Brightness exponent `γ_b`.**

* `γ_b = 1` reproduces CH13's near-affine 2012 slider: max/50 % ratio 2.0.
* `γ_b ≈ 2.9` is [DERIVED] from DH21. Dark mode saves 3–9 % at 30–50 % brightness but 39–47 % at 100 %, which
  requires the content power at 100 % to be ≈ 14× that at 40 %. Hence `γ = ln 14/ln 2.5`.
* With `γ_b = 2` and dark-UI `APL = 0.15`, the model gives a dark-mode saving of ≈ 6 % (b = 0.3) and ≈ 40 %
  (b = 1.0) of total power on the S4 social load [DERIVED], so it fits DH21 at both ends. It also gives full-white maximum-brightness display power `0.15 + 0.25·4/0.7 = 1.58 W`, inside
  the 1.0–2.0 W band of device_power §0.
* At `b = 0.4` (typical indoor auto-brightness, DH21 D11), `P_disp = 0.31 W`. This matches GU24.

**LCD variant (generalisation).** `P_LCD = s_scr·[P_panel + c1 b + c2 b²]`, content-independent and convex (CH10,
HE20).

### 3.3 SoC (CPU + GPU + memory) with DVFS

The Linux energy model gives `P = C·V²·f` (LX). With a utilisation-targeting governor, `f ∝ x_c`, and the OPP
voltage is affine, `V(f) = V0 + s·f`. For the MT8183 A73 OPPs:

* `s = 0.293 V/GHz`, `V0 = 0.468 V` at the zero-frequency intercept, `V_max = 1.05 V` [DERIVED].
* `v0 = V0/V_max = 0.446`.

```
P_SoC,dem(x_c) = P_SoC,max · x_c · [v0 + (1 − v0)·x_c]²                                               (E11)
```

* **Convexity.** The effective exponent is ≈ 1.9 near the top (device_power finding 4). The form makes
  "race-to-idle" versus "slow-and-steady" explicit.
* **Ceiling.** `P_SoC,max = 5.0 W` (3–7): CH13 CPU 2.85 W + GPU 1.42 W rails ×1.2.
* **Check values.**

| `x_c` | 0.05 | 0.10 | 0.25 | 0.30 | 0.80 |
|---|---|---|---|---|---|
| `P_SoC,dem` (W) | 0.056 | 0.126 | 0.427 | 0.562 | 3.16 |

`x_c = 0.25` matches GU24's 0.4–0.5 W of compute in Instagram. `x_c = 0.8` matches the ≈ 3 W sustainable SoC power.

* **Throttling.** E8 applies to this term. Leakage growth with temperature is left out of the baseline; it is
  structural variant M12.

### 3.4 Radios: Wi-Fi, LTE, 5G NSA (RRC tails, Poisson bursts, signal penalty)

```
P_net = w·P_idle^W + 1[c ≠ off]·P_idle^c + Φ_ρ                                                      (E12)
Φ_ρ (κ = burst) = m_ρ·[ (P_tail^ρ − P_idle^ρ)·(1 − e^(−λ T_tail^ρ)) + λ·e^(−λ T_tail^ρ)·E_pro^ρ + e_b^ρ·(d_dl + κ_ul·d_ul) ]
Φ_ρ (κ = rt)    = m_ρ·P_rt^ρ
m_ρ = 1 for Wi-Fi;   m_ρ = m_sig(RSRP) = min(m_max, 1 + γ_sig·max(0, (RSRP_ref − RSRP)/10 dB)) for cellular
```

**Derivation of the burst class** [our construction on HU12/NA21 state machines].

* Bursts arrive as a Poisson process of rate λ. Each burst resets an inactivity timer `T_tail`, and bursts are short
  compared with `T_tail`.
* The radio is in the tail state at time t if and only if a burst arrived in `(t − T_tail, t]`. This has probability
  `1 − e^(−λT_tail)`, exactly.
* A promotion happens only when a burst arrives after the timer has expired. That occurs at rate `λ·e^(−λT_tail)`.
* The marginal energy per bit `e_b` covers the extra active-state time during transfers.

**Real-time class.** `rt` covers voice and video calls. It keeps the modem in continuous Tx/Rx with no C-DRX sleep,
so the radio power is the measured connected-active level.

**Parameters (battery-side, per RAT).**

| Parameter | Wi-Fi | LTE | 5G NSA low-band | Sources |
|---|---|---|---|---|
| `P_idle` (beacon/paging) | 4 mW (2–50) | 12 mW (5–35) | 20 mW (10–50) | CH15 1.1 mA Wi-Fi; HU12 19.7 mW LTE DRX [DERIVED]; CH15 2.3–8.3 mA; GU24 5G ≈ 2× 4G standby (L) |
| `P_tail` | 0.12 W | 0.15 W (0.066–0.25) | 0.25 W (0.25–0.60) | HU12 (Wi-Fi 119 mW); NA21 Table 2 (4G 66–178, NSA-LB 249–260, SA 593 mW) |
| `T_tail` | 0.24 s | 10 s (5–11.6) | 10.3 s (10.2–21.4) | HU12; NA21 Table 7; XU20 chained 21.4 s |
| `E_pro` | 0.010 J | 0.31 J (0.1–0.5) | 1.0 J (0.5–2.0) | HU12 (124 mW·79 ms; 1211 mW·260 ms); NA21 4G→5G switch 0.70–0.80 W × 1.44 s promotion [DERIVED] |
| `e_b` (DL, W per Mbit/s) | 0.015 (0.005–0.06) | 0.030 (0.014–0.08) | 0.030 (0.014–0.08) | [DERIVED]: SU14 14–60 nJ/bit; Apple streamed − local video = 0.06–0.08 W (APL14P); NA21 DL slope ≈ 14 mW/Mbps; GU24 4G − Wi-Fi YouTube 0.12–0.24 W |
| `κ_ul` (UL/DL energy ratio) | 5 (2–8.4) | 5 (2–8.4) | 5 (2–8.4) | NA21 68/14 ≈ 4.9; HU12 438/52 = 8.4 |
| `P_rt` (continuous real-time) | 0.30 W (0.13–0.8) | 1.2 W (0.6–2.6) | 1.4 W (1.2–3.0) | device_power §0 Wi-Fi active 0.3 W; HU12 LTE β = 1.29 W; GU24 Meet cellular rail 2.58 W (upper); GU24 5G/4G 1.12–1.17; XU20 2–3× on a 2019 modem |
| `γ_sig`, `RSRP_ref`, `m_max` | – | 0.25 per 10 dB (0.10–0.60), −85 dBm, 3.5 | same | DI13: +52 % energy from −85 to −105 dBm; Tx/Rx/tail currents ×1.2–3.5 |

**Wi-Fi weak-signal cliff** (optional, scenario S11b):

```
m_W(RSSI) = 1 + 1.13·exp((−85 − RSSI)/2.54 dB)
```

[DERIVED] from DI13: +113 % at −85 dBm and +810 % at −90 dBm relative to −50 dBm.

**Check values (total `P_net`, idle terms of the associated radios included)** [PROTO].

| Case | Radio power |
|---|---|
| Video stream, 3 Mbit/s, λ = 0.25 s⁻¹ | Wi-Fi 0.070 W (LTE camping); LTE 0.235 W; 5G 0.342 W |
| Standby, 12 syncs/h over LTE (Wi-Fi associated) | 0.022 W |
| Standby, 60 syncs/h over 5G (Wi-Fi off) | 0.070 W |

* The Wi-Fi traffic term alone, `Φ`, is 0.054 W. Apple's streamed − local difference is 0.06–0.08 W.
* The LTE − Wi-Fi gap is 0.165 W. GU24 measured 0.12–0.24 W.

### 3.5 GNSS

`P_GNSS = g·P_G`, with `P_G = 0.25 W` (0.14–0.45).

* **Sources.** KW19 0.23 / 0.32 J per 1 Hz fix (single / dual frequency); CH13 0.39–0.43 W including the app; CH10
  0.14–0.17 W for the module only.
* **Cold fix.** A cold time-to-first-fix costs a one-off ≈ 154 J (CH13). It is optional, added as an energy impulse at
  GNSS start in S6/S9.
* **Known risk.** Bonn shows only 0.03 W on the iPhone 14 Pro (validation_data a.2). We do **not** tune to this value
  (see §7).

### 3.6 Camera

* **Rear pipeline (sensor, ISP, encoder; increment above the awake + display floor).** `c_r·ΔP_cam` with
  `ΔP_cam = 1.4 W` (1.0–2.0).
  * CH13 S3 still/video: 2.26–2.61 W whole device minus 0.805 W screen-on idle gives 1.45–1.81 W.
  * CX13: 1.40–1.65 W whole device while recording, including 0.17–0.18 W of AMOLED.
* **Front camera in a call.** `c_f·ΔP_camf` with `ΔP_camf = 0.30 W` (0.2–0.5): GU24 camera rail 0.26 W × 1.2.

### 3.7 Audio, Bluetooth, sensors, haptics

| Term | Value | Range | Source |
|---|---|---|---|
| Audio DSP/codec, headphones (`au·ΔP_aud`) | 0.05 W | 0.03–0.15 | CH13: 226 mW total, max volume adds only +5 mW; HE20 0.11 W incl. speaker |
| Bluetooth classic audio (`bt·ΔP_BT`) | 0.04 W | 0.03–0.05 | CH10 36–45 mW |
| Accelerometer / gyroscope / magnetometer | 5 / 30 / 12 mW | – | CH13 at 50 Hz |
| Haptics | omitted | – | no source; duty cycle tiny (device_power HP, unverified 50–200 mW while active) |

### 3.8 Burst amplitude ΔP_b [DERIVED]

| Activity class | ΔP_b | Basis |
|---|---|---|
| Standby wakeup | 1.0 W | LTE promotion 1.21 W (HU12); 4G→5G switch 0.70–1.49 W (NA21) |
| Interactive UI, navigation, call | 3.0 W | LEE20 peak currents 1–2 A ≈ 3.7–7.4 W total; CH13 CPU + GPU peak rails 4.3 W |
| Gaming, camera start | 4.0 W | CH13 device peak > 7 W |
| Video playback | 2.0 W | ASSUMPTION (decoder and network chunk peaks) |

In UQ every `ΔP_b` is multiplied by `U(0.5, 2.0)`. This parameter dominates the cold and aged brown-out (§9). It is
also the lever that an OS "peak-power governor" can pull (§11.2).

### 3.9 Deliberately omitted (with justification)

| Effect | Size | Why omitted |
|---|---|---|
| Rate-capacity (Peukert, KiBaM) | 0.2–1 % at 0.2–0.5C (RUN) | Phone C-rates are 0.01–0.6C. KiBaM is kept as structural variant M8 |
| Self-discharge | ≈ 0.1 %/month, ≤ 2 %/month with electronics (ZIL19, BU802b) | 2–3 orders below standby drain |
| OCV hysteresis | ≤ ~20 mV (unverified) | Absorbed in `δU` |
| PMIC light-load efficiency drop | unverified (device_power PM2) | Absorbed in the `P_susp` range |

---

## 4. Reference phone RP-24: full parameter table

**RP-24** is a generic 2024 Android flagship:

* 6.3-inch 19.5:9 OLED (97 cm² active area), LTPO 1–120 Hz;
* single LCO/graphite pouch: 4.70 Ah typical (the Pixel 9 spec-sheet value, PIX9), V_max 4.45 V, V_nom 3.87 V,
  18.2 Wh nominal, usable OCV energy 17.79 Wh at 25 °C [PROTO];
* mass 190 g (L, thermal_aging T20).

It is **not** a model of any specific product.

**Distributions.** U = uniform on the range; T = triangular (min, mode = default, max); LU = log-uniform. Log-uniform
is used where the range spans ≥ 3×.

### 4.1 Cell and pack

| Symbol | Value | Unit | Distribution / range | Source |
|---|---|---|---|---|
| `Q_n` | 4.70 | Ah | T(4.56, 4.70, 4.80) | PIX9: 4.70 typical / 4.558 minimum |
| OCV coefficients `b0, b1, a0..a3` | 0.60117, 36.426, 3.6193, 0.44311, −0.48605, 0.63127 | V, – | fixed | RUN on PYB-AI20 (CRM06 form) |
| `ΔV_hv`, `p` | 0.242, 3.86 | V, – | fixed (variant pairs 0.192/5.44 and 0.272/2.90) | RUN, ASSUMPTION |
| `δU` | 0 | mV | U(−20, +20) | OVEJ19 (mechanism), L |
| `∂U/∂T(z)` | table §2.5 | mV/K | multiplier U(0.5, 1.5) | RUN (Ai2020) |
| `R0,ref` (z = 0.5, 25 °C, pack incl. PCM) | 60 | mΩ | LU(35, 130) | PINE (≤ 100 mΩ at 1 kHz; PCM 20–60), BU802a (≈ 130 legacy), RUN (cell 19 mΩ) |
| `A_R`, `z0` | 0.80, 0.139 | –, – | A_R U(0.5, 2.0) | RUN |
| `E_a,R0` | 18 | kJ/mol | T(12, 18, 27) | RUN 18.0; BARC22 15.6; SUN23 15.8; LEE20-derived 12–27; KOL20 apparent 26.4 |
| `R1,ref` | 10 | mΩ | LU(3, 25) | RUN |
| `τ1,ref` | 30 | s | U(10, 60) | RUN 29–35 |
| `E_a,R1` | 30 | kJ/mol | T(25, 30, 55) | RUN 27.7–31.4; BARC22 29 |
| `E_a,τ` | 21 | kJ/mol | U(15, 27) | RUN |
| `V_cut` (fuel gauge, mean load) | 3.40 | V | U(3.30, 3.50) | QC-FG, LEE20 |
| `V_sys = V_cut − Δ_sys` (instantaneous) | 3.30 | V | `Δ_sys` U(0, 0.20) | QC-SMB |
| Usable-capacity spread (device to device) | – | – | via `Q_n` | Bonn CV 0.4–4.9 % (validation only) |

### 4.2 Thermal

| Symbol | Value | Unit | Distribution / range | Source |
|---|---|---|---|---|
| `C_th` | 150 | J/K | T(100, 150, 230) | MAL99 mass × c_p; HE19 τ·R; KANG19 energy/ΔT [DERIVED] |
| `R_th` | 3.0 | K/W | T(2, 3, 4) | HE19 Fig. 10 [DERIVED]; convection + radiation estimate |
| `R_skin` | 5.5 | K/W | T(4, 5.5, 9) | PU23 5 K/W (Pixel 5); KANG19 5–7; CHI15 |
| `T_lim` (so `T_thr = T_lim − 2 K`) | 45 | °C | U(41, 48) | CHI15; KANG19 pain at 42–45 °C; UL-TS 48 °C |
| `w_thr` | 1 | K | U(0.5, 3) | ASSUMPTION (graded AOSP-TH levels) |
| `P_fl` | 0.5 | W | U(0.3, 1.0) | ASSUMPTION |
| `k_env` | 1.0 | – | U(0.8, 1.5) | ASSUMPTION (KANG19: a case slowed heating by 54 %) |
| Thermal shutdown `T_shut` | 68 | °C | fixed | KANG19 Android example |

### 4.3 Device power

| Symbol | Value | Unit | Distribution / range | Source |
|---|---|---|---|---|
| `P_susp` | 0.020 | W | U(0.010, 0.030) | CH13 |
| `ΔP_awake` | 0.10 | W | LU(0.05, 0.25) | CH10, AOSP |
| `P_d0` | 0.15 | W | U(0.10, 0.25) | CH13, AOSP |
| `P_c,ref` | 0.25 | W | U(0.12, 0.50) | GU24, CH13 |
| `γ_b` | 2.0 | – | U(1.0, 2.9) | CH13, DH21 [DERIVED] |
| `κ_RR` | 2.8 | mW/Hz | U(1.7, 3.3) | AT20 [DERIVED] |
| `P_SoC,max` | 5.0 | W | U(3, 7) | CH13 |
| `v0` | 0.446 | – | U(0.35, 0.55) | LX [DERIVED] |
| `P_G` (GNSS) | 0.25 | W | U(0.14, 0.45) | KW19, CH13, CH10 |
| `ΔP_cam` / `ΔP_camf` | 1.4 / 0.30 | W | U(1.0, 2.0) / U(0.2, 0.5) | CH13, CX13 / GU24 |
| `ΔP_aud` / `ΔP_BT` | 0.05 / 0.04 | W | U(0.03, 0.15) / U(0.03, 0.05) | CH13, HE20 / CH10 |
| `ΔP_b` multiplier | 1 | – | U(0.5, 2.0) | §3.8 |

### 4.4 Radio

All radio parameters are in the §3.4 table, with the same U ranges. `γ_sig` is U(0.10, 0.60).

### 4.5 Aging

| Symbol | Value | Unit | Distribution / range | Source |
|---|---|---|---|---|
| `k_c` | 1.41e-3 | day^−½ | U(1.1e-3, 1.7e-3) | calibrated so that `q_cal(365 d, 25 °C, 50 % SOC) = 2.7 %`; SCH14 2.4 %, KEIL16 2–3 % after 10 months |
| `E_cal` | 58 | kJ/mol | U(20, 60) | SCH14 via CA21 (6976 K); SCM18 20.6 (LFP) |
| SOC stress `g(z)` levels | 0.7 / 1.0 / 1.75 | – | high level U(1.3, 1.9) | KEIL16 step near 60 % SOC |
| `k_y` (1000-cycle class) | 4.3e-3 | FCE^−½ | U(3.5e-3, 5.0e-3) | calibrated to 80 % at 1000 FCE (AP-BATT iPhone 15+; GOOG Pixel 8a+); 7.2e-3 for the 500-cycle class |
| `E_H` / `E_L` / `β_L` | 33 / 55.5 / 2.64 | kJ/mol, kJ/mol, h | E_H U(25, 45), E_L U(40, 70) | SCM18 (LFP; structure transferred) |
| `w_L` (plating share at 25 °C, 1C) | 0.3 | – | U(0.1, 0.5) | ASSUMPTION |
| `γ_R` (`r_growth = γ_R·q_loss`) | 2.5 | – | T(2, 2.5, 7) | SCH14 G8 (+41–65 % at 20 % fade); KEIL16 NCA up to 6.6 |
| Charge protocol | CC 1.0C to z = 0.8, then 0.3C to target | – | – | ASSUMPTION consistent with AP-WHY ("fast to 80 %, then trickle") |
| Charger loss `η_ch` | 0.90 | – | U(0.85, 0.95) | ASSUMPTION (thermal_aging §2.1) |
| JEITA rules | 0.5× current at 0–10 °C, none below 0 °C | – | – | TI-JEITA |

### 4.6 Validation device configurations (no re-tuning allowed)

**iPhone 14 Pro / 13 Pro (Bonn).**

* `Q = 3.32 / 3.16 Ah`: SCH26 measured initial capacities 3321–3326 / 3142–3181 mAh.
* `R0,ref` and `R1,ref` are scaled with `R·Q` held equal to RP-24's (60 mΩ·4.70/Q).
* All other parameters equal RP-24.
* Radios off (airplane mode); screen on and black: `s_scr = 1`, `APL = 0`, so only `P_d0` counts.

**LG HG2 cell (KOL20).**

* Cell only: no PCM, no `P_floor`; constant-current tests.
* Calibration on 25 °C data only:
  * NMC OCV re-fitted to the 25 °C C/20 curve;
  * `R0 + R1 = 0.038 Ω` from the 25 °C apparent resistance;
  * cell `C_th` and `R_th` from the 25 °C, 2C temperature rise.
* Everything at other temperatures is predicted with RP-24's activation energies.

---

## 5. Usage scenarios

### 5.1 Activity primitives (constant u)

Each column sets `u`. `—` means 0, off or not applicable. `Wi-Fi assoc.` means `w = 1`, and the cellular RAT is
camping (LTE) unless stated otherwise.

| Primitive | `s_scr` | `b` | `APL` | `f_r` | `x_c` | radio (`w`, `c`, `ρ`, `κ`) | `λ` (s⁻¹) | `d_dl` / `d_ul` (Mbit/s) | other | `ΔP_b` | P [PROTO] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SB: standby | 0 | — | — | — | 0 | 1, LTE, cell, burst | 12/3600 | — | — | 1.0 | 0.042 W |
| SB5: standby, 5G, heavy push | 0 | — | — | — | 0 | 0, NR, cell, burst | 60/3600 | — | — | 1.0 | 0.090 W |
| RD: reading | 1 | 0.4 | 0.7 | 60 | 0.05 | 1, LTE, wifi, burst | 1/30 | 0.05 / — | — | 3.0 | 0.504 W |
| SO: social scrolling | 1 | 0.4 | 0.7 | 120 | 0.25 | 1, LTE, wifi, burst | 0.5 | 2.0 / — | — | 3.0 | 1.089 W |
| SO-L: social over LTE | 1 | 0.4 | 0.7 | 120 | 0.25 | 0, LTE, cell, burst | 0.5 | 2.0 / — | — | 3.0 | 1.235 W |
| VW: video stream, Wi-Fi | 1 | 0.4 | 0.35 | 60 | 0.10 | 1, LTE, wifi, burst | 0.25 | 3.0 / — | `au = 1` | 2.0 | 0.596 W |
| VL / V5: video stream, LTE / 5G | 1 | 0.4 | 0.35 | 60 | 0.10 | 0, LTE or NR, cell, burst | 0.25 | 3.0 / — | `au = 1` | 2.0 | 0.761 / 0.867 W |
| NV: navigation | 1 | 0.8 | 0.6 | 60 | 0.20 | 0, LTE, cell, burst | 0.2 | 0.1 / — | `g = 1`, `P_sens = 47 mW` | 3.0 | 1.568 W |
| GM: 3D game, online | 1 | 0.6 | 0.5 | 120 | 0.80 | 1, LTE, wifi, burst | 20 | 0.2 / 0.05 | `au = 1` | 4.0 | 4.05 W demand (≈ 3.2 W capped) |
| VC: video call, LTE | 1 | 0.5 | 0.5 | 60 | 0.30 | 0, LTE, cell, rt | — | — | `c_f = 1`, `au = 1` | 3.0 | 2.573 W |
| CM: rear camera | 1 | 0.5 | 0.5 | 60 | 0.05 | 1, LTE, wifi, burst | 0 | — | `c_r = 1` | 4.0 | 1.921 W |

The assumptions in this table are labelled: λ values, bit rates (3 Mbit/s for a 720p stream), video APL 0.35,
navigation brightness 0.8 and holding times are [ASSUMPTION]s, each with UQ ranges in §8.

### 5.2 Scenario definitions

All scenarios use RP-24, start rested (`v1 = 0`) and run until the TTE event. Unless stated otherwise:
`T_amb = 25 °C`, `T0 = T_amb`, SOH = 1, `RSRP = −85 dBm`, and each scenario is run from `s0 ∈ {100, 50, 20} %`.

| ID | Name | u(t) | Environment / state |
|---|---|---|---|
| S1 | Standby (LTE + Wi-Fi) | SB, constant | – |
| S2 | Standby on 5G NSA with heavy push | SB5, constant; `T_tail` 10.3 s (21.4 s chained in UQ) | – |
| S3 | Light reading | RD, constant | – |
| S4 | Social scrolling, Wi-Fi, 120 Hz | SO, constant | – |
| S5a / S5b / S5c | Video streaming over Wi-Fi / 5G NSA / LTE | VW / V5 / VL, constant | – |
| S6 | Navigation (GNSS + LTE + bright screen) | NV, constant | – |
| S7 | Sustained 3D gaming | GM, constant | `T_amb ∈ {25, 35} °C` |
| S8 | Video call over LTE | VC, constant | – |
| S9 | Cold outdoor: navigation with photos | periodic, 30-min period: NV 10 min → CM 1 min (b = 0.9) → NV 19 min | `T_amb = −10 °C`, `T0 = 20 °C` (phone carried out from indoors); `k_env = 1` (in hand) |
| S10 | Aged battery | SO constant (S10a), and the S12 day (S10b) | `q_loss = 0.20`, `r_growth = 0.50` (SOH 0.8, R × 1.5); `T_amb ∈ {25, 0, −10, −20} °C` |
| S11 | Weak cellular signal | SO-L (S11a) and VC (S11c) at `RSRP ∈ {−85, −95, −105, −110, −120}` dBm; optional Wi-Fi cliff S11b at `RSSI = −85 / −90` dBm | – |
| S12 | Typical day (deterministic schedule from 07:00) | 07:00 RD 30 min · NV 45 · 7 × (SB 26 + SO 4) · VW 30 · SO 30 · 8 × (SB 26 + SO 4) · NV 45 · VC 30 · GM 60 · VW 30 · VW 90 · SO 60 · RD 60 · SB 8 h (23:00–07:00). Repeated until empty | – |
| S12-R | Stochastic day (CTMC) | 7 states {SB, RD, SO, VW, NV, GM, VC} | see below |

**S12 mean demand.** The schedule averages 0.552 W [PROTO]. Pixel 8a's published median-user profile implies
0.56 W (GPX), so S12 is a realistic "typical day". This is a plausibility check, not a fit.

**S12-R generator.**

* Mean holding times (min): SB 25, RD 6, SO 8, VW 20, NV 30, GM 25, VC 10 [ASSUMPTION].
* Jumps from SB go to (RD, SO, VW, NV, GM, VC) with probabilities (0.35, 0.35, 0.12, 0.08, 0.05, 0.05).
* Every active state returns to SB with probability 0.7. The other 0.3 is split over the remaining active states in
  proportion to the SB row.
* `Q_ij = (J_ij − δ_ij)/h_i`.
* Stationary distribution π = (0.597, 0.069, 0.092, 0.088, 0.089, 0.047, 0.019). This gives `P̄ = 0.590 W` [PROTO].

### 5.3 Baseline results [PROTO]

Fresh RP-24 from 100 %. The event is `cutoff` (`V ≤ V_cut`) or `brownout` (`V_pk ≤ V_sys`).

| Scenario | Mean P (W) | TTE | Event | s_strand | T at end |
|---|---|---|---|---|---|
| S1 standby | 0.042 | 428 h (17.8 d) | cutoff | 0.0 % | 25.1 °C |
| S2 5G standby, heavy push | 0.090 | 197 h (8.2 d) | cutoff | 0.0 % | 25.3 °C |
| S3 reading | 0.504 | 35.1 h | cutoff | 0.2 % | 26.6 °C |
| S4 social, Wi-Fi | 1.089 | 16.2 h | cutoff | 0.4 % | 28.5 °C |
| S5a video, Wi-Fi | 0.596 | 29.7 h | cutoff | 0.2 % | 26.9 °C |
| S5c video, LTE | 0.761 | 23.2 h | cutoff | 0.3 % | 27.4 °C |
| S5b video, 5G NSA | 0.867 | 20.4 h | cutoff | 0.3 % | 27.8 °C |
| S6 navigation | 1.568 | 11.2 h | cutoff | 0.6 % | 30.0 °C |
| S7 gaming, 25 °C | 4.05 → ≈ 3.1 (capped) | 5.6 h | cutoff | 1.1 % | 34.6 °C (skin ≈ 43 °C) |
| S7 gaming, 35 °C | 4.05 → ≈ 1.75 (capped) | 10.0 h | cutoff | 0.6 % | 40.5 °C |
| S8 video call, LTE | 2.573 | 6.8 h | cutoff | 0.9 % | 33.3 °C |
| CM rear camera | 1.921 | 9.1 h | brownout | 0.7 % | 31.2 °C |
| S9 cold outdoor (−10 °C) | ≈ 1.6 | 10.6 h (11.0 h at 25 °C) | brownout | 3.4 % | −4.9 °C |
| S10a aged, social, 25 / 0 / −10 °C | 1.089 | 12.9 / 12.5 / 11.8 h | brownout | 0.9 / 3.8 / 8.7 % | – |
| S10 aged, S9 profile at −10 °C | ≈ 1.6 | 7.8 h | brownout | 11.0 % | −4.8 °C |
| S11a social over LTE, −85 / −110 dBm | 1.235 / 1.359 | 14.3 / 13.0 h | cutoff | 0.5 % | – |
| S12 typical day, from 100 / 50 / 20 % | 0.552 | 34.7 / 11.7 / 7.7 h | cutoff | ≈ 1 % | – |
| S12 at 0 °C; aged at 25 °C; aged at 0 °C | – | 34.0 / 26.7 / 24.6 h | brownout / brownout / cutoff | – | – |
| S12-R CTMC (energy bucket, 4000 paths) | 0.590 | mean 30.6 h, sd 7.0 h, 5–95 % 19.3–42.1 h | – | – | – |

**Stranded-charge map** (S10a social load, `ΔP_b = 3 W`) [PROTO]. Each cell gives `s_strand` and, in parentheses,
TTE in hours:

| T_amb | SOH 1.0 | SOH 0.9 (R × 1.25) | SOH 0.8 (R × 1.5) |
|---|---|---|---|
| 25 °C | 0.4 % (16.2) | 0.5 % (14.5) | 0.9 % (12.9) |
| 10 °C | 0.7 % (16.1) | 1.4 % (14.4) | 2.1 % (12.7) |
| 0 °C | 1.4 % (16.0) | 2.5 % (14.2) | 3.8 % (12.5) |
| −10 °C | 2.8 % (15.8) | 5.0 % (13.8) | 8.7 % (11.8) |
| −20 °C | 6.5 % (15.1) | 13.3 % (12.6) | 23.6 % (9.9) |

**What the baseline already shows.**

* **Cold on a fresh modern pack costs little at phone-level power.** At −10 °C, 1.1 W, the loss is only −3 %.
  Self-heating keeps the battery 3.5 K above ambient, and R0 is only 60 mΩ.
* **Aging and cold interact.** The interaction raises stranded charge roughly tenfold at −20 °C.
* **Heat lengthens gaming TTE but costs performance.** At 35 °C throttling cuts gaming power and raises TTE by 78 %.
  The SoC then delivers only ≈ 50 % of the demanded throughput: `x_c` falls from 0.80 to ≈ 0.39, against ≈ 0.66
  when capped at 25 °C. These `x_c` values come from inverting E11 at the mean SoC power [DERIVED from PROTO].

---

## 6. TTE definition and numerical method

### 6.1 Definition

```
TTE(s0, scenario, θ) = inf{ t > 0 : min_i g_i(x(t), u(t)) ≤ 0 }
g1 = z                         (charge exhausted; never first in practice since z_e > 0)
g2 = V − V_cut                 (fuel gauge reports empty → orderly shutdown)
g3 = V_pk − V_sys              (instantaneous brown-out under a burst)
g4 = D                         (power infeasible; reachable only through demand steps)
g5 = T_shut − T                (thermal shutdown, 68 °C)
```

* **Conditional TTE.** "Remaining time" from an arbitrary state, `TTE(x_now)`, uses the same definition with the
  current state as the initial condition.
* **Event location at a switch.** `u` is piecewise constant, so g3 can jump below zero at a switching instant (for
  example, a camera start). The event is then located at the switching time.

### 6.2 Numerical method

**Stiffness.** Linearised eigenvalues:

| Mode | Eigenvalue |
|---|---|
| RC | `λ_v = −1/τ1`, from −1/18 s⁻¹ (45 °C) to −1/135 s⁻¹ (−20 °C) |
| Thermal | `λ_T ≈ −1/(k_env R_th C_th) ≈ −1/450 s⁻¹` |
| SOC | `λ_z ≈ +I·U'(z)/(Ũ·3600 Q) ~ +10⁻⁶ s⁻¹` (weakly unstable: a lower U means a higher I) |

The stiffness ratio is 10³–10⁵, which is mild. The `λ_z` growth is bounded by the events.

**Single trajectories (reference solver).**

* `scipy.integrate.solve_ivp(method="LSODA", rtol=1e-7, atol=[1e-9, 1e-7, 1e-5], max_step=60 s)`.
* Integrate segment by segment between input breakpoints, restarting the solver at each switch and carrying the
  state.
* Events g1–g5 are terminal with `direction = −1`.
* This is the [PROTO] configuration.

**Batches (Monte Carlo, Sobol).**

* Use the existing `src/phonebatt/integrate.py:rk4_batch`, vectorised over samples, with `h = 2 s`.
* Segment durations must be multiples of `h`.
* Event function: `g = min(z/0.01, (V − V_cut)/0.1, (V_pk − V_sys)/0.1, (T_shut − T)/10)`. It is NaN when `D < 0`;
  `rk4_batch` already treats NaN as an immediate event.
* The event type is the argmin component at the crossing.
* RK4 is stable on the real axis for `|hλ| < 2.785`, which allows `h ≤ 2.785·τ1,min ≈ 49 s` (τ1 ≈ 17.6 s at
  45 °C). We pick `h = 2 s` for accuracy at switches; the event-time error is O(h²) from linear interpolation.
* Standby scenarios may use `h = 10 s` after the convergence test below.

**Required numerical tests (CI unit tests).**

1. **Convergence.** `|TTE_RK4(h=2 s) − TTE_LSODA|/TTE < 0.2 %` on S4, S7, S9 and S10a at −10 °C, and the same
   tolerance between `h = 2 s` and `h = 1 s`.
2. **Energy closure.** With `v1(0) = 0`, the electrical energy identity is
   `∫ I·U dt = ∫ P dt + ∫ I²R0 dt + ∫ v1²/R1 dt + ½·C1·v1(t_end)²`. Here
   `∫ I·U dt = −3600·Q_max·∫ U dz` along the path. The identity must hold to a relative error below 10⁻⁴. The
   entropic term exchanges heat with the chemistry, not electrical energy, so it is excluded from this balance.
3. **Analytic agreement.** For constant P, isothermal (`C_th → ∞`) and `ΔP_b = 0`, the ODE TTE matches F2 to 0.5 %.
   [PROTO]: 0.00–0.40 % for P = 0.5–2.6 W.
4. **Physical sanity.** No NaN on feasible paths; `z` monotone; and `T → T_amb + k_env R_th(P + q_b)` at steady state.

### 6.3 Analytic and approximate TTE formulas

**F1. Energy bucket (level L0).** With `E_u(z0) = 3600·Q_max·∫_{z_e}^{z0} U(z,T_ref) dz`,

```
TTE0 = E_u / P̄
```

RP-24: `E_u(1) = 17.79 Wh`. Within 1 % of the ODE for P ≤ 1.1 W at 25 °C.

**F2. Exact quasi-steady quadrature.** This holds for constant P and T, with the RC branch equilibrated
(`v1 = R1·I`) and `R_Σ = R0 + R1`. Then `dz/dt = −P/(3600 Q V_qs(z))` integrates in closed form:

```
TTE = (3600·Q_max / P) · ∫_{z_e(P,T)}^{z0} V_qs(z) dz,     V_qs = [U + √(U² − 4 R_Σ P)]/2
z_e(P,T): V_qs(z_e) = V_cut,   or the brown-out root of F4 if that is larger
```

The inverse gives SOC(t) in closed form. [PROTO] errors against the full ODE are +0.01 % (0.50 W), −0.03 % (1.09 W),
−0.09 % (1.57 W) and −0.40 % (2.57 W). F2 is the cheap surrogate for Sobol runs (§9.3).

**F3. First-order resistive correction.**

```
TTE ≈ (E_u/P)·[1 − P·⟨R_Σ/U²⟩_E] − (3600 Q U(z_e)/P)·Δz_e,     Δz_e ≈ R_Σ(z_e) P /(V_cut·U'(z_e))
```

* `Π1 = P R_Σ/U²` is the resistive loss fraction: 0.5 % for RP-24 at 1 W and 25 °C, but ≈ 14 % at 4 W in the cold on
  an aged pack.
* Because `U'(z_e) ≈ 8.4 V` per unit SOC is so steep, the cut-off shift costs only ≈ 0.3 % at 25 °C.
* Together these explain why `R0` has a TTE elasticity of only −0.008 at room temperature (§9.1).

**F4. Brown-out threshold and stranded charge.** Set `V_pk = V_sys` in E9:

```
U(z*, T) = V_sys + v1 + R0(z*, T)·(P + ΔP_b)/V_sys,      v1 ≈ R1·P/U,      s_strand ≈ (z* − z_e)/(1 − z_e)
```

This is one scalar root.

* Check, fresh pack at −20 °C: `T_bat ≈ −16.5 °C`, `R0(0.1) = 264 mΩ`, so `U* = 3.65 V`, `z* ≈ 0.09` and
  `s_strand ≈ 6.8 %`. The ODE gives 6.5 % [PROTO].
* **Policy use.** The safe burst budget is `ΔP_b,safe(z,T,SOH) = V_sys(Ũ − V_sys)/R0 − P`. It drives the OS governor
  in §11.2.

**F5. Thermally capped load.** At thermal steady state, `T_s*` solves

```
T_s* = T_amb + R_skin·[P(T_s*) + q_b],     TTE_cap ≈ E_u / P(T_s*)
```

A useful bound is `P_cap ≈ max(P_nonSoC + P_fl, (T_thr − T_amb)/R_skin)`.

* 25 °C: `P_cap = 3.27 W`, so TTE ≈ 5.4 h. The ODE gives 5.6 h.
* The bound is only indicative when the sigmoid's width matters (35 °C).
* The step response of the bulk temperature is
  `T(t) = T_amb + R_th P + (T0 − T_amb − R_th P)·e^(−t/τ_th)` (HE19 Eq. 9).

**F6. Stochastic usage: renewal CLT for a CTMC power process.** The energy used up to time t is
`N(t) = ∫_0^t P_{X(s)} ds`, where `X` is the CTMC with generator Q, stationary distribution π and state powers
`P_i`. The first passage `TTE = inf{t : N(t) ≥ E_u}` satisfies, for large `E_u`,

```
TTE ≈ Normal( E_u/P̄ ,  σ²_∞·E_u/P̄³ ),       P̄ = Σ π_i P_i
σ²_∞ = 2 Σ_i π_i (P_i − P̄) g_i,   where   Q g = −(P − P̄·1),   π·g = 0          (Poisson equation)
Two-state on/off (rates α: off→on, β: on→off):  σ²_∞ = 2αβ(P_on − P_off)²/(α + β)³
```

[PROTO] for S12-R:

| | Mean | sd | Note |
|---|---|---|---|
| Analytic | 29.7 h | 7.3 h | CV 0.245 |
| Monte Carlo (4000 paths) | 30.6 h | 7.0 h | – |

The 3 % gap in the mean is the O(1) initial-state and overshoot correction.

* **Interpretation.** Usage randomness alone gives a TTE coefficient of variation of about 25 % when heavy states last
  tens of minutes.
* **Scaling.** The CV falls as `1/√E_u` and as `1/√(switching rate)`.
* **Parameter dependence.** σ²_∞ depends on θ only through the `P_i`, so `E[TTE|θ]` and `Var[TTE|θ]` cost almost
  nothing. This enables the variance decomposition in §8.

---

## 7. Validation plan

**Rule.** Every dataset below is used **only to test predictions**, never to set parameters, unless a step is
explicitly marked "calibration". All predictions with RP-24 defaults were computed before any comparison and are
**pre-registered** here [PROTO].

### V1. Phone-level TTE: Bonn iPhones (BONN26/SCH26, CC BY 4.0) — primary

**V1.1–V1.3. Blind predictions.** Configurations as in §4.6. Observed values are median runtimes to the last logged
1 %. The predicted TTEs below run to the model shutdown event. The official comparison uses the model time at
`s_iOS = 1 %`, which is about 1 % of the run earlier; the errors change by less than 1 percentage point.

| Device, scenario | Predicted P (W) | Predicted TTE (h) | Observed (h) | Error |
|---|---|---|---|---|
| 14 Pro, black-screen idle | 0.270 | 46.4 | 55.0 | −15.7 % |
| 14 Pro, rear camera | 1.670 | 7.38 | 8.43 | −12.4 % |
| 14 Pro, GPS best | 0.520 | 24.0 | 48.5 | −50.5 % |
| 14 Pro, camera + GPS (held-out additivity test) | 1.920 | 6.41 | 8.38 | −23.6 % |
| 13 Pro, idle | 0.270 | 44.2 | 50.8 | −13.1 % |
| 13 Pro, camera | 1.670 | 7.02 | 9.53 | −26.3 % |
| 13 Pro, GPS | 0.520 | 22.9 | 29.8 | −23.3 % |
| 13 Pro, camera + GPS | 1.920 | 6.09 | 8.33 | −26.9 % |

**Pass criteria.**

* (a) |ε_TTE| ≤ 15 % for idle and camera.
* (b) The observed device range lies inside the 90 % prediction interval from the §8 parameter UQ.
* (c) The model satisfies the additivity law `1/TTE_cxl = 1/TTE_cam + 1/TTE_GPS − 1/TTE_idle` to within 1 %. [PROTO]:
  0.3 %, since the ECM is nearly linear in P here.

The pre-registered expectation is that idle passes narrowly or fails narrowly, and GPS fails: 2019 Android GNSS
figures overstate the 2022 iPhone 14 Pro roughly 8×. We report this as a finding (GNSS "matters surprisingly
little" on modern hardware) and not as a reason to re-tune.

**V1.4. SOC-shape test (structural).** At constant power, coulombic SOC falls faster as V drops. The display map D1
predicts a drain slope of 31.1 → 28.8 → 27.9 min per % over the 99–50, 50–20 and 20–1 % windows. That follows from
the OCV window means 4.059, 3.756 and 3.645 V [PROTO]. D2 (state of energy) predicts a flat slope. Bonn observes
31.1 / 30.8 / 32.4 min per % (≤ 6 % variation). Compute the slope-error metric for D1 and D2 on all idle, GPS and
camera traces, and pick the display map that the data favour. The iOS top-hold (`h = 0.065`) is applied to both.

**V1.5. Calibration variant.** Fit only the idle floor (`ΔP_awake + P_d0`) to the 14 Pro idle median, then
re-predict camera, GPS and camera + GPS. This shows how far one baseline calibration carries.

### V2. Cell temperature behaviour: LG HG2 (KOL20, CC BY 4.0) — primary

* Calibration (25 °C only) as in §4.6.
* Predict, for T = −20, −10, 0, 10, 40 °C at 0.5C and 2C: delivered capacity ratio `f_Q`, energy ratio `f_E`, mean
  voltage, and the V(t) RMSE.
* Targets (validation_data plan):

| Quantity | Observed | Pass criterion |
|---|---|---|
| `f_Q(0 °C)`, 0.5C / 2C | 0.91 / 0.88 | |Δf_Q| ≤ 0.05 for T ≥ 0 °C |
| `f_Q(−20 °C)`, 0.5C / 2C | 0.73 / 0.49 | within ±0.10 |
| `f_Q(−20 °C)`, C/20 | 0.90 | must stay high (the loss is kinetic) |
| V(t) | – | RMSE ≤ 50 mV for T ≥ 0 °C |

* **Additional test.** Compare our effective activation energy with the apparent HG2 value of 26.4 kJ/mol. Our R0 + R1
  mix gives ≈ 20 kJ/mol between 25 and −20 °C.
* **Dynamic test.** One UDDS file at 25 °C and one at 0 °C: RMSE_V.

### V3. Aging-to-TTE map: NASA PCoE B0005 (NASA07; licence "not specified", documented in the paper) — secondary

Feed the *measured* fade (−28.7 %) and Re growth into `Q_max` and `R`, then predict the aged TTE at 2 A. Observed:
55.8 → 39.7 min. This tests the structural claim that aging enters as Q↓ and R↑, independent of the aging-rate law.
Cold checks use B0047 and B0055 (4 °C).

### V4. RC dynamics under pulsed load: Oxford ExampleDC (OX17, ODbL) — optional

Drive the cell ECM, scaled to 0.74 Ah at 40 °C, with the 1 Hz Artemis current (peaks of 6.8C). Report RMSE_V.

### V5. Physics reference: PyBaMM DFN (BSD-3; Ai2020 / Marquis2019 sets) — structural

Run constant-power discharges (0.1–1C equivalent) at 25, 0 and −20 °C, and 10 s pulses. Compare TTE and V(t) with the
1RC and 2RC ECM using the same OCV. Acceptance: TTE difference ≤ 3 % at ≥ 0 °C. This judges whether 1RC is enough.

### V6. Published-specification plausibility band (cited facts, not open data)

Predicted RP-24 power against the published implied power. The tolerance is ±25 % (undisclosed brightness and
network conditions).

| Published case | Published P | RP-24 predicted P | Within band? |
|---|---|---|---|
| Wi-Fi browsing at 150 nits (NBC22–24) | 1.06–1.31 W | social, Wi-Fi: 1.09 W | ✓ |
| Streamed video (APL14P) | 0.62 W | VW: 0.596 W | ✓ |
| Local video (NBC) | 0.59–0.78 W | – | ✓ |
| Typical day (GPX, Pixel 8a) | 0.56 W | S12: 0.552 W | ✓ |
| Max load or gaming (NBC) | 2.4–5.8 W | capped gaming: ≈ 3.1 W | ✓ |
| Idle at minimum brightness, no WLAN (NBC) | 0.29–0.71 W | 0.34 W | ✓ |
| Audio with screen off (APL14P) | 0.165 W | 0.226 W | ✗ (+37 %; expected — our awake overhead assumes a full AP wakelock, while iOS offloads audio to a DSP) |

### V7. Qualitative phenomenology (LEE20, HE19, AP-TEMP) — check

Configure a 2015-class pack with high resistance. Values are assumptions:
`Q = 2.7 Ah [unverified]`, `R0 = 0.12 Ω`, `R1 = 0.03 Ω`, ≈ 3.35 W UI load with 3 W bursts.

* The model must produce brown-out at 15–35 % displayed SOC at ≤ 0 °C. [PROTO]: 23 % at 0 °C.
* It must give a TTE loss at 0 °C of the same order as LEE20's −15 %. [PROTO]: −32 %; the device parameters are
  unknown, so this counts as order-of-magnitude only.
* The charge stranded at shutoff must be recovered after warming (LEE20 C12: +330 mAh). This happens automatically,
  because z is not reset.

---

## 8. Uncertainty quantification plan

**Sources of uncertainty and how each is treated.**

1. **Parameters θ.** About 50 uncertain inputs with the §4 distributions, treated as independent. `V_sys` is coupled
   to `V_cut` through `Δ_sys`.
2. **Usage.** The CTMC (S12-R) gives the path randomness. Holding times are scaled by `U(0.5, 2)` in an outer loop.
   Brightness, APL and λ are drawn per state visit from `±30 %` uniform bands [ASSUMPTION].
3. **Environment.** `T_amb` follows scenario-specific distributions. For "outdoor winter" that is `U(−15, 5) °C`.
   `k_env` follows its §4.2 distribution.
4. **Structure.** Model variants M1–M12 (§9.4), reported as a separate spread and not mixed into θ.
5. **Observation.** For Bonn only: the display map (D1 vs D2), 1 % quantisation, and the top hold.

**Propagation.**

* Scrambled Sobol' QMC (`scipy.stats.qmc.Sobol`), `N = 2048` per scenario, through the batched RK4. For S1/S2,
  standby uses `h = 10 s`.
* Percentile intervals: report the median, the 90 % PI (5–95 %), the coefficient of variation, `Pr(TTE < t*)` and the
  `s_strand` quantiles.
* The target-time event is `t* = 16 h` for S12 started at 07:00: "dies before 23:00".

**Variance decomposition, parameters vs usage.** Apply the law of total variance,
`Var(TTE) = Var_θ(E[TTE|θ]) + E_θ(Var[TTE|θ])`. Evaluate `E[TTE|θ]` and `Var[TTE|θ]` with F6, which is analytic per
θ sample, and confirm on 200 θ × 20 paths of the full ODE.

**Deliverables.**

* TTE table for scenarios × `s0` ∈ {100, 50, 20} % × `T_amb` ∈ {25, 0, −10} °C: median [5 %, 95 %].
* Violin plots.
* A statement of the form: "for a typical day, usage randomness contributes X % of the TTE variance and hardware
  parameters Y %".

---

## 9. Sensitivity analysis plan

### 9.1 Local elasticities

**Definition.** `ε_θ = ∂ln TTE/∂ln θ`, by central difference with Δ = ±5 %. Voltages are reported as % TTE per +10 mV.

**Applicability.** Use constant-load scenarios only. Under periodic schedules with brown-out, TTE is **non-smooth in
θ**: the event snaps to the next heavy segment. This produces spurious elasticities, for example `ε_Q < 0` in S9
[PROTO]. For those scenarios use the OAT ranges of §9.2 or the CTMC-stationary versions.

[PROTO] elasticities:

| Parameter | S4 social, 25 °C | S10a aged social, −10 °C | S7 gaming, 25 °C |
|---|---|---|---|
| `Q_n` | +1.00 | +1.00 | +1.01 |
| `P_SoC,max` | −0.39 | −0.41 | −0.15 (capped) |
| `P_c,ref` | −0.15 | – | −0.02 |
| `P_d0` | −0.14 | −0.15 | −0.01 |
| `ΔP_awake` | −0.09 | – | −0.01 |
| `R0,ref` | −0.008 | **−0.29** | −0.013 |
| `E_a,R0` | 0.00 | **−0.26** | 0.00 |
| `R1,ref` | −0.001 | −0.02 | 0.00 |
| `ΔP_b` | 0 (not binding) | **−0.20** | 0 |
| `R_skin` | 0 | 0 | **+0.79** |
| `R_th` | 0 | +0.03 (self-heating helps) | −0.01 |
| `V_cut` (per +10 mV) | −0.21 % (secant over +0.17 V) | 0 (not binding) | −0.33 % (secant over +0.17 V) |
| `V_sys` (per +10 mV) | 0 at default (V_cut binds) | **−0.70 %** | 0 at default (V_cut binds) |
| `δU` OCV offset (per +10 mV) | +0.37 % | +0.98 % | +0.41 % |

For S4 and S7, a +5 % step in `V_sys` (to 3.47 V) lifts it above `V_cut` and switches the binding event. Such
regime switches must be reported as ranges (§9.2), not as derivatives.

**Reading.** At room temperature TTE is an energy-budget problem: Q and the load parameters matter, and the cell
resistance does not. In the cold on an aged pack it becomes a voltage problem: R0, E_a, ΔP_b and V_sys. During
sustained gaming it is a thermal-design problem: R_skin.

### 9.2 One-at-a-time tornado

* Move each θ to its range endpoints with the others at default.
* Scenarios: S4, S6, S7, S10a (−10 °C), S12.
* Rank by |ΔTTE| and plot tornado diagrams (Figure F7).
* Add usage factors to the same chart, with these ranges:
  * brightness: `b ± 0.2`;
  * dark mode: `APL × 0.25`;
  * refresh: 60 vs 120 Hz;
  * background rate: `λ_bg × {0.5, 5}`;
  * RSRP: −85 → −115 dBm;
  * RAT: Wi-Fi / LTE / 5G;
  * `T_amb`: −10 → 35 °C.
* This puts hardware uncertainty and user behaviour on one axis.

### 9.3 Global variance-based (Sobol)

**Method.** Saltelli sampling with Jansen estimators (`SALib`, already in `requirements.txt`).

* `k = 24` grouped factors: radio parameters grouped per RAT, display {P_d0, P_c,ref, γ_b}, and the RC group
  {R1, τ1, E_a,R1}.
* `N = 1024`, giving `N(2k + 2) = 51,200` model evaluations per scenario.

**Surrogates.**

* For constant-load isothermal-dominated scenarios (S1–S6, S8, S11), use **F2 + F4 + F5** as an exact
  quasi-steady surrogate (≤ 0.4 % error).
* Confirm on S7, S9 and S10a with the full ODE at `N = 256`, with bootstrap 95 % CIs.

**Report.** First-order `S_i`, total `S_Ti`, and second-order `S_ij` for the triple {T_amb, SOH, ΔP_b}. The
interaction share `1 − Σ S_i` measures the cold × aging × burst coupling.

### 9.4 Structural sensitivity (model-form uncertainty)

Run every scenario under each variant and report ΔTTE and Δs_strand relative to M0 (Table T4).

| ID | Variant | What it tests |
|---|---|---|
| M0 | baseline (this spec) | – |
| M1 | energy bucket F1 (no ECM) | is a voltage model needed? |
| M2 | constant-current load (I = P/V_nom) instead of constant power | the PMIC constant-power assumption |
| M3 | isothermal at T_amb (no thermal ODE; no throttling) | self-heating in the cold; throttling in gaming |
| M4 | 2RC (2.5 mΩ/6.7 s + 14.1 mΩ/47 s, scaled by R·Q) | RC order |
| M5 | R0 split: 35 mΩ temperature-independent PCM + 25 mΩ Arrhenius cell | where the resistance lives |
| M6 | no burst event (mean-voltage cut-off only) | the peak-power mechanism |
| M7 | two-node thermal (skin node, τ_s = 3 min, KANG19) | quasi-steady skin assumption |
| M8 | KiBaM available/bound charge (c = 0.166, k′ = 0.122 min⁻¹; JH09, structure only) | rate-capacity and recovery |
| M9 | radio tail as ODE state `dr/dt = λ(1 − r) − r/T_tail` instead of the exact Poisson fraction | quasi-steady tail |
| M10 | brightness law γ_b = 1 (CH13) and γ_b = 2.9 (DH21-derived) | display law |
| M11 | display map D1 vs D2 | output definition |
| M12 | SoC leakage `P_leak = P_leak,25·exp(k_lk(T_s − 25 °C))`, `P_leak,25 = 0.1 W`, `k_lk ∈ [0, 0.04] K⁻¹` [unverified] | thermal runaway of power in heat |

### 9.5 Usage-fluctuation sensitivity

1. **Rate scaling.** With π fixed, scale all CTMC rates by `κ ∈ {0.25, 0.5, 1, 2, 4}`. By F6 the mean is invariant,
   while the sd scales as `κ^(−½)`. Verify by Monte Carlo.
2. **Burstiness at fixed mean.** Compare a smooth 1.2 W load with a 50/50 mix of 0.4 W and 2.0 W of the same mean,
   for switching periods of 10 s to 1 h. The energy-bucket TTE is unchanged. The ECM differs only through `I²R` and
   the brown-out; quantify this at 25 °C and at −10 °C on an aged pack.
3. **Order effect (path dependence).** Take the same S12 activities in reversed order, with heavy use at low SOC.
   This is irrelevant in the energy bucket but matters in the ECM in the cold, because of `R0(z)` and F4.
4. **Initial SOC.** TTE(s0) is nearly linear for s0 ≥ 0.2 at 25 °C. Report the non-linearity in the cold, where the
   effective empty point `z*` rises.

---

## 10. Aging plan

### 10.1 Slow-state equations

These follow thermal_aging §2.4(c). Days are the time unit for the calendar term, and FCE are the unit for the cycle
term.

```
dq_cal/dt = [k_c·A⁻¹(T; E_cal)·g(z)]² / (2 q_cal)                              (t in days)     (E13)
dq_cyc/dn = k_y²·[ w_H·A⁻¹(T; E_H) + w_L·A(T; E_L)·exp(β_L(c_ch − c_ref)) ]² / (2 q_cyc)          (E14)
dn/dt = |I| / (2·3600·Q_n)        (FCE, counting charge and discharge)
q_loss = q_cal + q_cyc;   r_growth = γ_R·q_loss;   Q_max = Q_n(1 − q_loss)
g(z) = 0.7 + 0.3·σ((z − 0.2)/0.02) + 0.75·σ((z − 0.6)/0.03)    (smooth KEIL16 step)
```

* **Arrhenius sign.** `A⁻¹(T;E) = exp[−(E/R_g)(1/T − 1/T_ref)]` grows with heat. The plating term uses
  `A(T;E_L)`, which grows with cold (SCM18 inverse Arrhenius).
* **Constants.** `w_H = 0.7`, `w_L = 0.3`, `c_ref = 1C`, and `c_ch` is the charge C-rate (0 while discharging).
* **Initialisation.** `q(0) = 10⁻⁴` avoids the singular initial rate.
* **Calibration** (25 °C, 1 FCE/day). The mean calendar stress is `ḡ ≈ 1.45` for a daily 100 → 30 % cycle
  [DERIVED estimate: ≈ 60 % of the time above z = 0.6]:

| Target | Calibrated value |
|---|---|
| `q_cal(1 yr, g = 1) = 2.7 %` | `k_c = 1.41e-3` |
| `q_cal(1000 d) ≈ 6.5 %` | consistent with thermal_aging §2.4 |
| 20 % total loss at N_80 = 1000 FCE | `k_y = 4.3e-3` |
| 20 % total loss at N_80 = 500 FCE | `k_y = 7.2e-3` |

### 10.2 Daily protocol for aging runs

* **Discharge.** The S12 day: about 13.3 Wh/day, ≈ 0.75 FCE/day.
* **Charge.** Plug in at 23:00. CC at 1C to z = 0.8, then 0.3C to the target, then hold until 07:00.
* **Charging heat.** Add `(1 − η_ch)·P_in` to E7. The fast-state ODE runs during charging with `I < 0`.
* **Compression.** Run one full day (≈ 4×10⁴ RK4 steps), compute the day's increments of `q_cal` and `q_cyc` by
  integrating E13 and E14 along it, then advance the slow state by one day.
* **Speed-up.** Recompute the fast day only when SOH has changed by more than 0.5 %: piecewise-frozen
  "cycle-averaged" forcing. This is exact to first order in the time-scale ratio.

### 10.3 Aging scenarios

| ID | Charging and environment | Question |
|---|---|---|
| A1 | charge to 100 %, plugged 8 h, 25 °C | baseline SOH(t) over 3 years |
| A2 | 80 % charge limit (AP-OBC, SAM) | calendar benefit vs a daily TTE cost of −20 % usable charge |
| A3 | optimized charging: hold 80 % until 1 h before unplug | how much of the A2 benefit A3 keeps at zero TTE cost |
| A4 | A1 at 35 °C ambient (hot climate) | heat as the dominant stress (×2.1 per +10 K) |
| A5 | 2 h/day at 45 °C (car dashboard), otherwise A1 | episodic heat |
| A6 | 2C fast charge at 10 °C vs 0.5C | the plating term; JEITA justification |
| A7 | 1 vs 2 FCE/day (heavy user) | throughput vs calendar |

**Outputs.**

* `SOH(t)` and `R(t)/R_new` for each scenario.
* TTE(S12) vs years.
* `s_strand(T_amb = 0 °C, S10 load)` vs years: when does a winter brown-out start?
* Years until SOH = 0.8.

**Coupling statement for the paper.** Aging reduces TTE through two channels:

* **Q↓**: proportional and benign.
* **R↑**: nonlinear, and invisible at 25 °C. At 25 °C, `ε_R0 = −0.008`, so a 50 % R increase costs < 1 % TTE. In
  the cold it dominates (§5.3 map).

The R↑ channel is why "battery health 80 %" phones die unexpectedly in winter.

---

## 11. Recommendations plan

### 11.1 Ranking user behaviours by TTE gain

**Method.** Apply each behaviour change to S12 (and to S12-R with UQ). Report ΔTTE as a median and 90 % PI.

[PROTO] first ranking, based on mean demand for S12 at 25 °C on a fresh pack (baseline P̄ = 0.552 W):

| Change | ΔTTE |
|---|---|
| Replace the 1 h game with standby | **+43 %** (gaming is 30 % of the day's energy) |
| Dark theme (APL × 0.25 on UI states) | +9.0 % |
| Brightness −0.1 on all screen-on states | +5.3 % |
| Wi-Fi (or offline maps) instead of LTE for the call and navigation data | +4.9 % |
| Brightness capped at 0.5 | +4.6 % |
| 60 Hz cap (no 120 Hz) | +4.6 % |
| Voice call instead of video call | +4.5 % |
| GNSS off when not navigating | +2.9 % |
| Halve background sync | +0.3 % (surprisingly small with modern 4G tails) |
| Weak signal all day (−110 dBm) | −4.0 % |
| Standby on 5G NSA instead of LTE + Wi-Fi | −5.1 % |

**Energy shares in S12 by component** [PROTO, demand-based]: SoC 40.5 %, display 30.3 %, radio 12.3 %, awake
overhead 7.2 %, suspend floor 3.6 %, GNSS 2.8 %, audio 1.5 %, camera 1.1 %, sensors 0.5 %.

**What drains most.** Sustained SoC load (gaming), bright light-themed screens, real-time cellular uplink (video calls,
worse at weak signal), 5G tails in standby, and cold plus aging plus bursts, which strands charge.

**What matters surprisingly little.** The standby floor (16 h of standby use only 4.5 % of the day's energy), halving
background sync, sensors, GPS on modern chips, battery self-heating, self-discharge, rate-capacity effects, and cell
resistance at room temperature.

**Cold-weather advice** (from §5.3 and F4):

* Keep the phone in an inner pocket, so `T_bat` stays high.
* Avoid camera or flash bursts below 20 % in the cold.
* Aged packs: at −10 °C an SOH-0.8 phone strands 8.7 % of displayed charge, against 2.8 % when new.

### 11.2 OS power-management strategies, each testable in the model

| # | Strategy | Model mechanism | How we quantify it |
|---|---|---|---|
| O1 | **Peak-power governor** (brown-out avoidance; cf. Apple performance management, AP-BATT) | cap `ΔP_b ≤ max(0, 0.9·P_cut,sys(x) − P)` using F4 online | [PROTO] aged social at −10 °C: `ΔP_b` 3 → 1.5 W cuts `s_strand` from 8.7 % to 2.6 % (TTE +6.1 %); at `ΔP_b` 1.0 W, 1.6 %. Report the TTE gain vs the "performance lost" (∫ capped burst energy) |
| O2 | Thermal-aware frame-rate and quality cap before throttling | lower `x_c` and `f_r` when `T_s → T_thr` | gaming TTE and frames delivered vs reactive throttling (E8) |
| O3 | Radio policy: prefer LTE for low-throughput traffic, use 5G only for bulk; batch background syncs into existing tails | lower `P_tail` and `E_pro` exposure; λ_eff ↓ | XU20 measured −24.8 % app energy; we report S2 and S12 gains |
| O4 | Adaptive dark theme and brightness cap at low SOC or high APL | `APL`, `b` | S12 from 20 % |
| O5 | Model-based remaining-time display with an uncertainty band | EKF on (z, v1) using measured V, I, T; online P̄ and σ²_∞ (F6) | "4 h 10 min (3 h 20 – 5 h 00)": calibrate the interval coverage on S12-R |
| O6 | Cold-aware fuel gauge | remap displayed 0 % to `z*(T, SOH, ΔP_b)` (F4) and warn "X % unavailable until warm" | removes surprise shutdowns at displayed > 0 |
| O7 | Charge-limit and optimized charging defaults; no fast charge below 10 °C | E13, E14 | A2, A3, A6 |

---

## 12. Generalisation to other portable devices

**Same equations, different parameters.** The model transfers through four scalings and four dimensionless groups.

**Scalings.**

* Pack resistance scales as `R ∝ 1/Q` at fixed chemistry. RP-24 has `R0·Q = 0.28 Ω·Ah`, including the PCM.
* `C_th ∝ m`.
* `R_th, R_skin ∝ 1/(h·A_surface)`.
* `P_disp ∝ A_display` (LCD has no APL term).
* Multi-cell packs with `N_s` cells in series and `N_p` in parallel: `U → N_s U`, `R → (N_s/N_p) R`, `Q → N_p Q`, and
  the constant-power load is shared.

**Dimensionless groups.**

| Group | Meaning | RP-24 values |
|---|---|---|
| `Π1 = P·R_Σ/Ū²` | resistive loss fraction | ≈ 0.005 at 1 W, 25 °C |
| `Π2 = P_dem/P_sust` | throttling index; capped if > 1 | 1.24 for gaming |
| `Π3 = τ_th/TTE` | thermal memory | 0.02 |
| `Π4 = ΔP_b/(P_cut,sys − P)` | brown-out margin; failure if → 1 | – |

* **Laptops** have small Π1 and large thermal mass.
* **Smartwatches and earbuds** have tiny `C_th`, so τ_th is about 1 minute [unverified scaling]. Their BLE and sensor
  floor dominates, and Π4 can be large in the cold because the cells are small.
* **E-readers** have a bistable display, with `P_disp` ≈ 0 between refreshes.
* **Drones** [unverified orders of magnitude] have Π1 ~ 0.05–0.2 and high Π4. Their hover load is close to constant
  power, so A1 holds well.

We place each device class on a Π1–Π4 map (Figure F12). The parameters needed per class are listed, and any product
numbers are marked unverified, since no device data were collected.

---

## 13. Figures and tables for the paper

| # | Type | Content |
|---|---|---|
| F1 | Diagram | Model architecture: u(t) → component powers → battery-side P → ECM (E2–E6) ↔ thermal (E7) ↔ throttle (E8); burst/brown-out (E9); slow SOH layer (E13–E14). Feedback arrows T→R, T_s→P_SoC, SOH→Q,R |
| T1 | Table | RP-24 parameter table, condensed from §4 (symbol, value, source, UQ range) |
| T2 | Table | Activity primitives and scenarios (§5.1–5.2) |
| F2 | Figure | Cell maps: U(z,T), R0(z,T) and the brown-out power surface `P_cut,sys(z,T)` with operating points of S4, S6 and S7 |
| F3 | Figure | Stacked battery-side power breakdown per scenario (§3) next to the published-spec band (V6) |
| F4 | Figure | Trajectories `s(t)`, `V(t)`, `T(t)`, `T_s(t)` for S1, S5a, S7 (throttle visible) and S9 (cold brown-out) |
| T3 | Table | TTE median [90 % PI] for scenarios × `s0` ∈ {100, 50, 20} % × `T_amb` ∈ {25, 0, −10} °C |
| F5 | Figure | Validation: (a) Bonn predicted vs observed TTE with PIs; (b) SOC-slope test D1 vs D2; (c) HG2 `f_Q(T)` at 0.5C and 2C, model vs data |
| F6 | Figure | UQ violins per scenario, with variance split into parameters vs usage (F6 decomposition) |
| F7 | Figure | Tornado diagrams for S4 and S10a at −10 °C, hardware and behaviour factors on one axis |
| F8 | Figure | Sobol `S_i` / `S_Ti` heatmap (factors × scenarios), with the {T_amb, SOH, ΔP_b} interaction share |
| T4 | Table | Structural sensitivity: ΔTTE and Δs_strand for variants M1–M12 |
| F9 | Figure | Stranded-charge heatmap `s_strand(T_amb, SOH)`, with and without the O1 peak-power governor |
| F10 | Figure | Aging: SOH(t) and TTE(S12)(t) for charging policies A1–A7 |
| F11 | Figure | Recommendation ranking: ΔTTE bars (median and PI) for user behaviours and OS strategies |
| F12 | Figure | Generalisation map of device classes in (Π1, Π4) and (Π2, Π3) |

---

## 14. Implementation notes for `src/phonebatt`

**Modules.** All are vectorised over a trailing sample axis of shape `(n_samples,)`.

* `params.py`: dataclasses `CellParams`, `ThermalParams`, `DeviceParams`, `RadioParams`, `AgingParams`, plus
  `RP24 = default()`. `UQ_SPEC` holds `{name: (dist, a, mode, b)}`. Constructors for `IPHONE14PRO`, `IPHONE13PRO`
  and `HG2`.
* `usage.py`: `UsageState` (the §1.4 fields), `PRIMITIVES` (§5.1), `Schedule` (list of `(duration_s, UsageState)`,
  repeating), `CTMC(states, hold_s, J)`.
* `power.py`: `power(u, Ts, p) -> (P_mean, P_peak, breakdown)`, implementing E1 and E8–E12. No state.
* `cell.py`: `ocv(z, T)`, `dUdT(z)`, `R0(z, T, rg)`, `R1(T, rg)`, `tau1(T)`, and
  `solve_current(U_tilde, R0, P) -> (I, V, D)` (E2), `peak_voltage(...)` (E9), `z_empty(p)`.
* `model.py`: `rhs(t, x, u, env, p)`, `events(...)`, `display_soc(z, mode)`, `run_scenario(...) -> Result(tte,
  event, s_strand, traj)`.
* `integrate.py`: the existing `rk4_batch`, plus a thin `solve_ivp` wrapper for segmented schedules.
* `analytic.py`: F1–F6 (`tte_bucket`, `tte_quadrature`, `brownout_threshold`, `thermal_cap_fixed_point`,
  `ctmc_tte_moments`).
* `aging.py`: E13–E14, the daily-protocol driver, and the charging model.

**Unit tests.**

* All four tests in §6.2.
* `R0(0.5, T_ref) == 60 mΩ`; `U(1) == 4.4496 V`; `z_e == 0.0263`.
* The `power()` values in the §5.1 "P [PROTO]" column, to ±0.5 mW.
* The §5.3 TTEs to ±1 %.
* The two-state F6 formula against Monte Carlo (5 %).
* `rk4_batch` NaN-event path when D < 0 (a forced 20 W step).

The [PROTO] scratch implementation is in `docs/design/scripts/pf_proto_model.py` and `pf_proto_checks.py`. It is
reference only, not production code.

---

## 15. Assumption register and limitations

| # | Assumption | Evidence or test |
|---|---|---|
| A1 | Constant-power load at the battery terminals (PMIC regulation) | standard (device_power §2.1); structural test M2 |
| A2 | One thermal node plus a quasi-steady skin offset | HE19 single-exponential settling (adj. R² ≈ 1); test M7 |
| A3 | All device power becomes heat in the phone | thermal_aging §2.1; light and RF are neglected |
| A4 | Bursts ≪ τ1, so only R0 acts during a burst | HU12 promotion 0.26 s, τ1 ≥ 18 s |
| A5 | Poisson traffic; exact tail fraction `1 − e^(−λT)` | our derivation; test M9 |
| A6 | Time-scale separation between SOH and discharge | weeks vs hours |
| A7 | No OCV hysteresis | absorbed in δU ±20 mV |
| A8 | 4.45 V OCV built by stretching a 4.2 V curve | no measured high-voltage OCV found (cell_ecm caveat) |
| A9 | Single Arrhenius on the whole pack R0 | matches LEE20 ratios; test M5 |
| A10 | Component powers add | Bonn 0.15 % additivity residual |
| A11 | Display SOC map D1, or D2 | tested in V1.4; output-only |
| A12 | No self-discharge in TTE runs | 0.1 %/month (ZIL19) |
| A13 | Aging laws use LFP/NMC structure transferred to LCO | no open phone-LCO aging model (thermal_aging caveat); E_cal 20–60 kJ/mol in UQ |

**Known limitations.**

* Absolute modern GNSS, idle-floor and 5G-modem powers carry the largest transfer risk from 2012–2020 devices.
  V1 quantifies this.
* No open dataset gives phone-level runs at controlled cold. Cold validation is therefore at cell level (V2) plus
  qualitative phone phenomenology (V7).
* Charging behaviour is modelled only enough for the aging layer.

---

## 16. Source keys (full records in the research notes)

### docs/research/device_power.md

| Key | Reference |
|---|---|
| CH10 | Carroll & Heiser, USENIX ATC 2010 |
| CH13 | Carroll & Heiser, APSys 2013, doi:10.1145/2500727.2500734 |
| HU12 | Huang et al., MobiSys 2012, doi:10.1145/2307636.2307658 |
| NA21 | Narayanan et al., SIGCOMM 2021, doi:10.1145/3452296.3472923 |
| XU20 | Xu et al., SIGCOMM 2020, doi:10.1145/3387514.3405882 |
| GU24 | Gupta et al., MobiCom 2024 poster, doi:10.1145/3636534.3695905 |
| DI13 | Ding et al., SIGMETRICS 2013 |
| CH15 | Chen et al., SIGMETRICS 2015, doi:10.1145/2745844.2745875 (via course slides) |
| DZ11 | Dong & Zhong, arXiv:1101.1240 |
| DH21 | Dash & Hu, MobiSys 2021 (via Purdue press release) |
| CX13 | Chen et al., HotMobile 2013 |
| HE20 | Herglotz et al., IEEE Access 2020, doi:10.1109/ACCESS.2020.2986580, CC BY 4.0 |
| KW19 | Karki & Won, arXiv:1910.13041 |
| SU14 | Sun et al., ICDCS 2014, doi:10.1109/ICDCS.2014.13 |
| JE12 | Jensen et al., VTC-Fall 2012 |
| BH19 | Bhat et al., DATE 2019, arXiv:1904.09814 |
| CHI15 / CHI23 | Chiriac et al., Electronics Cooling 2015 / 2023 |
| AT20 | AnandTech 2020, via 9to5Google and Ubergizmo (M/L) |
| LX | Linux kernel v6.6 energy-model documentation and device trees (GPL/BSD/MIT) |
| AOSP | Android power-profile documentation (illustrative values) |

### docs/research/cell_ecm.md

| Key | Reference |
|---|---|
| RUN | PyBaMM 26.10.1.1 runs in this project |
| PYB-AI20 | Ai et al., JES 2020, doi:10.1149/2.0122001JES, CC BY |
| CRM06 | Chen & Rincón-Mora, IEEE TEC 2006, doi:10.1109/TEC.2006.874229 |
| PLETT-N | Plett, ECE5710 lecture notes |
| BARC22 | Barcellona et al., J. Energy Storage 2022, doi:10.1016/j.est.2022.104658 |
| SUN23 | Sun et al., Energies 2023, doi:10.3390/en16207142 |
| PINE | PinePhone battery datasheet QZ01-396172 |
| BU802a | Battery University BU-802a |
| PIX9 | Google Pixel 9 spec sheet, 2024 |
| QNOVO17 | Qnovo blog, 2017 |
| LIN26 | Lin et al., Advanced Materials 2026, doi:10.1002/adma.202523570 |
| QC-FG / QC-SMB | Qualcomm/AOSP kernel device-tree bindings |
| OVEJ19 | Ovejas & Cuadras, Scientific Reports 2019, doi:10.1038/s41598-019-51474-5 |
| BERN85 | Bernardi et al., JES 1985, doi:10.1149/1.2113792 |
| JH09 | Jongerden & Haverkort, IET Software 2009, doi:10.1049/iet-sen.2009.0001 |
| RV01 | Rakhmatov & Vrudhula, ICCAD 2001 |

### docs/research/thermal_aging.md

| Key | Reference |
|---|---|
| KANG19 | Kang et al., MobiCom 2019, doi:10.1145/3300061.3300128 |
| PU23 | Pu & Wu, arXiv:2304.06512 |
| HE19 | He et al., ICCPS 2019 |
| LEE20 | Lee, He & Shin, MobiSys 2020, doi:10.1145/3386901.3389024 |
| MAL99 | Maleki et al., JES 1999, doi:10.1149/1.1391704 |
| SCH14 | Schmalstieg et al., J. Power Sources 2014, doi:10.1016/j.jpowsour.2014.02.012, via CA21 (Clarke & Alonso 2021, doi:10.2514/1.C036209) |
| SCM18 | Schimpe et al., JES 2018, doi:10.1149/2.1181714jes, CC BY 4.0 |
| KEIL16 | Keil et al., JES 2016, doi:10.1149/2.0411609jes, CC BY 4.0 |
| ZIL19 | Zilberman et al., J. Energy Storage 2019, doi:10.1016/j.est.2019.100900 |
| BU802b | Battery University BU-802b |
| AP-BATT / AP-TEMP / AP-OBC / AP-WHY | Apple support pages 101575, 118431, 108055 and Batteries – Why Lithium-ion |
| GOOG | Google Pixel Help 15738128 |
| SAM | Samsung Battery protection, One UI 6.1 |
| TI-JEITA | Qian, TI SLYT365, 2010 |
| UL-TS | UL brief on IEC 62368-1, 2010 |
| AOSP-TH | AOSP thermal mitigation documentation |

### docs/research/validation_data.md

| Key | Reference |
|---|---|
| BONN26 / SCH26 | Schweizer et al. 2026, dataset doi:10.60507/FK2/LSJGWW (CC BY 4.0); paper Sensors 26(10) 2923, doi:10.3390/s26102923 |
| KOL20 | Kollmeyer et al. 2020, LG HG2 data, doi:10.17632/cp3473x7xv.3, CC BY 4.0 |
| NASA07 | Saha & Goebel 2007, NASA PCoE Battery Data Set (licence not specified) |
| OX17 | Birkl & Howey 2017, doi:10.5287/bodleian:KO2kdmYGg, ODbL |
| APL14P | Apple iPhone 14 Pro technical specifications, support.apple.com/en-us/111849 |
| NBC22–24 | Notebookcheck reviews of iPhone 14 Pro Max, Pixel 8 and Galaxy S24 |
| GPX | Google Pixel hardware tech specs, support.google.com/pixelphone/answer/7158570 |
| N9M22 | 9to5Mac 2022 report of iPhone battery Wh figures |
