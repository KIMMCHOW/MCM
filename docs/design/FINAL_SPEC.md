# FINAL SPEC: coupled electro-thermal-aging model of smartphone battery drain

2026 MCM Problem A, *Modeling Smartphone Battery Drain*. Final model specification, written 2026-11-17 by the lead
modeller after judging the two design proposals
(`docs/design/proposal_physics_first.md` = **PF**, `docs/design/proposal_parsimony_stochastic.md` = **PS**).
This document supersedes both. The implementation in `src/phonebatt` follows this file, not the proposals.

## 0. Status, conventions and provenance

**Status.** Implementation-ready. Every equation is explicit. Every parameter has a value, a unit, a source (or an
engineering-estimate justification with a plausibility check) and an uncertainty-quantification (UQ) range. Every
scenario has an explicit input `u(t)`, and every analysis has named output files and JSON keys (§14).

**Tags.**

| Tag | Meaning |
|---|---|
| `[KEY]` | A source key, for example `CH13` or `LEE20`. Keys resolve to the research notes in `docs/research/` and to the reference list in §17. |
| **[RUN]** | A number computed from open data or PyBaMM in the research step. It is recorded in `docs/research/*.md` and `docs/research/data/`. |
| **[DERIVED]** | Our own arithmetic on cited numbers. The arithmetic is shown. |
| **[ASSUMPTION]** | A modelling choice with no measured source. It always has a UQ range and a plausibility check. |
| **[CAL]** | The single calibrated constant of the model (§4.3, V1.1). |
| **[unverified]** | Background knowledge that could not be checked in this project. It is given as a range. |
| **[PROTO]** | Computed in this session by `docs/design/scripts/final_proto.py`, which implements exactly the equations below (SciPy LSODA, rtol 1e-7). PROTO numbers are **regression targets** for `src/phonebatt` (§15.3). They are not validation results. Reproduce with `python3 -I docs/design/scripts/final_proto.py [section]`. |

**Units and signs.** SI units, with three exceptions: capacity `Q` is in A·h, energy is in W·h where stated, and
the factor 3600 appears explicitly. Temperature `T` is in kelvin inside equations and in °C in tables. Current
`I > 0` means discharge. Every power is a **battery-terminal** power. Constants: `R_g = 8.314 J mol⁻¹ K⁻¹`,
`T_ref = 298.15 K`.

**Two Arrhenius factors (named to avoid the sign clash found in PS).**

```
a_R(T; E)   = exp[ +(E/R_g)(1/T − 1/T_ref) ]     resistance-like: > 1 when colder than 25 °C
a_deg(T; E) = exp[ −(E/R_g)(1/T − 1/T_ref) ]     degradation-like: > 1 when hotter than 25 °C
```

### 0.1 What the final model takes from each proposal, and what was fixed

| Element | Taken from | Change in the final spec |
|---|---|---|
| Hierarchy L0 → L3 for structural sensitivity and closed forms | PS (M0–M4) | Renamed L0–L3. L2 is the full coupled ECM, and L3 adds aging. KiBaM is a structural variant (X9) |
| Main model: constant-power ECM with cancellation-free current root, Arrhenius `R(z,T)`, Bernardi heat with entropic term, skin-temperature throttling | PF (shared by PS M3) | Kept unchanged |
| Two voltage thresholds: fuel-gauge cut-off `V_cut` and brown-out `V_sys` under a sub-second burst envelope | PF | Kept. Explicit long bursts (≥ 10 s) are now part of `u(t)` and are checked against `V_cut` (from PS) |
| Mechanistic component loads (OLED with APL, DVFS SoC, radio with exact Poisson tail fraction) | PF | Three fixes: (i) a nested C-DRX "active" timer in the radio formula; (ii) display content power scales with panel area; (iii) camera increment re-derived from CX13 |
| Usage CTMC (hub-and-spoke × connectivity), PDMP mean via Dynkin, CLT variance, inverse-Gaussian law, hub variance attribution, energy-share elasticities | PS | Kept. Recomputed with the final state powers |
| Calendar SOC-stress step `g(z)` | PF (smooth KEIL16 step) | Multiplied above z = 0.8 by SCH14's linear voltage law (high-voltage LCO hold) |
| Cycle aging | both (each lumped two mechanisms in one state) | Split into SEI (`q_H`, total throughput) and plating (`q_L`, charge throughput) sub-states as in SCM18 Eq. 23. All aging states are integrated in squared form `d(q²)/dt = k²`, which needs no `q(0) = ε` hack |
| Bonn floor | PF (blind prior) and PS (tuned, but presented as validation) | One explicit calibrated constant, `ΔP_awake` [CAL]. The blind-prior run is reported too |
| SOC output | PF/PS (coulombic) | Coulombic displayed SOC for initial conditions and the shutdown SOC. The **energy-fraction SOC** of the problem glossary ("energy remaining … compared to its full capacity") is reported as well. The two differ by ≤ 2.4 points [PROTO] |
| Recommendation metric | PF used mean demand; PS used the M0 CTMC mean | Formalised: levers are ranked on smooth metrics (`E_u/P̄` of S12 and the Dynkin CTMC mean). The L2 TTE of a fixed schedule is non-smooth in parameters (§12.1) |

The full critique with scores is in Appendix B.

---

## 1. Model overview

### 1.1 Hierarchy (structural ladder)

All levels share the same load model `P(t) = P(u(t), T_s(t); θ)` (§3). Each level adds one physical mechanism.

| Level | States | Adds | Role |
|---|---|---|---|
| **L0** energy bucket | `E` (W·h) | `dE/dt = −P/3600` | closed forms; elasticity = −(energy share) |
| **L1** quasi-static ECM | `z` | OCV `U(z)`, series `R_Σ = R0 + R1` at fixed `T_amb`, constant-power current, cut-off `V_cut` | exact TTE quadrature; ohmic and cut-off losses |
| **L2** dynamic electro-thermal ECM | `z, v1, T` | RC polarisation, lumped thermal node, Arrhenius `R(z,T)`, entropic heat, skin throttling, burst envelope and brown-out | **main TTE model** (fresh or at fixed SOH) |
| **L3** L2 + aging | `z, v1, T` (fast); `Q_cal, Q_H, Q_L` (slow) | calendar fade, SEI cycle fade, plating fade, resistance growth | **full main model**; SOH frozen within one discharge |

"Main model" in the paper means L3. For a single discharge it is L2 evaluated with `(Q_max, R)` from the slow
states.

### 1.2 Time scales

| Process | Time scale | Source |
|---|---|---|
| Load bursts (radio promotion, CPU boost, camera start) | 0.01–1 s | HU12 (LTE promotion 260 ms), LEE20 |
| C-DRX inactivity timer | 0.08–0.2 s | XU20 (80/100 ms) |
| RC polarisation `τ1` | 20 s (40 °C) to 135 s (−20 °C) | [RUN] cell_ecm a.3 |
| RRC tail | 0.24 s (Wi-Fi), 5–21 s (cellular) | HU12, NA21, XU20 |
| Thermal `τ_th = R_th·C_th` | 7.5 min (5–15) | HE19 |
| Usage switching | minutes | [ASSUMPTION] §5.3 |
| Discharge (TTE) | 4 h to 18 days | this model |
| Aging | weeks to years | AP-BATT, GOOG, SCH14 |

Consequences: bursts are handled **algebraically** (they see `R0` only, with `v1` frozen); RRC tails are handled by
an exact **time-average** (Poisson formula); RC and thermal dynamics are **ODE states**; SOH is a **slow state**,
frozen within a discharge.

### 1.3 State vector

| State | Symbol | Unit | Range | Initial condition |
|---|---|---|---|---|
| Coulombic SOC relative to the present maximum capacity `Q_max` | `z` | – | [0, 1] | `z0 = z_e + s0(1 − z_e)` from the displayed SOC `s0` (E11) |
| RC polarisation voltage | `v1` | V | ≥ 0 on discharge | 0 (rested) |
| Lumped phone/battery temperature | `T` | K | – | `T0 = T_amb` (soaked) unless the scenario sets it |
| Squared calendar loss | `Q_cal = q_cal²` | – | ≥ 0 | 0 (new) or the scenario value |
| Squared SEI cycle loss | `Q_H = q_H²` | – | ≥ 0 | 0 |
| Squared plating loss | `Q_L = q_L²` | – | ≥ 0 | 0 |

Derived quantities: `q_loss = √Q_cal + √Q_H + √Q_L`, `SOH = 1 − q_loss`, `Q_max = Q_n(1 − q_loss)`, and
`r_growth = γ_R·q_loss`.

### 1.4 Inputs

`u(t)` is piecewise constant. It is either a deterministic schedule of activity primitives (§5.1) or a CTMC path
(§5.3).

| Field | Symbol | Domain | Meaning |
|---|---|---|---|
| Screen on | `s_scr` | {0, 1} | display powered |
| Awake | `a` | {0, 1} | application processor held awake; `a = max(s_scr, a_flag)` |
| Brightness slider | `b` | [0, 1] | OS slider fraction |
| Average picture level | `APL` | [0, 1] | linear-light mean pixel level |
| Dark-theme multiplier | `m_APL` | (0, 1] | 1 by default; 0.25 in the dark-theme lever |
| Panel refresh rate | `f_r` | {60, 90, 120} Hz | effective rate |
| SoC throughput demand | `x_c` | [0, 1] | fraction of maximum CPU+GPU throughput |
| Wi-Fi associated | `w` | {0, 1} | |
| Cellular camping RAT | `c` | {none, LTE, NR} | NR = 5G NSA low band |
| Data path | `ρ` | {wifi, cell, none} | RAT that carries the traffic |
| Traffic class | `κ` | {burst, rt} | Poisson bursts with timers, or continuous real-time media |
| Burst arrival rate | `λ` | s⁻¹ | traffic bursts or background wake-ups |
| Mean DL / UL throughput | `d_dl`, `d_ul` | Mbit/s | |
| Signal | `RSRP` | dBm | cellular reference signal received power |
| GNSS duty | `g` | [0, 1] | |
| Rear / front camera | `c_r`, `c_f` | {0, 1} | |
| Audio / Bluetooth audio | `au`, `bt` | {0, 1} | |
| Sensor power | `P_sens` | W | |
| Burst-envelope amplitude | `ΔP_b` | W | sub-second peak above the mean (§3.7) |
| Ambient temperature | `T_amb(t)` | K | environment |
| Enclosure factor | `k_env` | – | multiplies `R_th` (case, hand, pocket) |

### 1.5 Outputs

* `z(t)`, the displayed SOC `s(t)` (D1, E11), and the energy-fraction SOC `ς(t)` (D2, E11), which matches the
  problem glossary.
* `V(t)`, `I(t)`, `P(t)` and the component breakdown `P_k(t)`; `T(t)`, `T_s(t)`, the throttle factor `σ(t)` and the
  delivered SoC performance ratio `ρ_perf(t)` (E13).
* `TTE`, the event type ∈ {`charge`, `cutoff`, `brownout`, `infeasible`, `thermal`}, and the displayed SOC at
  shutdown `s_off = s(TTE)` ("stranded charge"; it comes back after warming, LEE20 C12).
* For stochastic usage: the TTE distribution, its quantiles and `Pr(TTE < t*)`.
* For aging: `SOH(t)`, `R/R_new`, years to SOH = 0.8, TTE and `s_off` against age.

### 1.6 Assumption register

| # | Assumption | Evidence or test |
|---|---|---|
| A1 | The PMIC makes the phone a constant-**power** load at the battery terminals | device_power §2.1; structural variant X4 (constant current) |
| A2 | Component powers add | Bonn additivity residual 0.15 % / −0.10 % (BONN26 [RUN]); V1.3 hold-out |
| A3 | One thermal node; skin hot spot is quasi-steady, `T_s − T_amb = (R_skin/R_th)(T − T_amb)` | HE19 single-exponential settling (adj. R² ≈ 1); X7 two-node variant |
| A4 | All device power ends up as heat inside the phone | thermal_aging §2.1 (emitted light and RF neglected) |
| A5 | Sub-second bursts see only `R0`; `v1` is frozen during a burst | burst ≤ 1 s ≪ `τ1` ≥ 20 s |
| A6 | Burst envelope: in an activity with `ΔP_b > 0`, a burst of that size occurs often enough that the first time it *would* brown out, it does | interactive bursts every few seconds (HU12 promotions, LEE20 peak currents); X7/X8 variants remove it |
| A7 | Traffic bursts are Poisson; radio timers reset on each burst | our derivation on HU12/NA21 state machines; X10 explicit RRC ODE |
| A8 | Time-scale separation: SOH frozen within one discharge | weeks vs hours |
| A9 | OCV hysteresis neglected; absorbed in an offset `δU` ±20 mV | cell_ecm a.2 (OVEJ19 mechanism; magnitude unverified) |
| A10 | 4.45 V OCV = Ai2020 4.2 V curve + 0.242·z^3.86 stretch | cell_ecm a.2 (construction; no measured 4.45 V curve found) |
| A11 | One Arrhenius factor on the whole pack `R0` (cell + PCM + contacts) | matches LEE20 phone-level R(−5)/R(25) = 2.0–3.3; X6 split variant |
| A12 | Pack resistance scales as `R·Q = const` across phones of the same chemistry | cell_ecm a.3 (DFN `R0·Q` ≈ 85–95 mΩ·Ah) |
| A13 | Usage is a Markov chain (exponential sojourns) | parsimony and closed forms; X15 semi-Markov variant |
| A14 | Aging laws: SCM18 mechanism structure, SCH14/KEIL16 stress factors, transferred to LCO | no open phone-LCO aging model (thermal_aging §4); wide UQ on `E_cal`, `β_L`, `γ_R` |
| A15 | Self-discharge neglected in TTE runs | 0.1 %/month (ZIL19) vs ≥ 0.2 %/h standby drain |

---

## 2. Governing equations

### 2.1 Power demand at the battery terminals

```
P(t) = P_susp + a·ΔP_awake + (1 − a)·λ·E_wake + P_disp + P_SoC(T_s) + P_net
       + g·P_G + c_r·ΔP_cam + c_f·ΔP_camf + au·ΔP_aud + bt·ΔP_BT + P_sens                         (E1)
P_pk(t) = P(t) + ΔP_b(u(t))                                                                       (E1b)
```

* All parameters are battery-side values, so they already contain the PMIC loss. `η_PMIC = 0.83` (CH13) is used
  only to convert rail-level sources (LX, CH13 rails) and is never applied twice.
* The term `(1 − a)·λ·E_wake` is the CPU energy of background wake-ups while the processor is otherwise asleep.

### 2.2 Current under constant power, and feasibility

With `Ũ = U(z,T) − v1`, the power balance `P = V·I` and `V = Ũ − R0·I` give `R0 I² − Ũ I + P = 0`. The physical
(high-voltage) root, written without cancellation, is

```
D = Ũ² − 4·R0·P,     I = 2P / (Ũ + √D),     V = (Ũ + √D)/2                     (requires D > 0)    (E2)
```

* `P_max = Ũ²/(4R0)` is the largest deliverable power, at `V = Ũ/2`.
* Because `V_cut = 3.40 V > Ũ/2` (≈ 1.8–2.2 V), the cut-off event is always reached before `D < 0` under smooth
  demand. `D ≤ 0` can only be reached by a demand step. The RHS then uses the clamped `I = Ũ/(2R0)` and event g4
  fires.
* The binding power limit is the cut-off-limited power `P_cut(V*) = V*(Ũ − V*)/R0` (thermal_aging T3).

### 2.3 Charge conservation

```
dz/dt = −I / (3600·Q_max),     Q_max = Q_n·(1 − q_loss)                                           (E3)
```

Coulombic efficiency is 1 on discharge and `η_c = 0.99` on charge (PLETT-N). Charging appears only in the aging
layer (§2.12).

### 2.4 Polarisation (1RC Thevenin)

```
dv1/dt = (R1(T)·I − v1) / τ1(T)                                                                   (E4)
```

The 2RC form (2.5 mΩ / 6.7 s + 14.1 mΩ / 47 s at 25 °C, [RUN] cell_ecm a.3) is structural variant X5.

### 2.5 Open-circuit voltage

```
U(z,T) = −b0·e^(−b1 z) + a0 + a1 z + a2 z² + a3 z³ + ΔV_hv·z^p + δU + (T − T_ref)·k_S·∂U/∂T(z)     (E5)
```

* Chen–Rincón-Mora form (CRM06) fitted to the Ai2020 LCO/graphite equilibrium OCV ([RUN], RMSE 7.0 mV, max
  18 mV): `b0 = 0.60117 V`, `b1 = 36.426`, `a0 = 3.6193 V`, `a1 = 0.44311 V`, `a2 = −0.48605 V`, `a3 = 0.63127 V`.
* High-voltage stretch `ΔV_hv = 0.242 V`, `p = 3.86` (cell_ecm a.2 construction): `U(1) = 4.4496 V`, mean OCV
  3.8699 V [PROTO]. Nearly no change below z = 0.3.
* Entropic coefficient `∂U/∂T(z)`: linear interpolation in the 21-point [RUN] table
  (`docs/research/data/cell_ecm_pybamm/pybamm_ocv_fit.json`, key `dUdT_cell_mV_per_K`), clamped outside [0, 1]:

| z | 0 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ∂U/∂T (mV/K) | −0.4485 | −0.7057 | −0.5343 | −0.4982 | −0.4787 | −0.4213 | −0.3722 | −0.3439 | −0.3068 | −0.2655 | −0.2482 |

| z | 0.55 | 0.60 | 0.65 | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 | 0.95 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|
| ∂U/∂T (mV/K) | −0.2412 | −0.2482 | −0.2281 | −0.2047 | −0.1798 | −0.1528 | −0.1246 | −0.0966 | −0.0705 | −0.0459 |

  Values below z ≈ 0.45 extrapolate the Ai2020 graphite fit outside its stated validity (cell_ecm caveat). They
  are covered by the multiplier `k_S ~ U(0.5, 1.5)`.
* Check values [PROTO]: `U(0.1) = 3.644`, `U(0.2) = 3.694`, `U(0.3) = 3.728`, `U(0.5) = 3.815`,
  `U(0.8) = 4.088`, `U(1) = 4.450 V`; `z_e = U⁻¹(V_cut; T_ref) = 0.0263`.

### 2.6 Resistances

```
R0(z,T) = R0∞·[1 + A_R·e^(−z/z0)]·a_R(T; E_a,R0)·(1 + r_growth),   R0∞ = R0,ref / (1 + A_R·e^(−0.5/z0))   (E6a)
R1(T)   = R1,ref·a_R(T; E_a,R1)·(1 + r_growth)                                                        (E6b)
τ1(T)   = τ1,ref·a_R(T; E_a,τ)                                                                         (E6c)
```

`R0(0.5, T_ref) = R0,ref` by construction. Check values [PROTO]:

| T (°C) | 40 | 25 | 10 | 0 | −10 | −20 |
|---|---|---|---|---|---|---|
| R0 at z = 0.5 (mΩ) | 42.4 | 60.0 | 88.1 | 116.6 | 157.6 | 218.1 |
| R1 (mΩ) | 5.6 | 10.0 | 19.0 | 30.3 | 50.0 | 86.0 |
| τ1 (s) | 20.0 | 30.0 | 47.0 | 65.1 | 92.6 | 135.3 |

Plausibility: `(R0 + R1)(−5 °C)/(R0 + R1)(25 °C) ≈ 2.5`, inside LEE20's phone-level 2.0–3.3.

### 2.7 Lumped thermal balance (Bernardi heat, Newton cooling)

```
C_th·dT/dt = P + I²R0 + v1²/R1 − I·T·k_S·∂U/∂T(z) − (T − T_amb)/(k_env·R_th)                      (E7)
T_s = T_amb + (R_skin/R_th)·(T − T_amb)                                                            (E8)
```

* `P = V·I` is the power delivered to the phone; all of it becomes heat inside the phone (A4).
* `I²R0 + v1²/R1` is the irreversible battery heat. The capacitor term `C1·v1·dv1/dt` is stored, not dissipated.
* `−I·T·∂U/∂T` is the reversible (entropic) heat. It is positive on discharge, because `∂U/∂T < 0` for LCO, and
  comparable to the Joule heat at ≈ 1 A (cell_ecm finding 8).
* At steady state, `T − T_amb = k_env·R_th·q` and `T_s − T_amb = k_env·R_skin·q`, where `q` is the total heat.

### 2.8 Thermal throttling

SoC demand `P_SoC,dem(x_c)` comes from E13. Above the floor `P_fl` it is scaled continuously by a logistic factor of
the skin temperature:

```
P_SoC = P_SoC,dem                                                   if P_SoC,dem ≤ P_fl
P_SoC = P_fl + (P_SoC,dem − P_fl)·σ((T_lim − ΔT_thr − T_s)/w_thr)   otherwise;   σ(y) = 1/(1 + e^(−y))   (E9)
```

The chain `T → T_s → P → I` is explicit, so there is no algebraic loop.

### 2.9 Burst envelope and brown-out

During a sub-second burst `v1` is frozen (A5):

```
D_pk = Ũ² − 4R0·P_pk,   V_pk = (Ũ + √D_pk)/2   (V_pk := 0 if D_pk ≤ 0)
brown-out  ⇔  V_pk ≤ V_sys = V_cut − Δ_sys   ⇔   P + ΔP_b ≥ P_cut,sys = V_sys(Ũ − V_sys)/R0              (E10)
```

* `V_cut` (3.40 V) is the fuel-gauge cut-off on the load-averaged voltage (QC-FG example; LEE20 "≈ 3.4 V").
* `V_sys` (3.30 V) is the instantaneous minimum system voltage (QC-SMB example). Mapping it to a brown-out under
  sub-second peaks is an [ASSUMPTION] with `Δ_sys ~ U(0, 0.20) V`.
* Explicit bursts longer than about 10 s (for example the 30 s camera-and-upload burst of S9) are part of `u(t)`.
  They are checked against `V_cut` like any other load.

### 2.10 Aging coupling (within one discharge)

`q_loss` and `r_growth` are constants within a TTE run. They scale `Q_max` (E3) and multiply `R0` and `R1` (E6).

### 2.11 SOC outputs and initial conditions

```
D1 (coulombic displayed SOC, used for s0 and s_off):  s = (z − z_e)/(1 − z_e),   z_e = U⁻¹(V_cut; T_ref) = 0.0263
D2 (energy-fraction SOC, problem glossary):           ς = ∫_{z_e}^{z} U(ζ,T_ref)dζ / ∫_{z_e}^{1} U(ζ,T_ref)dζ    (E11)
iOS top-hold map (Bonn comparison only):              s_iOS = floor_1%( min(1, s/(1 − h)) ),   h = 0.065
```

* `max|s − ς| = 0.024` at s ≈ 0.6 [PROTO]. At `s = 0.5 / 0.2 / 0.1` the energy fraction is
  `ς = 0.477 / 0.187 / 0.092`.
* `h = 0.065` is the iOS "100 %" hold measured on the Bonn traces (validation_data b.4, [RUN]).
* The shutdown SOC is `s_off = s(z(TTE))`.

### 2.12 Aging layer (L3 slow states), squared-loss form

Every √-law loss `q = k·√(stress-weighted exposure)` is integrated as the **linear** ODE `d(q²)/dt = k²·(stress)²`.
This is the rate-based integral of SCM18 (after Thomas et al.). It is exact for time-varying stress and has no
singularity at `q = 0`.

```
dQ_cal/dt_d = [ k_c·a_deg(T; E_cal)·g(z) ]²                                     (t_d in days)          (E17)
dQ_H/dn     = [ k_H·a_deg(T; E_H) ]²,          dn/dt = |I| / (2·3600·Q_n)        (FCE, charge + discharge) (E18)
dQ_L/dn_ch  = [ k_L·a_R(T; E_L)·e^(β_L(c_ch − c_ref)) ]²,   dn_ch/dt = max(−I, 0)/(3600·Q_n),   c_ch = max(−I,0)/Q_n  (E19)
q_loss = √Q_cal + √Q_H + √Q_L,    r_growth = γ_R·q_loss,    Q_max = Q_n(1 − q_loss)                       (E20)
```

**SOC stress of calendar aging** (KEIL16 graphite step × SCH14 voltage law above z = 0.8):

```
g(z) = [0.7 + 0.3·σ((z − 0.2)/0.02) + 0.75·σ((z − 0.6)/0.03)] · h_HV(z)^κ_HV
h_HV(z) = max{1, (7.543·U(z,T_ref) − 23.75) / (7.543·U(0.8,T_ref) − 23.75)}                          (E16)
```

* The step part reproduces KEIL16's plateaus (≈ 0.7 below 20 % SOC, 1.0 at 20–60 %, ≈ 1.75 above 60 %).
* `h_HV` is SCH14's calendar voltage factor (G1, `α_cap ∝ 7.543V − 23.75`), applied only above the 80 % plateau
  where the 4.45 V LCO cell leaves the voltage window of KEIL16. It gives `h_HV(1) = 1.385` [DERIVED]. Applying it
  above 4.2 V extrapolates SCH14 beyond its tested window, so `κ_HV ~ U(0, 1.5)`. NEL18 reports qualitatively that
  LCO cells held at 4.35–4.4 V fade severely.
* [PROTO] `g(0.1, 0.5, 0.7, 0.8, 0.9, 1.0) = 0.702, 1.026, 1.724, 1.749, 2.043, 2.423`.

**Charging model (aging layer only).**

```
CC at c_cc = 1.0C up to z = 0.8, then c_cv = 0.3C up to the target z_top ∈ {1.0, 0.8}                  (E21)
JEITA (TI-JEITA): c ≤ 0.5C when the cell is at 0–10 °C; no charging below 0 °C or above 60 °C
Cell temperature while charging: T_ch = T_amb + ΔT_ch;  in use: T_use = T_amb + k_env·R_th·P̄_day          (E22)
```

* The CC/CV split follows Apple's "fast to 80 %, then trickle" (AP-WHY); the rates are an [ASSUMPTION].
* `ΔT_ch = 8 K` is an [ASSUMPTION]; KANG19 measured ≈ 37 °C on the surface of an idle phone while fast charging.
* The full fast model may be run during charging (`I < 0`, charger heat `(1 − η_charger)·P_in` added to E7).
  The baseline aging driver uses the quasi-static day profile of §11.2, which is cheaper and was used for [PROTO].

### 2.13 Complete fast system (L2/L3 within one discharge) and well-posedness

```
x = (z, v1, T);   given u(t), T_amb(t), θ, q_loss:
dz/dt  = −I(x,u)/(3600·Q_max)
dv1/dt = (R1(T)·I(x,u) − v1)/τ1(T)
dT/dt  = [P(u,T_s) + I²R0 + v1²/R1 − I·T·k_S·∂U/∂T − (T − T_amb)/(k_env R_th)] / C_th
with I from E2 using P(u, T_s(T)) (E1, E8, E9)
```

* On the feasible set `F = {D > 0, 0 ≤ z ≤ 1}` the right-hand side is locally Lipschitz in `x`. The logistic, the
  exponentials and `√D` are smooth for `D > 0`, and the `∂U/∂T` interpolant is piecewise linear. It is piecewise
  continuous in `t`, because `u` is piecewise constant. A unique solution therefore exists until the first exit
  from `F` or the first event.
* `z` is strictly decreasing whenever `P > 0`. `T` stays within `T_amb + k_env·R_th·(P_max + q_batt,max)`.

---

## 3. Component power sub-models (battery side)

### 3.1 Floor

* Suspend floor `P_susp` (SoC suspend, RAM self-refresh, PMIC, no radios): always on.
* Awake overhead `a·ΔP_awake`: the application processor is awake and the display is off. It also applies to
  screen-off activities that hold a wakelock (audio, tracking).
* Background wake-ups while asleep: `(1 − a)·λ·E_wake`, with `E_wake = 0.10 J` [DERIVED: ΔP_awake ≈ 0.1 W × an
  assumed ≈ 1 s wake].

### 3.2 Display (OLED)

```
P_disp = s_scr·[ P_d0 + p_c·A_disp·(b/0.5)^γ_b·(m_APL·APL)/APL_ref + κ_RR·max(f_r − 60 Hz, 0) ]          (E12)
```

* **Physics.** OLED pixel power is linear in linear-light RGB (DZ11, R² > 0.95), so the content term scales with
  APL. Black pixels cost only the fixed overhead `P_d0` (CH10, CH13).
* **Brightness.** With `γ_b = 2`, the luminance fraction is `ℓ = b²`. Panel power is proportional to luminance, so
  a test "at 150 nits" on a panel with typical maximum `L_max` uses `b = (150/L_max)^(1/γ_b)`.
* **Calibration points.** At `b = 0.5`, `APL = 0.7`, `A = 97.4 cm²`, the content power is
  `p_c·A = 0.25 W` (GU24 display rail ≈ 0.30 W Pixel 7a, ≈ 0.45 W S23+, totals including `P_d0`). Full-white
  maximum-brightness power is `0.12 + 0.25·4/0.7 = 1.55 W` [PROTO], inside device_power's 1.0–2.0 W.
* **Dark-mode check [PROTO].** In the S4 social state, dark theme (APL 0.55 → 0.15) saves 5.1 / 8.6 / 12.6 / 33.2 %
  of state power at b = 0.3 / 0.4 / 0.5 / 1.0. DH21 reports 3–9 % at 30–50 % brightness and 39–47 % at 100 %. The
  model sits slightly above the band at b = 0.5 and slightly below at b = 1. `γ_b ~ U(1.0, 2.9)` spans CH13's
  near-affine slider (γ ≈ 1) and the DH21-derived γ ≈ 2.9 (PF [DERIVED]).
* **LCD variant** (generalisation only): `P_LCD = s_scr·[P_panel + c1·b + c2·b²]`, content-independent and convex
  (CH10, HE20).

### 3.3 SoC (CPU + GPU + memory) with DVFS

The Linux energy model gives `P = C·V²·f` (LX). With a utilisation-tracking governor, `f ∝ x_c`, and the OPP
voltage is affine in `f`, so `V = V_max·(v0 + (1 − v0)x_c)`. Hence

```
P_SoC,dem(x_c) = P_SoC,max · x_c · [v0 + (1 − v0)·x_c]²                                               (E13)
ρ_perf = x_eff/x_c,   x_eff solves P_SoC,dem(x_eff) = P_SoC   (delivered/demanded throughput under throttling)
```

* `v0 = 0.446` is [DERIVED] from the MT8183 A73 OPP table: slope 0.293 V/GHz, intercept 0.468 V, `V_max` 1.05 V
  (device_power S3).
* `P_SoC,max = 5.0 W` is CH13's CPU 2.85 W + GPU 1.42 W rails × 1.2 PMIC factor.
* Check values: `P_SoC,dem(0.05, 0.10, 0.25, 0.30, 0.50, 0.80) = 0.056, 0.126, 0.427, 0.562, 1.307, 3.163 W`
  [PROTO]. `x_c = 0.25` matches GU24's 0.4–0.5 W of compute in Instagram. `x_c = 0.8` sits just above CHI23's ≈ 3 W
  sustainable SoC power.
* Leakage growth with temperature is left out of the baseline (structural variant X13).

### 3.4 Radios: Wi-Fi, LTE, 5G NSA

```
P_net = w·P_idle^W + 1[c ≠ none]·P_idle^c + Φ_ρ                                                        (E14)
Φ_ρ(κ = burst) = m_ρ·[ (P_tail − P_idle)(1 − e^(−λT_tail)) + (P_act − P_tail)(1 − e^(−λT_act))
                       + λ·e^(−λT_tail)·E_pro + e_b·(d_dl + κ_ul·d_ul) ]
Φ_ρ(κ = rt)    = m_ρ·P_rt
m_ρ = 1 for Wi-Fi;   m_ρ = m_sig(RSRP) = min{m_max, 1 + γ_sig·max(0, (RSRP_ref − RSRP)/10 dB)} for cellular
```

**Derivation (Appendix A.3).** Bursts arrive as a Poisson process of rate λ and are short compared with the timers.
The radio is connected-active if and only if a burst arrived within the C-DRX inactivity timer `T_act`. It is in
the connected tail (DRX) state if the last burst arrived between `T_act` and `T_tail` ago. Hence
`Pr(active) = 1 − e^(−λT_act)` and `Pr(active or tail) = 1 − e^(−λT_tail)`, exactly. A promotion is paid only when
a burst finds the radio idle, at rate `λe^(−λT_tail)`. The per-bit term covers the extra energy of transfers.

* **Why `P_act` was added (fix to PF).** Without it, a game exchanging 20 packets/s over LTE costs the same radio
  power as over Wi-Fi (0.13 W). With it, GM-L costs 0.55 W [PROTO], consistent with GU24's ≈ 0.73 W cellular rail
  in a game on an Exynos-modem Pixel 7a.
* **Real-time class.** Voice and video calls keep the modem in continuous Tx/Rx; their measured connected-active
  power is used directly.
* **Weak signal.** `m_sig` applies to the cellular traffic, tail and promotion terms, not to idle paging (no
  source). `γ_sig = 0.26` gives `m_sig(−105) = 1.52`, which reproduces DI13's +52 % energy from −85 to −105 dBm.
* **Optional Wi-Fi cliff** (scenario S11b): `m_W(RSSI) = 1 + 1.13·exp((−85 − RSSI)/2.54 dB)` [DERIVED from DI13:
  +113 % at −85 dBm and +810 % at −90 dBm relative to −50 dBm].
* **Check values [PROTO].** Radio power (idle of the associated radios included): video at 3 Mbit/s with
  λ = 0.25 s⁻¹ costs 0.074 W over Wi-Fi, 0.220 W over LTE and 0.363 W over 5G NSA. The Wi-Fi traffic term alone is
  0.055 W; Apple's streamed − local video difference (both Wi-Fi associated) is 0.06–0.08 W. The LTE − Wi-Fi gap of
  0.146 W lies in GU24's 0.12–0.24 W.

### 3.5 GNSS, camera, audio, Bluetooth, sensors

```
P_GNSS = g·P_G;   P_cam = c_r·ΔP_cam + c_f·ΔP_camf;   P_aud = au·ΔP_aud + bt·ΔP_BT;   P_sens per primitive   (E15)
```

* Optional cold-fix impulse: 154 J at GNSS start (CH13), which is 0.24 % of `E_u`.
* `ΔP_cam` is the rear pipeline above the awake floor and display (sensor, ISP, encoder, CPU). It is [DERIVED]
  from CX13 (Galaxy S2 recording 1.40–1.65 W total, AMOLED 0.17–0.18 W, ≈ 0.2 W floor), giving 1.03–1.28 W. CH13's
  S3 totals (2.26–2.61 W) are the upper end. Plausibility: Bonn camera increments are 1.22 W (14 Pro) and 1.04 W
  (13 Pro) (validation_data a.2).
* Sensors at 50 Hz: accelerometer 5 mW, gyroscope 30 mW, magnetometer 12 mW (CH13). Navigation and tracking use
  `P_sens = 0.047 W`.
* Omitted, with justification: haptics (no source; tiny duty cycle), flash storage (< 8 mW, CH13), loudspeaker at
  high volume (no source), BLE scanning (no source), self-discharge (A15), PMIC light-load efficiency drop
  (unverified, absorbed in the `P_susp` range).

### 3.6 SoC power under throttling and the performance ratio

The delivered SoC throughput under throttling is `x_eff = P_SoC,dem⁻¹(P_SoC)`, found by a scalar root of E13. The
performance ratio is `ρ_perf = x_eff/x_c`, which is a frequency ratio, since `f ∝ x`.

### 3.7 Burst-envelope amplitudes [DERIVED / ASSUMPTION]

| Activity class | ΔP_b (W) | Basis |
|---|---|---|
| Standby, audio, tracking (radio wake-ups) | 1.0 | LTE promotion 1.21 W (HU12); 4G→5G switch 0.70–1.49 W (NA21) |
| Interactive UI, navigation, calls | 3.0 | LEE20 peak currents 1–2 A ≈ 3.7–7.4 W total; CH13 CPU + GPU peak rails 4.3 W |
| Gaming, camera | 4.0 | CH13 device peak > 7 W |
| Video playback | 2.0 | [ASSUMPTION] decoder and chunk-download peaks |
| Explicit long burst `BU` (S9) | 0 | its power is part of the mean load `u(t)` |

All amplitudes are multiplied by a common factor `m_b ~ U(0.5, 2.0)` in UQ. This factor dominates the cold, aged
brown-out (§10). It is also the lever of the OS peak-power governor (O1, §12.2).

---

## 4. Parameter tables (reference phone RP-63)

**RP-63** is a generic 2024-class 6.3-inch phone: 19.5:9 LTPO OLED (active area `A = 2.455·d² cm²` for diagonal
`d` in inches, so 97.4 cm²), a 5G-capable modem, and a single LCO/graphite pouch of 4.70 Ah, `V_max` 4.45 V, mean
OCV 3.870 V. It is **not** a model of any specific product. Derived energies [PROTO]: `E_n = 18.19 Wh` (OCV integral
over 0 ≤ z ≤ 1) and `E_u = 17.79 Wh` (from `z_e` to 1 at vanishing load).

**Distributions.** `U(a,b)` uniform; `T(a,m,b)` triangular with mode `m` = default; `LU(a,b)` log-uniform (used when
`b/a ≥ 3`). All factors are independent except `V_sys = V_cut − Δ_sys`. "Fixed" means the quantity is a
definition or normalisation (e.g. `APL_ref`, `RSRP_ref`, `c_ref`), a fitted shape constant whose uncertainty is
carried by a companion factor (OCV coefficients by `δU`; `z0` by `A_R`), or a term too small to matter (Wi-Fi
tail and promotion energies: < 0.01 W at any burst rate used here). Model-form choices behind fixed constants are
covered by the structural variants of §10.4. Confidence (Conf.): H primary/official, M figure read-off or secondary
or model output, L our estimate.

### 4.1 Cell and pack

| Symbol | Value | Unit | UQ | Source / justification | Conf. |
|---|---|---|---|---|---|
| `Q_n` | 4.70 | Ah | T(4.56, 4.70, 4.80) | PIX9: Pixel 9 4.700 typical / 4.558 minimum; iPhone 16 Plus/Pro Max 4.67–4.69 (MR24) | H |
| `b0, b1, a0, a1, a2, a3` | 0.60117, 36.426, 3.6193, 0.44311, −0.48605, 0.63127 | V, –, V, V, V, V | fixed | [RUN] CRM06 form fitted to PYB-AI20 OCV (RMSE 7.0 mV) | M |
| `ΔV_hv, p` | 0.242, 3.86 | V, – | fixed; X18 uses (0.192, 5.44) and (0.272, 2.90) | cell_ecm a.2 construction for 4.45 V / 3.87 V (4.40/3.85 and 4.48/3.89 variants) | L |
| `δU` | 0 | mV | U(−20, 20) | hysteresis and cell-to-cell offset; OVEJ19 mechanism, magnitude unverified | L |
| `k_S` (multiplier on `∂U/∂T`) | 1 | – | U(0.5, 1.5) | [RUN] Ai2020 entropic fits; extrapolated below z ≈ 0.45 | M–L |
| `R0,ref` (pack: cell + PCM + flex, z = 0.5, 25 °C) | 0.060 | Ω | LU(0.035, 0.130) | cell_ecm a.3: PINE ≤ 100 mΩ at 1 kHz incl. PCM, PCM 20–60 mΩ, BU802a legacy 130 mΩ, [RUN] cell 19 mΩ. Plausibility: under `R·Q` scaling, `R0 + R1` = 0.12 Ω for a 2.7 Ah pack; HE19 measured 0.1 Ω (Galaxy S3), LEE20 0.13–0.33 Ω | L |
| `A_R` | 0.80 | – | U(0.5, 2.0) | [RUN] Ai2020 R0(0.1)/R0(0.5) = 1.39 | M |
| `z0` | 0.139 | – | fixed | [RUN] | M |
| `E_a,R0` | 18 | kJ/mol | T(12, 18, 27) | [RUN] 18.0; BARC22 15.6; SUN23 15.8; LEE20-derived 12–27 (thermal_aging C6) | M |
| `R1,ref` | 0.010 | Ω | LU(0.003, 0.025) | [RUN] Ai2020 14.4 mΩ·(2.28/4.7) ≈ 7; real cells 5–25 (cell_ecm a.3) | M/L |
| `τ1,ref` | 30 | s | U(10, 60) | [RUN] 29–35 s | M |
| `E_a,R1` | 30 | kJ/mol | T(25, 30, 55) | [RUN] 27.7–31.4; BARC22 29.2 (55 exponential part) | M |
| `E_a,τ` | 21 | kJ/mol | U(15, 27) | [RUN] Marquis2019 τ1 29 → 131 s from 25 to −20 °C | M |
| `V_cut` | 3.40 | V | U(3.30, 3.50) | QC-FG example 3400 mV ("FG will steer SOC to zero"); LEE20 "≈ 3.4 V" | M/H |
| `Δ_sys` (`V_sys = V_cut − Δ_sys`) | 0.10 | V | U(0, 0.20) | QC-SMB example 3300 mV "minimum allowed system voltage"; its use as brown-out level is an [ASSUMPTION] | M/L |
| `T_shut` (battery) | 68 | °C | fixed | KANG19 Android shutdown example | M |
| `η_c` (charge coulombic efficiency) | 0.99 | – | U(0.98, 1.00) | PLETT-N | M |

### 4.2 Thermal

| Symbol | Value | Unit | UQ | Source / justification | Conf. |
|---|---|---|---|---|---|
| `C_th` | 150 | J/K | T(100, 150, 230) | thermal_aging T16–T22: m·c_p (MAL99), τ·R (HE19), energy/ΔT (KANG19) agree | M |
| `R_th` | 3.0 | K/W | T(2, 3, 4) | HE19 Fig. 10 [DERIVED]; convection + radiation estimate 2.5–3.6 K/W (T26) | M/L |
| `R_skin` | 5.5 | K/W | T(4, 5.5, 9) | PU23 Pixel 5: 5 K/W; KANG19 5–7 K/W | M |
| `T_lim` | 45 | °C | U(41, 48) | CHI15 design limit; KANG19 pain at 42–45 °C; IEC 62368-1 TS1 48 °C (UL-TS) | H |
| `ΔT_thr` | 2 | K | fixed | [ASSUMPTION] throttling starts just before the limit | L |
| `w_thr` | 1.0 | K | U(0.5, 3) | [ASSUMPTION] graded AOSP-TH levels. Plausibility: capped gaming 2.9–3.1 W at 25 °C [PROTO], inside CHI15's 2.2–3.8 W | L |
| `P_fl` | 0.5 | W | U(0.3, 1.0) | [ASSUMPTION] minimum SoC power when throttled | L |
| `k_env` | 1.0 | – | U(0.8, 1.5) | [ASSUMPTION]; KANG19: a case slowed the rise to 45 °C by 54 % | L |
| `ΔT_ch` | 8 | K | U(5, 12) | [ASSUMPTION]; KANG19 H9: ≈ 37 °C surface while fast-charging an idle phone | L |
| `η_charger` | 0.90 | – | U(0.85, 0.95) | [ASSUMPTION] (thermal_aging §2.1) | L |

### 4.3 Device power (battery side)

| Symbol | Value | Unit | UQ | Source / justification | Conf. |
|---|---|---|---|---|---|
| `P_susp` | 0.020 | W | U(0.010, 0.030) | CH13 suspend 24–27 mW (2008–12 phones); CH15 SoC suspend 3.8–5.1 mA; HU12 base 11.4 mW | H/M |
| `ΔP_awake` | **0.090** [CAL] | W | LU(0.05, 0.25) (prior) | Calibrated once, so that the Bonn iPhone 14 Pro black-screen idle TTE is 55.0 h with the other floor terms at their defaults (exact value 0.0877 W [PROTO], V1.1). Prior: CH10 Freerunner awake − suspend = 200 mW; AOSP cpu.awake 50 mA (illustrative). Blind prior 0.10 W gives −5.1 % on the same run | M/L |
| `E_wake` | 0.10 | J | LU(0.05, 0.30) | [DERIVED] ≈ ΔP_awake × 1 s wake ([ASSUMPTION] duration) | L |
| `P_d0` | 0.12 | W | U(0.08, 0.25) | CH13 S3 constant term 223 mW (2012 panel, upper); GU24 display rail ≈ 0.30 W total (Pixel 7a) | M |
| `p_c` (content power density at b = 0.5, APL_ref) | 2.567 (0.25 W at 97.4 cm²) | mW/cm² | U(1.23, 5.13) (0.12–0.50 W) | GU24 0.30 / 0.45 W total display (6.1" / 6.6"); CH13 0.3–0.6 W at 50 % on 4.8" | M |
| `APL_ref` | 0.70 | – | fixed (normalisation) | DZ11: 65 % white pixels on mobile web | H |
| `γ_b` | 2.0 | – | U(1.0, 2.9) | CH13 near-affine slider (≈ 1); DH21-derived 2.9 (PF) | M/L |
| `κ_RR` | 2.8 | mW/Hz | U(1.0, 3.3) | AT20 160–180 mW from 60 to 120 Hz on a fixed-rate panel; LTPO lowers it (unverified) | M/L |
| `P_SoC,max` | 5.0 | W | U(3, 7) | CH13 CPU 2.85 + GPU 1.42 W rails × 1.2 | M |
| `v0` | 0.446 | – | U(0.35, 0.55) | LX MT8183 OPP table [DERIVED] | M |
| `P_G` | 0.27 | W | T(0.14, 0.27, 0.45) | KW19 0.23 / 0.32 J per 1 Hz fix; CH13 0.39–0.43 W incl. app; CH10 0.14–0.17 W module | H |
| `ΔP_cam` | 1.15 | W | U(0.9, 2.0) | CX13 [DERIVED] 1.03–1.28 W; CH13 upper | M |
| `ΔP_camf` | 0.30 | W | U(0.2, 0.5) | GU24 camera rail ≈ 0.26 W × 1.2 | M |
| `ΔP_aud` | 0.05 | W | U(0.02, 0.15) | CH13 audio 226 mW total, volume +5 mW; HE20 0.11 W incl. speaker | M |
| `ΔP_BT` | 0.04 | W | U(0.03, 0.05) | CH10 +36–45 mW | H |
| `P_sens` (nav/tracking) | 0.047 | W | fixed | CH13 accelerometer 5 + gyroscope 30 + magnetometer 12 mW | H |
| `m_b` (burst multiplier) | 1 | – | U(0.5, 2.0) | §3.7 | L |

### 4.4 Radio

| Parameter | Wi-Fi | LTE | 5G NSA low band | Sources |
|---|---|---|---|---|
| `P_idle` (W) | 0.004, U(0.002, 0.010) | 0.015, U(0.005, 0.035) | 0.030, U(0.015, 0.060) | HU12 Wi-Fi beacon 1.6 mW, LTE DRX 19.7 mW [DERIVED]; CH15 1.1 mA Wi-Fi, 2.3–8.3 mA paging; GU24 5G ≈ 2× 4G standby |
| `P_tail` (W) | 0.12, fixed | 0.12, U(0.066, 0.178) | 0.255, U(0.249, 0.593) | HU12 Wi-Fi 119 mW; NA21 Table 2 (4G 66–178, NSA-LB 249–260, SA 593 mW) |
| `T_tail` (s) | 0.24, fixed | 10.0, U(5.0, 11.6) | 10.3, U(10.2, 21.4) | HU12; NA21 Table 7; XU20 chained NSA tail 21.4 s |
| `P_act` (W) | 0.13, U(0.10, 0.30) | 0.60, LU(0.30, 1.30) | 0.90, LU(0.50, 1.80) | HU12 Wi-Fi β 133 mW; LTE: GU24 Pixel 7a cellular rail 0.73–0.95 W (Exynos, upper), S23+ ≤ 0.60 W, HU12 β 1.29 W (2012, upper) [DERIVED]; NR: NA21 low-band intercept ≈ +0.5 W over 4G (Fig. 11 read-off), GU24 5G − 4G +0.07–0.15 W, XU20 "2–3×" |
| `T_act` (s) | 0.24, fixed | 0.10, U(0.08, 0.20) | 0.10, U(0.08, 0.20) | XU20 inactivity timers 80 / 100 ms |
| `E_pro` (J) | 0.0098, fixed | 0.31, U(0.1, 0.5) | 1.0, U(0.5, 2.0) | HU12 124.4 mW × 79.1 ms; HU12 1210.7 mW × 260 ms; NA21 4G→5G switch 0.70–0.80 W × 1.44 s [DERIVED] |
| `e_b` (W per Mbit/s) | 0.015, LU(0.005, 0.06) | 0.030, LU(0.014, 0.08) | 0.030, LU(0.014, 0.08) | SU14 14–60 nJ/bit (median 60); NA21 4G DL slope ≈ 14 mW/Mbps; HU12 LTE α_d 52 mW/Mbps; Apple streamed − local 0.06–0.08 W [DERIVED] |
| `κ_ul` (UL/DL energy per bit) | 5, U(2, 8.4) | same | same | NA21 68/14 ≈ 4.9; HU12 438/52 = 8.4 |
| `P_rt` (W) | 0.30, U(0.13, 0.8) | 1.20, U(0.6, 2.6) | 1.40, U(1.2, 3.0) | device_power §0 Wi-Fi active 0.3 W (SU14 0.24–0.80 W); HU12 LTE β 1.29 W; GU24 Meet cellular rail 2.58 W (upper); XU20 5G 2–3× 4G (2019 modem) |
| `γ_sig`, `RSRP_ref`, `m_max` | – | 0.26 per 10 dB, LU(0.10, 1.0); −85 dBm; 3.5 | same | DI13: +52 % energy −85 → −105 dBm [DERIVED]; Tx/Rx current ratios 1.2–3.5× |

### 4.5 Aging

| Symbol | Value | Unit | UQ | Source / justification | Conf. |
|---|---|---|---|---|---|
| `k_c` | 1.413e-3 | d^−½ | U(1.05e-3, 2.09e-3) | so that `q_cal(365 d, 25 °C, g = 1) = 2.7 %` (range 2–4 %): SCH14 2.4 % at 3.7 V (G6), KEIL16 2–3 % after 10 months at 50 % SOC | M |
| `E_cal` | 58 | kJ/mol | U(20, 60) | SCH14 via CA21 (6976 K × R_g); SCM18 LFP 20.6 | M |
| `g(z)` step levels | 0.7 / 1.0 / 1.75 | – | high level U(1.3, 1.9) | KEIL16 NCA at 25 °C (read-off; step near 60 % SOC) | M |
| `κ_HV` | 1.0 | – | U(0, 1.5) | SCH14 voltage law extrapolated above 4.2 V; NEL18 qualitative | L |
| `k_H` | 2.697e-3 | FCE^−½ | follows calibration | calibrated (§11.1): SOH = 0.80 at 1000 FCE (AP-BATT iPhone 15+, GOOG Pixel 8a+) with plating share `w_L` | M/L |
| `k_L` | 1.288e-3 | FCE^−½ | follows calibration | same calibration | M/L |
| `w_L` (plating share of cycle loss at 25 °C, 1C) | 0.3 | – | U(0.1, 0.5) | [ASSUMPTION]; WAL14/WU17: cycle aging is mildest near 25 °C, with plating below and SEI above | L |
| `E_H` | 33 | kJ/mol | U(25, 45) | SCM18 32.7 (LFP; structure transferred) | M |
| `E_L` | 55.5 | kJ/mol | U(40, 70) | SCM18 55.5 (inverse Arrhenius) | M |
| `β_L` | 2.64 | h | U(0.5, 2.64) | SCM18 (LFP). Only the upper end is sourced; informal phone fast-charge tests suggest a smaller penalty (thermal_aging H12, L) | M/L |
| `c_ref` | 1.0 | h⁻¹ | fixed | SCM18 1C reference | M |
| `γ_R` | 2.5 | – | LU(1.5, 7.0) | SCH14 G8 (+41–65 % at 20 % fade, γ ≈ 2.0–3.3); KEIL16 NMC +48 % at 23 % (2.1), NCA +83 % at 12.5 % (6.6) | L |
| `N_80` | 1000 | FCE | {500, 800, 1000} | AP-BATT (≤ iPhone 14: 500; iPhone 15+: 1000), GOOG (Pixel 3–8 Pro: 800; 8a+: 1000) | H |
| Self-discharge | 0.1 | %/month | 0.05–2 | ZIL19 (storage scenarios only) | M/L |

### 4.6 Usage model (CTMC, §5.3)

| Symbol | Typical | Heavy | UQ | Source / justification |
|---|---|---|---|---|
| Session start rate `λ_on` | 4 h⁻¹ | 4 h⁻¹ | × LU(0.7, 1.3) | [ASSUMPTION] (from PS). Plausibility: screen-on 6.5 h/day (typical) and 11.5 h/day (heavy) [PROTO]; typical `P̄` = 0.440 W is consistent with Pixel 8a's "at least 31 h" median-user claim (`P̄ ≤ 0.56 W`, GPX) |
| Activity probabilities `p_j` (social, video, game, nav, vcall, camera, audio) | 0.55, 0.12, 0.05, 0.04, 0.04, 0.08, 0.12 | 0.45, 0.20, 0.15, 0.05, 0.05, 0.05, 0.05 | Dirichlet(50·p) | [ASSUMPTION] |
| Mean durations `d_j` (min) | 5, 15, 20, 20, 10, 2, 30 | 8, 25, 30, 25, 15, 3, 30 | each × LU(0.7, 1.3) | [ASSUMPTION] |
| Connectivity Wi-Fi ⇄ LTE, mean sojourns | 2.5 h / 2.5 h | same | × LU(0.5, 2) | [ASSUMPTION]; 5G variant swaps LTE for NR |
| Per-visit jitter of `b`, `APL`, `λ`, `d_dl` | ±30 % | ±30 % | U(0.7, 1.3) | [ASSUMPTION]; only in the usage-fluctuation analysis (§10.5) |

### 4.7 Validation device configurations (same physics, no re-tuning)

Each device uses the RP-63 parameters except those listed. Pack resistances are scaled by `R·Q = const` (A12):
`R0,ref, R1,ref ← (4.70/Q)·(RP-63 values)`. The display area follows the diagonal at 19.5:9. Brightness for
"150 nit" tests is `b = (150/L_max)^(1/γ_b)` with `L_max = 1000 cd/m²` (iPhone 14 Pro typical maximum, APL14P;
an [ASSUMPTION] for the other phones).

| Device | Q (Ah) | Diagonal | Source of Q | Notes |
|---|---|---|---|---|
| iPhone 14 Pro (Bonn) | 3.32 | 6.1" | SCH26 measured initial 3321–3326 mAh | airplane mode; OLED on and black (`APL = 0`); `E_u` = 12.565 Wh [PROTO] |
| iPhone 13 Pro (Bonn) | 3.16 | 6.1" | SCH26 measured 3142–3181 mAh (pre-used) | same |
| iPhone 14 Pro (Apple, LMAG22) | 3.32 | 6.1" | as above | Apple tests: Wi-Fi associated, BT headphones |
| iPhone 14 Pro Max | 4.323 | 6.7" | N9M22 (16.68 Wh) | NBC22 |
| Pixel 8 | 4.575 | 6.2" | GPX typical | NBC23 |
| Galaxy S24 (Exynos) | 4.000 | 6.2" | NBC24 typical | NBC24 |
| Pixel 8a | 4.492 | 6.1" | GPX typical | median-user claim |
| LG HG2 cell | 3.0 (C/20: 2.78 measured) | – | KOL20 | cell only: no PCM, no device loads. CAL at 25 °C only: NMC OCV from the C/20 curve, `R0 + R1 = 0.038 Ω` split 6:1, cell `C_th`, `R_th` from the 25 °C 2C temperature rise. All other temperatures predicted with RP-63 activation energies |
| NASA B0005 cell | 1.86 (fresh) | – | NASA07 | measured fade and EIS `R_e`, `R_ct` growth inserted directly |

**[REVIEW 2026-11-17] Cell-only devices and the heat balance.** Assumption A4 ("all device power ends up as heat
inside the phone") holds for a phone, not for a bare cell on a tester, whose load is external. For the HG2 and NASA
cells the term `P` is therefore dropped from E7: `C_th·dT/dt = I²R0 + v1²/R1 − I·T·k_S·∂U/∂T − (T − T_amb)/(k_env·R_th)`
(`Params.load_heat = False`). This is also the heat definition used by the HG2 thermal [CAL] (`I·(U − V)`). With the
phone form of E7, the calibrated HG2 cell at 2C (20.9 W delivered) would heat to the 68 °C shutdown in 0.31 h,
whereas the data show a 3.7 K rise over the 0.45 h discharge; with `P` dropped the model gives 0.455 h and +4.7 K.

---

## 5. Usage: activity primitives, scenarios and CTMC profiles

### 5.1 Activity primitives (constant `u`)

`—` means 0 or not applicable. `λ` in s⁻¹; throughputs in Mbit/s. "P [PROTO]" is the battery-side mean power on
RP-63 at 25 °C before throttling.

| ID | s_scr | a | b | APL | f_r | x_c | w | c | ρ | κ | λ | d_dl / d_ul | other | ΔP_b (W) | P [PROTO] (W) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SB-W standby, home Wi-Fi | 0 | 0 | — | — | — | 0 | 1 | LTE | wifi | burst | 1/120 | 0 | — | 1.0 | 0.040 |
| SB-L standby, cellular only | 0 | 0 | — | — | — | 0 | 0 | LTE | cell | burst | 1/120 | 0 | — | 1.0 | 0.047 |
| SB-5 standby, 5G NSA, heavy push | 0 | 0 | — | — | — | 0 | 0 | NR | cell | burst | 1/60 | 0 | — | 1.0 | 0.102 |
| RD reading / web | 1 | 1 | 0.4 | 0.65 | 60 | 0.05 | 1 | LTE | wifi | burst | 1/30 | 0.05 / 0 | — | 3.0 | 0.456 |
| SO social feed, Wi-Fi | 1 | 1 | 0.4 | 0.55 | 120 | 0.25 | 1 | LTE | wifi | burst | 2.0 | 2.0 / 0 | — | 3.0 | 1.060 |
| SO-L social feed, LTE | 1 | 1 | 0.4 | 0.55 | 120 | 0.25 | 0 | LTE | cell | burst | 2.0 | 2.0 / 0 | — | 3.0 | 1.218 |
| VW video stream, Wi-Fi | 1 | 1 | 0.4 | 0.30 | 60 | 0.10 | 1 | LTE | wifi | burst | 0.25 | 3.0 / 0 | au | 2.0 | 0.548 |
| VL video stream, LTE | 1 | 1 | 0.4 | 0.30 | 60 | 0.10 | 0 | LTE | cell | burst | 0.25 | 3.0 / 0 | au | 2.0 | 0.694 |
| V5 video stream, 5G NSA | 1 | 1 | 0.4 | 0.30 | 60 | 0.10 | 0 | NR | cell | burst | 0.25 | 3.0 / 0 | au | 2.0 | 0.837 |
| NV navigation | 1 | 1 | 0.7 | 0.55 | 60 | 0.25 | 0 | LTE | cell | burst | 0.2 | 0.1 / 0 | g = 1, P_sens | 3.0 | 1.486 |
| GM 3D game, Wi-Fi | 1 | 1 | 0.5 | 0.40 | 120 | 0.80 | 1 | LTE | wifi | burst | 20 | 0.2 / 0.05 | au | 4.0 | 3.906 (capped ≈ 3.1) |
| GM-L 3D game, LTE | 1 | 1 | 0.5 | 0.40 | 120 | 0.80 | 0 | LTE | cell | burst | 20 | 0.2 / 0.05 | au | 4.0 | 4.302 |
| VC video call, LTE | 1 | 1 | 0.5 | 0.45 | 60 | 0.30 | 0 | LTE | cell | rt | — | — | c_f, au | 3.0 | 2.518 |
| VC-W video call, Wi-Fi | 1 | 1 | 0.5 | 0.45 | 60 | 0.30 | 1 | LTE | wifi | rt | — | — | c_f, au | 3.0 | 1.622 |
| CM rear camera | 1 | 1 | 0.6 | 0.45 | 60 | 0.05 | 1 | LTE | wifi | burst | 0 | — | c_r | 4.0 | 1.687 |
| AU audio, screen off, BT, Wi-Fi | 0 | 1 | — | — | — | 0 | 1 | LTE | wifi | burst | 1/30 | 0.256 / 0 | au, bt | 1.0 | 0.224 |
| AU-L audio, screen off, BT, LTE | 0 | 1 | — | — | — | 0 | 0 | LTE | cell | burst | 1/30 | 0.256 / 0 | au, bt | 1.0 | 0.261 |
| TRK screen-off GNSS tracking | 0 | 1 | — | — | — | 0 | 0 | LTE | cell | burst | 1/60 | 0.01 / 0 | g = 1, P_sens | 1.0 | 0.464 |
| BU explicit burst: camera + HDR + upload | 1 | 1 | 1.0 | 0.45 | 60 | 0.50 | 0 | LTE | cell | burst | 1.0 | 0 / 4.0 | c_r | 0 | 4.095 |

Labelled assumptions: every `λ`, throughput, APL, `x_c` and brightness above is an [ASSUMPTION] with the ±30 %
jitter of §4.6 in usage-fluctuation runs. Anchors: social `x_c` (GU24 compute 0.4–0.5 W); video APL (dark film
content); navigation brightness 0.7 (outdoor, above DH21's 30–40 % indoor auto-brightness); reading APL 0.65
(DZ11); streaming 3 Mbit/s (a 720p-class stream).

### 5.2 Scenario definitions

Common settings unless stated: RP-63, rested start (`v1 = 0`), `T0 = T_amb = 25 °C`, SOH = 1, RSRP = −85 dBm,
`k_env = 1`, `s0 = 100 %`. A schedule is a list of `(duration, primitive)` that repeats until the TTE event. Every
constant scenario is also run from `s0 ∈ {50, 20} %` (S15).

| ID | Name | u(t) | Environment / state |
|---|---|---|---|
| S1 | Standby at home | SB-W, constant | – |
| S2a / S2b | Standby on cellular / 5G NSA with heavy push | SB-L / SB-5, constant | – |
| S3 | Reading | RD, constant | – |
| S4 | Social feed on Wi-Fi | SO, constant | – |
| S5a / S5b / S5c | Video streaming over Wi-Fi / LTE / 5G NSA | VW / VL / V5, constant | – |
| S6 | Navigation | NV, constant | – |
| S7 | Sustained 3D gaming | GM, constant | `T_amb ∈ {15, 25, 35} °C`; variant with throttling disabled |
| S8a / S8b | Video call over LTE / Wi-Fi | VC / VC-W, constant | – |
| S9 | Cold outdoor hike with photos | periodic 600 s: TRK 570 s → BU 30 s | `T_amb ∈ {25, −10, −20} °C`, `T0 = 20 °C` (carried out from indoors); SOH ∈ {1.0, 0.8} |
| S10 | Aged battery | S4 at `T_amb ∈ {25, 10, 0, −10, −20} °C` (S10a); S12 at 25 and 0 °C (S10b) | `q_loss ∈ {0, 0.1, 0.2}`, `r_growth = γ_R·q_loss` |
| S11 | Weak cellular signal | SO-L (S11a) and VC (S11c) at RSRP ∈ {−85, −95, −105, −110, −120} dBm; optional Wi-Fi cliff S11b (SO at RSSI −85 / −90 dBm) | – |
| S12 | Scheduled full day (9.5 h screen-on) | from 07:00, minutes: RD 30 · NV 45 · 7 × (SB-W 26 + SO 4) · VW 30 · SO 30 · 8 × (SB-W 26 + SO 4) · NV 45 · VC 30 · GM 60 · VW 120 · SO 60 · RD 60 · SB-W 480 (23:00–07:00); repeated | – |
| S13 | Typical day, stochastic | CTMC "typical" (§5.3), start in (idle, Wi-Fi) | 25 °C |
| S14 | Heavy day, stochastic | CTMC "heavy" | 25 °C |
| S15 | Initial-SOC sweep | S4, S5a, S7, S12, S13 from `s0 ∈ {100, 80, 50, 20, 10} %` | – |
| S16 | Audio, screen off | AU, constant | – |

### 5.3 CTMC usage model

**Activity chain (hub and spoke).** From `idle`, a session starts at rate `λ_on`. It is activity `j` with
probability `p_j`, lasts Exp(mean `d_j`), and returns to idle. The generator is `Λ_A[idle, j] = λ_on·p_j`,
`Λ_A[j, idle] = 1/d_j`.

**Connectivity chain.** Wi-Fi ⇄ LTE with mean sojourns 2.5 h each. The product generator is the Kronecker sum
`Λ = Λ_A ⊗ I₂ + I₈ ⊗ Λ_C` (16 states).

**State → primitive map.**

| Activity | Wi-Fi state | LTE state |
|---|---|---|
| idle | SB-W | SB-L |
| social | SO | SO-L |
| video | VW | VL |
| game | GM | GM-L |
| nav | NV | NV |
| vcall | VC-W | VC |
| camera | CM | CM |
| audio | AU | AU-L |

**State powers for the analytic formulas.** These are the quasi-steady thermal fixed points of F5 at `T_amb`. Only
the game states are throttled; the others equal their demand. In the L2 Monte Carlo the throttle evolves
dynamically.

### 5.4 Baseline results [PROTO]

RP-63, fresh, 25 °C, from 100 %, L2 model with LSODA (§6.2). `s_off` is in D1 percent.

| Scenario | Mean P (W) | TTE | Event | s_off | T_end / T_s,end (°C) |
|---|---|---|---|---|---|
| S1 standby, home | 0.040 | 442.7 h (18.4 d) | cutoff | 0.02 % | 25.1 / 25.2 |
| S2a standby, cellular | 0.047 | 378.3 h | cutoff | 0.02 % | 25.1 / 25.3 |
| S2b standby, 5G heavy push | 0.102 | 173.8 h | cutoff | 0.04 % | 25.3 / 25.6 |
| S3 reading | 0.456 | 38.88 h | cutoff | 0.18 % | 26.4 / 27.7 |
| S4 social, Wi-Fi | 1.060 | 16.64 h | cutoff | 0.40 % | 28.4 / 31.2 |
| S5a video, Wi-Fi | 0.548 | 32.32 h | cutoff | 0.21 % | 26.7 / 28.2 |
| S5b video, LTE | 0.694 | 25.49 h | cutoff | 0.27 % | 27.2 / 29.0 |
| S5c video, 5G NSA | 0.837 | 21.10 h | cutoff | 0.32 % | 27.7 / 29.9 |
| S6 navigation | 1.486 | 11.83 h | cutoff | 0.56 % | 29.8 / 33.7 |
| S7 gaming, 25 °C | 3.906 → 2.94 (throttled) | 5.68 h | cutoff | 1.07 % | 34.5 / 42.4 |
| S7 gaming, 35 °C | → 1.65 | 10.29 h | cutoff | 0.57 % | 40.3 / 44.7 |
| S7 gaming, 15 °C | → 3.89 | 4.40 h | brownout | 1.89 % | 27.6 / 38.2 |
| S7 gaming, 25 °C, no throttling | 3.906 | 4.43 h | – | – | – |
| S8a video call, LTE | 2.518 | 6.93 h | cutoff | 0.93 % | 33.1 / 39.9 |
| S8b video call, Wi-Fi | 1.622 | 10.83 h | cutoff | 0.61 % | 30.2 / 34.5 |
| CM camera (constant) | 1.687 | 10.40 h | brownout | 0.68 % | 30.4 / 34.9 |
| S12 scheduled day, from 100 / 50 / 20 % | 0.528 (demand mean) | 34.92 / 11.79 / 8.70 h | cutoff | 0.93 % (from 100 %) | – |
| S16 audio, screen off | 0.224 | 79.20 h | cutoff | 0.09 % | 25.7 / 26.3 |

| S9 cold hike (TRK 0.464 W + BU 4.095 W bursts) | 25 °C | −10 °C | −20 °C |
|---|---|---|---|
| SOH 1.0: TTE, s_off | 26.99 h, 1.7 % | 25.33 h, 7.2 % | 22.67 h, 16.9 % |
| SOH 0.8: TTE, s_off | 21.32 h, 2.7 % | 17.00 h, 22.2 % | 12.33 h, **44.2 %** |

All S9 events are `cutoff` during a BU burst. Cold-soaking the phone at start (`T0 = −10 °C`) does not change the
−10 °C aged result (17.00 h): the thermal transient (7.5 min) is short compared with the TTE.

**Stranded-charge map (S10a social with envelope bursts `ΔP_b = 3 W`)** [PROTO]: `s_off` and, in parentheses,
TTE. The event is `brownout` in every cell except 25 °C / SOH 1.0 (`cutoff`).

| T_amb | SOH 1.0 | SOH 0.9 (R × 1.25) | SOH 0.8 (R × 1.5) |
|---|---|---|---|
| 25 °C | 0.4 % (16.64 h) | 0.5 % (14.94 h) | 0.9 % (13.22 h) |
| 10 °C | 0.7 % (16.57 h) | 1.3 % (14.80 h) | 2.1 % (13.05 h) |
| 0 °C | 1.4 % (16.43 h) | 2.4 % (14.62 h) | 3.8 % (12.80 h) |
| −10 °C | 2.8 % (16.18 h) | 4.9 % (14.22 h) | 8.5 % (12.15 h) |
| −20 °C | 6.4 % (15.54 h) | 13.1 % (12.99 h) | 23.3 % (10.23 h) |

**Other [PROTO] results.**

* S10b scheduled day: 25 °C SOH 0.8 gives 29.19 h (−16.4 %, cutoff). At 0 °C: SOH 1.0 gives 34.50 h (brownout,
  `s_off` 1.9 %), SOH 0.8 gives 24.78 h (−29.0 %, brownout, `s_off` 4.7 %).
* S11: SO-L at −85 / −95 / −105 / −110 / −120 dBm lasts 14.46 / 13.72 / 13.05 / 12.73 / 12.15 h (−12.0 % at
  −110 dBm). VC lasts 6.93 / 6.17 / 5.58 / 5.33 / 4.87 h (−23 % at −110 dBm).
* Navigation at −10 °C, fresh pack, cold-soaked: 11.43 h vs 11.83 h at 25 °C (−3.4 %). Self-heating keeps
  `T ≈ T_amb + 4.5 K`, and `R0` is small.

**What the baseline shows.**

1. **25 °C, fresh: an energy problem.** L2 TTE is within 1 % of `E_u/P` for P ≤ 1.1 W. The internal resistance
   barely matters (`ε_R0 = −0.007`, §10.1).
2. **Cold plus aged plus bursts: a voltage problem.** Stranded charge rises superadditively. At −20 °C, aging
   multiplies it by 3.6 (6.4 % → 23.3 %). With explicit 4 W bursts (S9) the shutdown happens at 44 % displayed.
3. **Gaming: a thermal problem.** Hot ambient lengthens TTE (5.68 h at 25 °C → 10.29 h at 35 °C) at the cost of
   performance. Cool ambient shortens it (4.40 h at 15 °C).

---

## 6. Time-to-empty: definition and numerical method

### 6.1 Definition and events

```
TTE(s0, u, θ, q_loss) = inf{ t > 0 : min_i g_i(x(t), u(t)) ≤ 0 }
g1 = z                       "charge"      (never first in practice, since z_e > 0)
g2 = V − V_cut               "cutoff"      fuel gauge reports empty; orderly shutdown
g3 = V_pk − V_sys            "brownout"    only while ΔP_b(u(t)) > 0 and the burst envelope is enabled
g4 = D                       "infeasible"  demand above P_max; reachable only through a demand step
g5 = T_shut − T              "thermal"     battery at 68 °C
```

* **Switching instants.** `u` is piecewise constant. At each switch, evaluate all `g_i` with the new `u` before
  integrating. If any is ≤ 0, `TTE` is the switch time and the event is the first violated one (for example a
  camera start in the cold).
* **Other outputs.** Report the time to displayed 0 % (`t(s = 0)`), the event type, `s_off`, `ς_off`, the end
  temperatures and the delivered SoC performance.
* **Conditional TTE.** "Remaining time" from any state `x_now` uses the same definition with `x_now` as the initial
  condition.
* **Horizon.** If no event occurs before 2000 h, return `TTE = ∞` with event `none`.

### 6.2 Reference solver (single trajectories)

* `scipy.integrate.solve_ivp(method="LSODA", rtol=1e-7, atol=[1e-9, 1e-7, 1e-5], max_step=60 s)`.
* Integrate segment by segment between input breakpoints, restarting the solver at each switch and carrying the
  state.
* Events `g1`–`g5` are terminal with `direction = −1`.
* This is exactly the [PROTO] configuration. A 30 h run takes ≈ 0.1–0.5 s.

### 6.3 Batch solver (Monte Carlo, Sobol, UQ)

* Use the existing `src/phonebatt/integrate.py:rk4_batch`, vectorised over samples (state shape `(3, n)`).
* Step `h = 2 s` by default. `h = 10 s` is allowed for CTMC Monte Carlo and standby after the convergence test T1
  passes. All scenario segment durations are multiples of 10 s.
* Combined event function (one call per step):
  `g = min(z/0.01, (V − V_cut)/0.1, (V_pk − V_sys)/0.1 [if active], (T_shut − T)/10, D/Ũ²)`. It is NaN when
  `D < 0`, which `rk4_batch` already treats as an immediate event. The event type is the argmin at the crossing.
* **[REVIEW 2026-11-17] Steps that end at an input switch need two evaluations of `g`.** First with the input the
  step was integrated with: a crossing inside the step (for example the cut-off at the very end of a burst that
  is followed by a lighter segment). Then with the new input (§6.1): an event that only the new input causes is
  placed at the switching instant itself, not interpolated back into the step. Evaluating only the new input, as
  a literal reading of "one call per step" suggests (and as `final_proto.py:mc_L2` does), misses the first kind
  whenever the next input is lighter: on S9 (−10 °C, SOH 0.8) 5 of 48 parameter samples ended one full 600 s
  period late (+1.1 to +2.4 % TTE, against the 0.2 % of T1).
* CTMC inputs: sample each path as jump times, map them to the step grid, and hold the usage state constant within
  an RK4 step (all four stages use the state at the step start). The error is at most `h` per sojourn and has no
  bias.

**Stiffness and stability.**

| Mode | Eigenvalue |
|---|---|
| RC | `−1/τ1`: −1/20 s⁻¹ (40 °C) to −1/135 s⁻¹ (−20 °C); −1/17.6 s⁻¹ at 45 °C |
| Thermal | `−1/(k_env R_th C_th)` ≈ −1/450 s⁻¹ |
| SOC | `+I·U′(z)/(Ũ·3600·Q)` ~ +10⁻⁶ s⁻¹ (weakly unstable: lower U → higher I; bounded by the events) |

The stiffness ratio is 10³–10⁵, which is mild. Classical RK4 is stable on the real axis for `|hλ| < 2.785`, so
`h ≤ 49 s`; `h = 2–10 s` is chosen for accuracy at switches.

### 6.4 Numerical acceptance tests (pytest)

| ID | Test | Tolerance |
|---|---|---|
| T1 | Convergence: `\|TTE_RK4(h=2 s) − TTE_LSODA\|/TTE` on S4, S7 (25 °C), S9 (−10 °C, SOH 0.8), S10a (−10 °C, SOH 0.8); and h = 2 s vs 1 s | < 0.2 % ; < 0.1 % |
| T1b | CTMC: mean TTE on 100 paths at h = 10 s vs 5 s | < 0.5 % |
| T2 | Electrical energy closure (v1(0) = 0): `∫I·U dt = ∫P dt + ∫I²R0 dt + ∫v1²/R1 dt + ½·C1·v1(t_end)²` with `C1 = τ1/R1`, isothermal run | relative error < 1e-4 |
| T3 | Analytic agreement: isothermal (`C_th = 1e12`), no bursts, constant P ∈ {0.5, 1.06, 1.49, 2.52, 3.5} W: L2 vs F2 | < 0.05 % ([PROTO] ≤ 0.001 %) |
| T4 | Sanity: no NaN on feasible paths; `z` monotone; `T → T_amb + k_env R_th q` at steady state | exact / 1e-3 K |
| T5 | Forced 20 W step: event `infeasible` (NaN path in `rk4_batch`) | exact |

---

## 7. Analytic and approximate results (insight layer)

**F1. Energy bucket (L0).**

```
E_u(z0) = Q_max·∫_{z_e}^{z0} U(z,T_ref) dz   [Wh];      TTE0 = E_u/P̄      (piecewise u: Σ_k P_k Δt_k = E_u)
∂ln TTE0/∂ln P_i = −s_i  (s_i = energy share of component i);     TTE_without_i / TTE = 1/(1 − s_i)
```

`E_u(1) = 17.79 Wh` [PROTO]. The pie chart of energy shares **is** the elasticity chart (from PS).

**F2. Exact quasi-static quadrature (L1).** For constant P and T, with the RC branch equilibrated (`R_Σ = R0 + R1`),
`dt = 3600·Q·dz/I` and `1/I = V/P` give

```
TTE = (Q_max/P)·∫_{z_c(P,T)}^{z0} V_qs(z) dz  [h],     V_qs = [U + √(U² − 4R_Σ P)]/2
z_c solves  U(z_c,T) = V_cut + P·R_Σ(z_c,T)/V_cut     (scalar root, brentq on [1e-6, z0])
```

Inverting the monotone quadrature gives SOC(t) in closed form. [PROTO]: F2 equals the isothermal L2 TTE to
≤ 0.001 % for P = 0.5–3.5 W. F2 at the thermal steady state `T = T_amb + k_env R_th P`, with the throttled
`P*` of F5 where relevant, is the Sobol surrogate (§10.3).

**F3. First-order resistive correction.** Expanding F2 in P:

```
E_u(P) ≈ E_u0·(1 − a·P),     a = [ ∫_{z_e}^{1} R_Σ/U dz + R_Σ(z_e)/U′(z_e) ] / ∫_{z_e}^{1} U dz
∂ln TTE/∂ln P = −1 − aP/(1 − aP)
```

[PROTO] at 25 °C: ohmic part 0.00495 W⁻¹ + cut-off part 0.00322 W⁻¹ = **a = 0.00818 W⁻¹**, with
`U′(z_e) = 8.82 V` per unit SOC. Battery losses therefore cost 0.4 % at 0.5 W, 0.9 % at 1.06 W and 2.9 % at
3.5 W (L1 vs L0 [PROTO]: −0.4 %, −0.9 %, −3.1 %). `a` scales roughly with `a_R(T)` and `(1 + r_growth)`, so it
grows ×3–5 in the cold and on aged packs.

**F4. Brown-out threshold, stranded charge and safe burst budget.** Setting `V_pk = V_sys` in E10:

```
U(z*,T) = V_sys + v1 + R0(z*,T)·(P + ΔP_b)/V_sys,   v1 ≈ R1·P/U,   s_off ≈ (z* − z_e)/(1 − z_e)
ΔP_b,safe(z, T, SOH) = V_sys·(Ũ − V_sys)/R0 − P
```

[PROTO] check: fresh pack, social load at −20 °C (battery at −16.6 °C from the ODE): `z* = 0.0881`,
`s_off = 6.35 %`, against 6.38 % from the full ODE.

Cut-off-limited power `P_cut(V_cut, R_Σ)` for loads ≥ 3τ1, and burst-limited power `P_cut,sys(V_sys, R0)` for
sub-second peaks, in W [PROTO]:

| Condition | P_cut at z = 0.5 / 0.2 / 0.1 | P_cut,sys at z = 0.5 / 0.2 / 0.1 |
|---|---|---|
| 25 °C, new | 20.2 / 12.5 / 9.0 | 28.3 / 18.6 / 13.9 |
| −10 °C, new | 6.9 / 4.5 / 3.4 | 11.0 / 7.4 / 5.6 |
| −10 °C, SOH 0.8 | 4.6 / 3.0 / 2.2 | 7.3 / 4.9 / 3.7 |
| −20 °C, SOH 0.8 | 3.2 / 2.1 / 1.6 | 5.3 / 3.6 / 2.7 |

So a 4 W, 30 s burst cannot be served below z ≈ 0.5 by a −20 °C, SOH-0.8 phone. That is the 44 % shutdown of S9.

**F5. Thermal cap (gaming).** At thermal steady state

```
T_s* = T_amb + k_env·R_skin·[P(T_s*) + q_batt],        TTE_cap ≈ E_u(P*)/P*,   P* = P(T_s*)
sharp-cap bound: P* ≈ max(P_nonSoC + P_fl, (T_lim − ΔT_thr − T_amb)/(k_env R_skin))
dTTE/dT_amb = E_u/(k_env·R_skin·P*²) > 0      whenever the cap binds
step response:  T(t) = T_amb + R_th q + (T0 − T_amb − R_th q)·e^(−t/τ_th)   (HE19 Eq. 9)
```

[PROTO] (with `q_batt` neglected): `P* = 3.116 W` at 25 °C (TTE_cap 5.71 h vs ODE 5.68 h), delivered performance
`ρ_perf = 0.864`; `P* = 1.728 W` at 35 °C (TTE_cap 10.30 h vs ODE 10.29 h), `ρ_perf = 0.531`. The sharp-cap
bound gives 3.27 W and 1.45 W. The game state's non-SoC power is 0.743 W.

**F6. CTMC mean (PDMP).** Let `X(t)` be the CTMC with generator Λ, stationary distribution π and state powers `P_x`.
In L0 the energy is `dE/dt = −P_{X(t)}`.

```
fluid:   μ_fluid = E_u/P̄,   P̄ = πᵀP
exact (L0, Dynkin; Appendix A.6):   E[TTE | X0 = x0] = ( E_u − h(x0) + E[h(X_TTE)] )/P̄
         h = Dm·(P − P̄1),   Dm = (1πᵀ − Λ)⁻¹ − 1πᵀ  (deviation matrix),   E[h(X_TTE)] ≈ Σ_x π_x P_x h_x / P̄
```

The last line is the inspection-paradox correction: depletion happens preferentially in high-power states.

**F7. CLT variance and inverse-Gaussian law.**

```
σ² = 2·Σ_x π_x (P_x − P̄)·h_x   [W²h],     Var(TTE) ≈ σ²·E_u/P̄³
CV_TTE = √(2τ_c/μ)·CV_P,   CV_P² = Var_π(P)/P̄²,   τ_c = σ²/(2·Var_π P)
TTE ≈ IG(mean μ_Dynkin, shape λ_IG = E_u²/σ²)
two-state on/off check:   σ² = 2αβ(P_on − P_off)²/(α + β)³
```

**F8. Hub-and-spoke closed form and variance attribution** (activity chain alone; Appendix A.7).

```
E[L] = 1/λ_on + Σ_j p_j d_j,   P̄ = (P_idle/λ_on + Σ_j p_j d_j P_j)/E[L],   σ² = 2·Σ_j p_j d_j² (P_j − P̄)² / E[L]
share_j = 2 p_j d_j² (P_j − P̄)² / (σ² E[L])
∂ln μ_fluid/∂ln d_j = −π_j (P_j − P̄)/P̄;    ∂ln μ_fluid/∂ln λ_on = −π_idle(1 − π_idle)(P̄_session − P_idle)/P̄
```

**F9. Semi-Markov sessions** (squared coefficient of variation `CV_j²`, `CV_0²` for idle):
`σ² = [Σ_j p_j(1 + CV_j²) d_j² (P_j − P̄)² + (P_idle − P̄)²(CV_0² − 1)/λ_on²]/E[L]`. Lognormal sessions with
`CV_j = 1.5` multiply the session terms by 1.625.

**F10. Burstiness (Jensen) penalty in L1.** `E[TTE]_L1 ≈ E[TTE]_L0·(1 − a·P̄·(1 + CV_P²))`. Bursts near empty also
raise the effective cut-off.

**[PROTO] results for F6–F10.**

| Quantity | Typical (S13) | Heavy (S14) |
|---|---|---|
| P̄ (W); screen-on time | 0.4400; 6.5 h/day | 0.8224; 11.5 h/day |
| μ_fluid / μ_Dynkin / Monte Carlo L0 (20 000 paths) (h) | 40.43 / 41.25 / 41.34 | 21.63 / 22.50 / 22.52 |
| σ² (W²h); CLT SD / Monte Carlo SD (h); CV | 0.2492; 7.21 / 7.04; 0.175 | 0.8388; 5.18 / 5.11; 0.230 |
| τ_c (h); CV_P | 0.234; 1.659 | 0.362; 1.310 |
| 5 / 50 / 95 % quantiles: IG vs Monte Carlo L0 (h) | 30.25 / 40.59 / 54.49 vs 29.58 / 41.46 / 52.69 | 14.72 / 21.85 / 32.50 vs 14.36 / 22.36 / 31.16 |
| Pr(TTE < 24 h): IG vs Monte Carlo L0 | 0.002 vs 0.009 | 0.650 vs 0.623 |
| Monte Carlo L2 (400 paths, h = 10 s): mean, SD, 5/50/95 % (h) | 39.53, 7.05, 28.36 / 39.38 / 51.34 | 21.22, 4.77, 14.13 / 20.87 / 29.94 |
| L2/L0 mean ratio; mean `s_off`; brown-out share | 0.956; 0.70 %; 7 % of paths | 0.942; 0.91 %; 11 % |
| Energy share by activity | social 0.279, game 0.282, nav 0.106, video 0.100, audio 0.078, vcall 0.074, idle 0.058, camera 0.024 | game 0.545, social 0.157, video 0.119, nav 0.071, vcall 0.059 |
| Energy share by component (= −elasticity in L0) | SoC 0.374, display 0.221, radio 0.182, awake 0.084, suspend 0.045, audio/BT 0.043, camera 0.027, GNSS 0.019 | SoC 0.502, display 0.202, radio 0.158 |
| Variance share (hub, Wi-Fi only; σ² 0.236495 by both formulas) | **game 0.807**, nav 0.102, social 0.032, vcall 0.032, audio 0.021 | game 0.94 |
| Usage elasticities of μ_fluid (Wi-Fi only) | `d_game` −0.257, `d_social` −0.169, `d_audio` +0.065, `λ_on` −0.531 | – |

* The L2/L0 ratio (0.956) is lower than the Jensen estimate (`1 − aP̄(1 + CV_P²) = 0.987`). The difference comes
  from game sessions that start cool and run unthrottled for the first minutes, and from 7 % of paths ending in a
  brown-out.
* A longer music session *extends* TTE (`d_audio` elasticity > 0), because it displaces higher-power use.
* Rate scaling (all sojourns × κ, π fixed): Dynkin mean 40.63 / 41.25 / 43.71 h and SD 3.61 / 7.21 / 14.43 h at
  κ = 0.25 / 1 / 4. The SD scales exactly as √κ.

**F11. Aging closed forms.** Under constant conditions, `q_cal = k_c·a_deg(T;E_cal)·g(z)·√t_d`,
`q_H = k_H·a_deg(T;E_H)·√n` and `q_L = k_L·a_R(T;E_L)·e^(β_L(c − c_ref))·√n_ch`. Calendar fade accelerates ×2.14
per +10 K near 25 °C (thermal_aging G5). At 25 °C, `TTE(q) ≈ TTE_new·(1 − q)·(1 − a·γ_R·q·P)`. The second factor is
< 1 % at q = 0.2 and P = 1.5 W, but large in the cold (§5.4 map).

---

## 8. Validation plan (pre-registered)

**Rule.** The datasets below test predictions. They never set parameters, except the one step marked [CAL].
Metrics: `ε_TTE = (TTE_model − TTE_obs)/TTE_obs`; `RMSE_V`; `|Δf_Q|`; coverage of the 90 % UQ interval. All [PROTO]
predictions below were computed with the §4 defaults before any comparison beyond V1.1.

### V1. Phone-level TTE: Bonn iPhones (BONN26 / SCH26, CC BY 4.0), primary

* **Configuration** (§4.7): airplane mode (`w = 0`, `c = none`); OLED on and black (`s_scr = 1`, `APL = 0`,
  `f_r = 60`), so only `P_d0` counts; `x_c = 0`; 25 °C; no bursts. Camera: `c_r = 1`. GPS: `g = 1`.
* **Observed** values are median runtimes to the last logged 1 % (validation_data a.2).
* **V1.1 [CAL]**: 14 Pro idle sets `ΔP_awake` (exact 0.0877 W; default 0.090 W).
* **V1.2 hold-outs**: everything else. **V1.3**: camera + GPS additivity hold-out.

| Device, scenario | Model P (W) | Model TTE (h) | Observed (h) | ε (default, ΔP_awake = 0.090) | ε (blind prior, ΔP_awake = 0.10) |
|---|---|---|---|---|---|
| 14 Pro idle [CAL] | 0.230 | 54.48 | 55.03 | −1.0 % | −5.1 % |
| 14 Pro rear camera | 1.380 | 8.96 | 8.43 | +6.3 % | +5.6 % |
| 14 Pro GPS "best" | 0.500 | 24.98 | 48.49 | **−48.5 %** | −49.5 % |
| 14 Pro camera + GPS | 1.650 | 7.48 | 8.38 | −10.8 % | −11.3 % |
| 13 Pro idle | 0.230 | 51.85 | 50.82 | +2.0 % | −2.2 % |
| 13 Pro camera | 1.380 | 8.52 | 9.53 | −10.5 % | −11.2 % |
| 13 Pro GPS "best" | 0.500 | 23.77 | 29.76 | −20.1 % | −21.7 % |
| 13 Pro camera + GPS | 1.650 | 7.11 | 8.33 | −14.7 % | −15.2 % |

**Pass criteria.** (a) |ε| ≤ 15 % on hold-outs; (b) the observed device min–max inside the 90 % parameter-UQ
interval; (c) model additivity `1/TTE_cxl = 1/TTE_cam + 1/TTE_GPS − 1/TTE_idle` within 1 % ([PROTO] +0.41 %; data
+0.15 %).

**Pre-registered expectation.** Idle and camera pass on both phones. GPS **fails** on the 14 Pro: the data imply only
+0.03 W for GPS "best", against the 0.27 W KW19 prior. The phones were in a closed box in airplane mode, so the
receiver was probably duty-cycled or had no fix [unverified]. We report this as a finding ("GPS matters surprisingly
little on modern iPhones when stationary") and do **not** re-tune. Structural variant: duty `g ∈ [0.1, 1]`.

**V1.4 SOC-shape test (display-map test).** At constant P, the coulombic map D1 predicts drain slopes proportional
to the mean OCV in each display window: 4.059 / 3.756 / 3.645 V over 99–50 / 50–20 / 20–1 % [PROTO], i.e. slopes in
the ratio 1 : 0.925 : 0.898 (≈ 10 % faster near empty). D2 (energy fraction) predicts a constant slope. Bonn
observes 31.1 / 30.8 / 32.4 min per % (≤ 6 % variation, slightly *slower* near empty). Compute the slope-error
metric for D1 and D2 (both with the top hold `h = 0.065`) on every Bonn idle, GPS and camera trace in
`docs/research/data/validation/bonn_iphone_summary_soc_curves.csv`. Expected: the data favour D2. That supports
reporting the glossary's energy-fraction SOC.

### V2. Cell temperature behaviour: LG HG2 (KOL20, CC BY 4.0), primary

* **CAL at 25 °C only** (§4.7). Then predict 0.5C and 2C discharges at −20, −10, 0, 10 and 40 °C with RP-63's
  activation energies. Report `f_Q`, `f_E`, mean voltage and `RMSE_V`.
* **Targets** (validation_data a.1): `f_Q` at 0.5C = 0.73 / 0.86 / 0.91 / 0.93 / 0.97, and at 2C =
  0.49 / 0.73 / 0.88 / 0.93 / 0.98 (−20 / −10 / 0 / 10 / 40 °C); at C/20, −20 °C, `f_Q` = 0.90.
* **Pass**: `|Δf_Q| ≤ 0.05` for T ≥ 0 °C; ±0.10 at −20 °C; `RMSE_V ≤ 50 mV` for T ≥ 0 °C.
* **Pre-registered expectation.** With a 6:1 split, our activation energies give `R_Σ(−20 °C) ≈ 0.165 Ω` against the
  apparent 0.25 Ω (effective E_a ≈ 20.5 kJ/mol vs HG2's apparent 26.4 kJ/mol) [DERIVED]. So L2 alone will
  **under-predict** the cold loss at −20 °C. The KiBaM variant X9 (`κ_K` Arrhenius, `E_a,K = 35 kJ/mol`) tests
  whether diffusion-limited availability closes the gap. Report both.
* **Dynamic check**: one UDDS file at 25 °C and one at 0 °C (RMSE_V).

### V3. Aging → TTE: NASA PCoE B0005 (NASA07; licence "not specified", documented), secondary

Insert the measured fade (1.862 → 1.331 Ah, −28.7 %) and the measured `R_e`/`R_ct` growth into `Q_max` and `R`.
Predict the aged TTE at 2 A, observed 55.8 → 39.7 min (−28.9 %). Capacity fade alone gives −28.5 %; the remainder
must come from resistance growth through the cut-off. This tests the structural claim that aging enters as Q↓ and
R↑, independent of the rate law. Temperature sanity checks: B0055 (4 °C, 2 A, ratio 0.71 to 24 °C) and B0029
(43 °C, 4 A).

### V4. RC dynamics under pulsed load: Oxford ExampleDC (OX17, ODbL), optional

Drive the cell ECM, scaled to 0.74 Ah at 40 °C, with the 1 Hz Artemis current (peaks 5.0 A = 6.8C). Take the OCV
from low-current segments (CAL). Report `RMSE_V`.

### V5. Physics reference: PyBaMM DFN (BSD-3; Ai2020 and Marquis2019 sets), structural

Constant-power discharges at 0.1–1C-equivalent, at 25, 0 and −10 °C, plus 10 s pulses. Compare TTE and V(t) with L1,
L2 (1RC), X5 (2RC) and X9 (KiBaM), all using the same OCV. Acceptance: TTE difference ≤ 3 % at ≥ 0 °C. In the cold,
report the difference as structural error.

### V6. Published-specification plausibility band (cited facts, not open data)

Tolerance ±25 % (undisclosed brightness, network and vendor differences). Workloads: web browsing = SO or SO-L at the
150-nit brightness; local video = VW without network traffic (`λ = 0`, `d_dl = 0`); Apple tests add BT headphones
(`bt = 1`, `au = 0`); "idle at minimum brightness, WLAN off" = screen on, `b = 0`, LTE camping, no traffic.

| Published case (source) | Observed | Model [PROTO] | ε | Pass |
|---|---|---|---|---|
| iPhone 14 Pro local video (APL14P) | 23 h | 26.08 h (0.479 W) | +13.4 % | ✓ |
| iPhone 14 Pro streamed video (APL14P) | 20 h | 23.40 h (0.534 W) | +17.0 % | ✓ |
| iPhone 14 Pro audio, BT (APL14P) | 75 h | 55.90 h (0.224 W) | −25.5 % | ✗ (marginal) |
| iPhone 14 Pro Max Wi-Fi web, 150 nits (NBC22) | 15.8 h | 15.18 h (1.068 W) | −3.9 % | ✓ |
| iPhone 14 Pro Max local video, 150 nits (NBC22) | 28.13 h | 33.01 h | +17.3 % | ✓ |
| iPhone 14 Pro Max idle, minimum brightness (NBC22) | 57.9 h | 66.63 h (0.245 W) | +15.1 % | ✓ |
| iPhone 14 Pro cellular web, 150 nits (LMAG22) | 10.37 h | 10.29 h (1.203 W) | −0.8 % | ✓ |
| Pixel 8 Wi-Fi web (NBC23) | 13.57 h | 16.37 h | +20.6 % | ✓ |
| Pixel 8 local video (NBC23) | 22.72 h | 35.70 h | +57.1 % | ✗ |
| Pixel 8 idle, minimum brightness (NBC23) | 24.88 h | 70.52 h | +183 % | ✗ |
| Galaxy S24 Wi-Fi web (NBC24) | 14.67 h | 14.30 h | −2.5 % | ✓ |
| Galaxy S24 local video (NBC24) | 20.45 h | 31.20 h | +52.6 % | ✗ |
| Galaxy S24 idle, minimum brightness (NBC24) | 39.97 h | 61.64 h | +54.2 % | ✗ |
| Pixel 8a typical mixed use (GPX, a lower bound) | ≥ 31 h | 34.67 h (S12) | ≥ bound | ✓ |
| Max load or gaming (NBC22–24) | 2.4–5.8 W | capped GM 3.12 W | inside | ✓ |
| Standby drain band (device_power BG9) | 0.15–0.5 %/h | 0.23 %/h (S1) | inside | ✓ |

**Reading.** 11 of 16 pass. Every iPhone case passes except audio, which is marginal. The model's full-wakelock
audio state costs 0.224 W, against 0.165 W implied by Apple's figure. Variant X14 (audio offloaded to a DSP,
`a = 0`) gives 0.138 W and +21.6 %, so the truth lies between the two. The Android idle and video cases fail: their
screen-on floors are 1.6–2.9× ours (NBC idle 0.39–0.71 W vs 0.29 W on the iPhone 14 Pro Max). This vendor spread
of the floor is the largest cross-device uncertainty. The prior range of `ΔP_awake` covers part of it, and we do not
tune to it.

### V7. Published phone-level cold phenomenology (LEE20, HE19, AP-TEMP), check

Configure a 2015-class phone: `Q = 2.7 Ah` [unverified], `R0 + R1 ≈ 0.25 Ω` at 20 °C (LEE20 C3), UI load with 1–2 A
peaks (`ΔP_b = 3 W`). Requirements: brown-out at 15–35 % displayed SOC at ≤ 0 °C (LEE20 C9/C10); a 0 °C TTE loss
of the same order as LEE20's −15 % (159 → 135 min); the stranded charge comes back after warming (C12: +330 mAh),
which happens automatically because `z` is not reset.

---

## 9. Uncertainty quantification

**Sources of uncertainty and how each is treated.**

| # | Source | Treatment |
|---|---|---|
| U1 | Parameters θ (epistemic) | the §4 distributions, sampled jointly |
| U2 | Usage path (aleatory) | CTMC paths at fixed θ (S13, S14) |
| U3 | Usage-model parameters | `λ_on`, `p`, `d`, connectivity means (§4.6) |
| U4 | Environment | `T_amb` per scenario (e.g. "winter outdoor" `U(−15, 5) °C`); `k_env` |
| U5 | Initial-SOC reading | ±2 % displayed [ASSUMPTION], through E11 |
| U6 | Model structure | the X-variants of §10.4, reported as a separate spread and never mixed into θ |
| U7 | Observation (Bonn only) | display map D1/D2, 1 % quantisation, top hold `h` |

**Propagation.**

1. **Deterministic scenarios** (S1–S12, S16; and S15 at s0 = 50 / 20 %). Scrambled Sobol' QMC
   (`scipy.stats.qmc.Sobol`, `scramble=True`, `seed=2026`), `N = 2048`, mapped through the inverse CDFs, run with
   `rk4_batch` (h = 2 s; S1/S2 h = 10 s). Report median, 5–95 % interval, CV, `Pr(TTE < t*)` and `s_off`
   quantiles.
2. **First-order check (L0 delta method).** `CV_TTE² ≈ CV_Eu² + Σ_i s_i²·CV_{P_i}²`, compared with item 1.
3. **Stochastic scenarios** (S13, S14). Law of total variance:
   `Var(TTE) = E_θ[Var(TTE|θ)] + Var_θ(E[TTE|θ])` (first term = usage, second = parameters). The inner moments
   come from F6–F7 per θ sample (state powers recomputed from θ; milliseconds each), with `N_θ = 2048`. Check on
   50 θ × 200 paths of full L2: mean within 3 % and SD within 10 % after applying the L2/L0 mean ratio.
4. **Bonn-constrained posterior (optional ABC).** Keep θ samples whose 14 Pro idle TTE is within ±5 % of 55.0 h.
   Report how much the RP-63 intervals narrow. This makes the single calibration explicit and quantitative.
5. **Target event.** `t* = 16 h` for S12 and S13 started at 07:00, i.e. "dies before 23:00".

**Deliverables.** Table T3: TTE median [5 %, 95 %] for the constant scenarios × `s0` ∈ {100, 50, 20} % × `T_amb`
∈ {25, 0, −10} °C. Violin figure. A statement of the form "for a typical day, usage randomness contributes X % of
the TTE variance and hardware parameters Y %".

---

## 10. Sensitivity analysis

### 10.1 Local elasticities

`ε_θ = ∂ln TTE/∂ln θ` by central differences at ±5 % (voltages: % TTE per +10 mV). Use constant-load scenarios and
the smooth CTMC mean only. Under periodic schedules TTE is **non-smooth** in θ: the event snaps to the next heavy
segment, which produces spurious elasticities. For S9 and S12 use the OAT ranges of §10.2 or the smooth metrics of
§12.1.

[PROTO] elasticities (RP-63; ΔP_b multiplier = `m_b`):

| Parameter | S4 social, 25 °C, fresh | S4 social, −10 °C, SOH 0.8 | S7 gaming, 25 °C |
|---|---|---|---|
| `Q_n` | +1.000 | +1.000 | +1.009 |
| `P_SoC,max` | −0.406 | −0.431 | −0.180 (capped) |
| `p_c` (display content) | −0.119 | −0.127 | −0.012 |
| `P_d0` | −0.114 | −0.121 | −0.010 |
| `ΔP_awake` | −0.086 | −0.091 | −0.007 |
| `R0,ref` | −0.007 | **−0.274** | −0.012 |
| `E_a,R0` | +0.001 | **−0.236** | +0.003 |
| `R1,ref` | −0.001 | −0.018 | −0.001 |
| `m_b` | 0.000 (not binding) | **−0.190** | 0.000 |
| `R_skin` | 0.000 | 0.000 | **+0.767** |
| `R_th` | 0.000 | +0.030 (self-heating helps) | −0.007 |
| `C_th` | 0.000 | 0.000 | −0.009 |

**Reading.** At room temperature TTE is an energy-budget problem: Q and the load parameters matter, the cell
resistance does not. In the cold on an aged pack it becomes a voltage problem: `R0`, `E_a,R0` and the burst size
dominate. Under sustained gaming it is a thermal-design problem: `R_skin`.

### 10.2 One-at-a-time tornado

Move each θ to its range endpoints with the others at default. Scenarios: S4, S7, S10a (−10 °C, SOH 0.8), S12 (on
the mean-power metric) and S13 (Dynkin mean). Put the usage factors on the same chart:

* brightness `b ± 0.2`; dark theme `m_APL = 0.25`; refresh 60 vs 120 Hz;
* background rate `λ_bg × {0.5, 5}`; RSRP −85 → −115 dBm; RAT Wi-Fi / LTE / 5G;
* `T_amb` −10 → 35 °C; `k_env` 0.8 → 1.5.

This puts hardware uncertainty and user behaviour on one axis.

### 10.3 Global variance-based sensitivity (Sobol)

* SALib Saltelli sampling with Jansen estimators. 95 % bootstrap CIs from 1000 resamples.
* **24 factor groups** (each group moved jointly by its own §4 distributions):
  1 `Q_n`; 2 `R0,ref`; 3 `A_R`; 4 `E_a,R0`; 5 RC {`R1,ref`, `τ1,ref`, `E_a,R1`, `E_a,τ`}; 6 OCV {`δU`, `k_S`};
  7 `V_cut`; 8 `Δ_sys`; 9 `C_th`; 10 `R_th`; 11 `R_skin`; 12 throttle {`T_lim`, `w_thr`, `P_fl`}; 13 `k_env`;
  14 floor {`P_susp`, `ΔP_awake`, `E_wake`}; 15 display {`P_d0`, `p_c`, `γ_b`}; 16 `κ_RR`;
  17 SoC {`P_SoC,max`, `v0`}; 18 Wi-Fi radio; 19 LTE radio; 20 NR radio; 21 `γ_sig`;
  22 peripherals {`P_G`, `ΔP_cam`, `ΔP_camf`, `ΔP_aud`, `ΔP_BT`}; 23 `m_b`; 24 `γ_R`.
* **Surrogate for constant-load scenarios** (S3–S6, S8, S11): F2 at the thermal steady state, plus F4 for the
  brown-out root and F5 for the throttled power. `N = 1024`, i.e. `N(2k + 2) = 51,200` evaluations per scenario.
  Check against full L2 on 200 random points; the error must be ≤ 1 %.
* **Full L2** for S7, S9 and S10a: `N = 256` (12,800 runs each, vectorised in `rk4_batch`).
* **Interaction study** (S10a-G): add `T_amb ~ U(−20, 10) °C`, `SOH ~ U(0.75, 1.0)` and `m_b` as explicit factors
  (k = 26 with group 23 removed). Report first-order `S_i`, total `S_Ti`, second-order `S_ij` within
  {T_amb, SOH, m_b}, and the interaction share `1 − Σ S_i`.

### 10.4 Structural sensitivity (model-form uncertainty)

Run every scenario under each variant. Report ΔTTE (%) and Δ`s_off` relative to X0 (Table T4).

| ID | Variant | Question tested |
|---|---|---|
| X0 | baseline L2/L3 (this spec) | – |
| X1 | L0 energy bucket | is a voltage model needed? |
| X2 | L1 quasi-static ECM, isothermal at `T_amb` | are RC and thermal dynamics needed? |
| X3 | L2 isothermal (no thermal ODE, no throttling) | self-heating in the cold; throttling in gaming |
| X4 | constant-current load `I = P/3.87 V` | the constant-power (PMIC) assumption A1 |
| X5 | 2RC (2.5 mΩ / 6.7 s + 14.1 mΩ / 47 s at 25 °C, scaled by `R·Q`) | RC order |
| X6 | `R0` split: 35 mΩ temperature-independent PCM + 25 mΩ Arrhenius cell | where the resistance lives (A11) |
| X7 | no burst envelope (cut-off only) | the brown-out mechanism |
| X8 | single threshold `V_sys = V_cut` (PS-style) | the two-threshold assumption |
| X9 | KiBaM instead of RC: `c = 0.166` (JH09), `κ_K = 40 h⁻¹` (33–47, upper bound [DERIVED] from DFN rate loss), `E_a,K = 35 kJ/mol` (25–55) | rate-capacity and diffusion in the cold |
| X10 | RRC tail as an ODE state `dr/dt = λ(1 − r) − r/T_tail` | quasi-steady tail formula A7 |
| X11 | brightness law `γ_b = 1` (CH13) and `γ_b = 2.9` (DH21-derived) | display law |
| X12 | reported SOC D2 instead of D1 for `s0` and `s_off` | output definition |
| X13 | SoC leakage `P_leak = 0.1 W·exp(k_lk(T_s − 25 °C))`, `k_lk ∈ [0, 0.04] K⁻¹` [unverified] | heat-driven power growth |
| X14 | audio offloaded to a DSP (`a = 0` in AU) | V6 audio mismatch |
| X15 | lognormal sessions with CV 1.5 (semi-Markov) | usage law A13 |
| X16 | diurnal CTMC: `λ_on = 0` from 00:00 to 07:00, 5.5 h⁻¹ otherwise | time-homogeneity |
| X17 | calendar law `t^0.75` (SCH14) instead of √t | aging-law shape |
| X18 | OCV stretch pairs (0.192 V, 5.44) and (0.272 V, 2.90) | constructed 4.45 V OCV (A10) |

### 10.5 Usage-fluctuation sensitivity

1. **Rate scaling.** Scale all sojourns by κ ∈ {0.25, 0.5, 1, 2, 4} with π fixed. F7 predicts SD ∝ √κ; verify by
   Monte Carlo ([PROTO] analytic: SD 3.61 / 7.21 / 14.43 h at κ = 0.25 / 1 / 4).
2. **Burstiness at fixed mean.** A constant 1.2 W load vs a 50/50 alternation of 0.4 W and 2.0 W with periods from
   10 s to 1 h, at 25 °C fresh and at −10 °C on an SOH-0.8 pack. L0 TTE is unchanged; L2 differs through `I²R`, the
   cut-off and the brown-out.
3. **Order effect.** S12 activities in reversed order (heavy use at low SOC), at 25 °C and at −10 °C aged.
4. **Initial SOC.** TTE(`s0`) is nearly linear for `s0 ≥ 0.2` at 25 °C. Report the curvature in the cold, where the
   effective empty point `z*` rises (F4).
5. **Heavy-tailed sessions** (F9, X15) and **per-visit jitter** (±30 % on `b`, APL, `λ`, `d_dl`).
6. **Usage-model parameters.** `λ_on`, `p`, `d` with the §4.6 ranges, through F6–F8.

---

## 11. Aging analysis (L3)

### 11.1 Calibration (protocol C0)

C0: cell at a constant 25 °C. Each day: discharge 100 → 0 % linearly over 16 h; charge per E21 (CC 1C to 80 %,
0.3C to 100 %); hold at 100 % until the next day. That is 1 FCE/day.

1. `k_c = 0.027/√365 = 1.413e-3 d^−½` (2.7 %/yr at g = 1).
2. Under C0, `q_cal(1000 d) = 0.0782` (mean `g` = 1.749) [PROTO].
3. Set `q_H(1000) + q_L(1000) = 0.20 − 0.0782`, with `q_L` share `w_L = 0.3`. This gives
   `k_H = 2.697e-3` and `k_L = 1.288e-3` FCE^−½ [PROTO], and SOH(1000 FCE) = 0.8000 by construction.
4. Variants: `N_80` = 500 or 800 FCE (older iPhones and Pixels).

### 11.2 Daily driver (quasi-static day profile, production default)

* Segments of ≤ 1 h. Days start at 07:00 unplugged.
* 07:00–23:00: SOC falls linearly by `ΔSOC = min(0.95, 24 h·P̄_day/(E_n(1 − q_loss)))`, with
  `P̄_day = 0.440 W` (S13 typical). Cell temperature `T_use = T_amb + k_env·R_th·P̄_day` (1.32 K above ambient).
* 23:00: plug in. Charge per E21 at `T_ch = T_amb + ΔT_ch`, then hold at the target at `T_amb` until 07:00.
  Optimised charging (OBC): charge to 80 %, hold, resume 1 h before 07:00.
* Throughput per segment: `Δn = ΔSOC/2` (discharge) and `c·Δt/2` (charge); `Δn_ch = c·Δt`.
* Update `Q_cal`, `Q_H`, `Q_L` with E17–E19, which is exact for piecewise-constant segments. Update `Q_max` and `R`
  once per day.
* Optional exact variant: re-run one fast L2 day (with charging) whenever SOH has changed by 0.5 %.

### 11.3 Aging scenarios and [PROTO] results (25 °C unless stated)

| ID | Charging and environment | SOH at 1 / 2 / 3 years | Years to SOH 0.8 |
|---|---|---|---|
| A1 | charge to 100 %, held overnight | 0.883 / 0.834 / 0.795 | 2.86 |
| A2 | 80 % charge limit (AP-OBC, SAM) | 0.897 / 0.853 / 0.818 | > 3 |
| A3 | optimised charging (hold 80 %, finish 1 h before unplug) | 0.889 / 0.841 / 0.804 | > 3 |
| A4 | hot climate, `T_amb = 35 °C` | 0.796 / 0.709 / 0.641 | 0.96 |
| A5 | 2 h/day at 45 °C (car dashboard), otherwise A1 | 0.841 / 0.774 / 0.721 | 1.58 |
| A6a | fast charge: CC 1.5C to 80 % | 0.858 / 0.796 / 0.747 | 1.93 |
| A6b | charging in a 0 °C room at 1C (not JEITA-compliant) | 0.849 / 0.783 / 0.730 | 1.72 |
| A6c | charging in a 0 °C room at 0.5C (JEITA rule) | 0.896 / 0.853 / 0.819 | > 3 |
| A7 | 1 FCE/day (heavy user, ΔSOC 0.95) | 0.872 / 0.820 / 0.779 | 2.46 |

At 2 years, A2 gains +1.9 pp of SOH over A1 and A3 gains +0.7 pp, while heat (A4) costs −12.5 pp. The A6 numbers
rest on SCM18's LFP `β_L`, which is only an upper bound (§4.5); report them with the `β_L ~ U(0.5, 2.64)` band.

**Outputs.** `SOH(t)` and `R(t)/R_new` per scenario; TTE of S13 (Dynkin mean) and S12 against age; `s_off` at 0 °C
and −10 °C (S4) against age ("when does a winter shutdown start?"); years to 0.8. Figure `fig_aging.pdf`.

**Coupling statement for the paper.** Aging reduces TTE through two channels. Q↓ is proportional and benign. R↑ is
nonlinear and invisible at 25 °C (`ε_R0 = −0.007`), but it dominates in the cold (§5.4 map, F4). That is why
phones at "80 % battery health" shut down unexpectedly in winter (AP-BATT, LEE20).

---

## 12. Recommendations

### 12.1 User levers (ranking metric)

**Metric.** TTE of a fixed schedule is non-smooth: removing 1 % of the energy can move the shutdown across a
low-power night (a large shift) or within a video call (a small shift). Levers are therefore ranked on two smooth
metrics:

* `E_u/P̄` of the S12 day, with quasi-steady state powers (game at its thermal fixed point);
* the Dynkin mean of the S13 CTMC (F6).

Both are in % TTE gain. The production code also reports the L2 Monte Carlo median with a 90 % interval (§9).
Why not the L2 TTE of S12? Its shutdown lands inside the evening video call of day 2 (P ≈ 2.5 W), so the L2 TTE
moves by only +1.4 % for dark theme, +0.8 % for brightness −0.1 and +0.0 % for halved background traffic
[PROTO], although the daily energy falls by 7.1 %, 4.2 % and 0.1 %. Skipping the game moves it by +37.7 %.

| Lever (applied to all states where it makes sense) | S12 `E_u/P̄` | S13 Dynkin mean |
|---|---|---|
| Skip 3D gaming (GM → standby) | **+34.9 %** | **+36.9 %** |
| Dark theme (APL × 0.25 on RD, SO, SO-L, NV) | +7.6 % | +4.6 % |
| Voice call instead of video call | +4.8 % | +4.0 % |
| Brightness −0.1 on all screen-on states | +4.4 % | +3.1 % |
| Cap refresh at 60 Hz | +4.1 % | +4.7 % |
| Video call over Wi-Fi instead of LTE | +3.9 % | +1.6 % |
| Prefer Wi-Fi: Wi-Fi share 50 → 80 % | – | +3.7 % (41.25 → 42.78 h) |
| Halve background wake-ups | +0.1 % | +0.4 % |
| *Avoid*: standby on 5G NSA without Wi-Fi | −4.4 % | −4.5 % |
| *Avoid*: weak signal (−110 dBm) all day | −3.9 % | −5.7 % |
| *Avoid*: brightness +0.2 (outdoors) | −10.3 % | −7.5 % |

**Scenario-specific levers [PROTO].** Video over Wi-Fi instead of 5G NSA: 32.32 vs 21.10 h (+53 %). A video call over
Wi-Fi instead of LTE: 10.83 vs 6.93 h (+56 %). In the cold, keeping the phone warm prevents shutdowns at 7–44 %
displayed (S9). Gaming in a hot car keeps playing time but cuts performance to `ρ_perf` = 0.53 (F5).

**Energy shares of S12 by component** [PROTO]: display 32.1 %, SoC 29.4 %, radio 17.0 %, awake overhead 8.5 %,
suspend floor 4.8 %, GNSS 4.0 %, audio 2.0 %, camera 1.5 %, sensors 0.7 %.

**What drains most.** Sustained SoC load (gaming: 28 % of typical-day energy and 81 % of its TTE variance); bright,
light-themed, high-refresh screens; real-time cellular uplink (video calls, worse at weak signal: −23 % at
−110 dBm); 5G NSA tails in standby and streaming; and, as a condition, cold + aging + bursts, which strands up to
44 % of the displayed charge.

**What matters surprisingly little.**

* Halving background sync (+0.1–0.4 % per day), because modern 4G tails are cheap (NA21 vs HU12).
* GNSS: 1.9 % of typical-day energy; Bonn shows +0.03 W on a stationary iPhone 14 Pro.
* Bluetooth audio (0.04 W) and sensors (47 mW only while navigating).
* Cell resistance and rate-capacity effects at 25 °C: −0.4 % to −3 % (F3); KiBaM < 0.5 % (X9).
* Cold on a **fresh** pack at moderate load: −3.4 % for navigation at −10 °C (self-heating).
* Self-discharge (0.1 %/month) and optimised charging (+0.7 pp SOH at 2 years).

### 12.2 OS power-management strategies (each quantified in the model)

| # | Strategy | Model mechanism | Quantification [PROTO where given] |
|---|---|---|---|
| O1 | **Peak-power governor** (brown-out avoidance; cf. Apple performance management, AP-BATT) | cap `ΔP_b ≤ max(0, 0.9·ΔP_b,safe)` with F4 evaluated online | aged pack at −10 °C, social: `m_b` 1 → 0.5 cuts `s_off` from 8.5 % to 2.6 % and gains +6.0 % TTE (12.15 → 12.88 h); `m_b = 0` gives 1.5 % and 13.00 h. Report the performance cost (capped burst energy) |
| O2 | Thermal-aware frame cap instead of reactive throttling | game `x_c` 0.80 → 0.64 (−20 % frequency) | steady game power 3.116 → 2.722 W; +3.4 % (S12) / +2.7 % (S13) TTE; constant performance instead of throttle oscillation |
| O3 | Tail-aware coalescing of background traffic (×5) | `λ → λ/5` in E14 for standby states | standby TTE: cellular +25.4 %, 5G NSA +66.7 %, home Wi-Fi +2.4 %; day level +0.1–0.7 % |
| O4 | RAT steering: LTE or Wi-Fi for low-throughput states, 5G only for bulk (XU20: −24.8 %) | `P_tail`, `E_pro`, `P_act` | avoided losses: 5G standby −4.5 %; streaming over 5G −35 % vs Wi-Fi |
| O5 | Signal-aware deferral of bulk transfers while RSRP < −100 dBm | `m_sig` | avoided: −5.7 % (S13 at −110 dBm); −12 % for social over LTE |
| O6 | Adaptive refresh (LTPO): 120 Hz only while scrolling | `f_r` | +4.1–4.7 % |
| O7 | Cold-aware fuel gauge | remap displayed 0 % to `z*(T, SOH, ΔP_b)` (F4); warn "X % unavailable until warm" | removes surprise shutdowns at displayed > 0 (S9, §5.4 map) |
| O8 | Model-based remaining-time display with an uncertainty band | EKF on (`z`, `v1`) from V, I, T; online `P̄`, σ² and the IG law (F6–F7) | "4 h 10 min (3 h 20 min – 5 h 00 min)"; calibrate the interval coverage on S13 paths |
| O9 | Charging defaults | E17–E21 | 80 % limit +1.9 pp SOH at 2 years; OBC +0.7 pp; keep the phone cool while charging (A4/A5: −6 to −12.5 pp); JEITA cold rule (A6c vs A6b: +7.0 pp) |

---

## 13. Generalisation to other portable devices

**Same equations, different constants.**

| Quantity | Scaling rule |
|---|---|
| Energy | `E_u = N_s·Q·Ū` for `N_s` cells in series (laptops); `TTE = E_u/P̄` still holds |
| Pack resistance | `R ∝ (N_s/N_p)/Q_cell` (A12); RP-63 has `R0·Q = 0.28 Ω·Ah` including the PCM |
| Thermal | `C_th ∝ mass`; `R_th, R_skin ∝ 1/(h·A_surface)`, so `P_sust ∝ (T_lim − T_amb)·h·A_surface` |
| Display | OLED `∝ A·g(b)·APL`; LCD `P_panel + c1·b + c2·b²`; e-paper: energy per refresh |
| Usage | a device-specific activity set; the same PDMP formulas F6–F9 |

**Dimensionless groups** (RP-63 values [DERIVED from PROTO]).

| Group | Meaning | RP-63 |
|---|---|---|
| `Π1 = P·R_Σ/Ū²` | resistive loss fraction | 0.005 at 1 W, 25 °C |
| `Π2 = P_dem/P_sust` | throttling index (capped if > 1) | 1.19 for gaming at 25 °C |
| `Π3 = τ_th/TTE` | thermal memory | 0.02 (gaming) |
| `Π4 = ΔP_b/(P_cut,sys − P)` | brown-out margin (failure as it → 1) | 0.11 (S4, 25 °C, z = 0.5); 1.2 (S4, −20 °C, SOH 0.8, z = 0.2) |
| `Π5 = CV_TTE = √(2τ_c/μ)·CV_P` | usage predictability | 0.175 (typical day) |

Device classes are placed on a (Π1, Π4) and (Π2, Π5) map in `fig_generalization.pdf`. The values are orders of
magnitude [unverified], since no device data were collected for other classes:

* **Tablets**: display-dominated; larger area relaxes throttling (smaller Π2).
* **Laptops**: multi-cell, small Π1, large thermal mass, CPU-dominated.
* **Smartwatches and earbuds**: tiny `C_th` (τ_th ≈ a minute); standby and BLE floor dominate, so the floor's
  elasticity → −1; small cells give a large Π4 in the cold.
* **E-readers**: `TTE ≈ E_u/(P_sleep + ν_page·E_page)` with page turns as a Poisson process.
* **Drones**: near-constant hover power; large Π1 and Π4.

---

## 14. Analysis registry: what each analysis produces

All analysis scripts live in `analysis/`, import `phonebatt`, write one JSON file to `results/` and PDF figures to
`paper/figures/`. `analysis/run_all.py` runs them in the order below. JSON numbers are SI except where the key
says otherwise (`_h`, `_pct`, `_C`, `_W`, `_Wh`). Every JSON file has a top-level `"meta"`: git SHA, spec version
`"FINAL_SPEC 2026-11-17"`, seed, solver and h.

| ID | Analysis | Script | Results file and keys | Figures / tables |
|---|---|---|---|---|
| AN0 | Resolved parameter table | `run_params.py` | `params.json`: `{symbol: {value, unit, dist, lo, mode, hi, source, conf}}` | `paper/tables/tab_params.tex` (T1) |
| AN1 | Primitive powers and component breakdown | `run_power.py` | `power.json`: `primitives.{ID}.{P_W, P_capped_W, components{...}}` | `fig_power_breakdown.pdf` (stacked bars per primitive with V6 markers); `tab_scenarios.tex` (T2) |
| AN2 | Cell maps | `run_cell_maps.py` | `cell_maps.json`: `ocv_table`, `R0_grid_mOhm`, `Pcut_W`, `Pcut_sys_W` | `fig_cell_maps.pdf` (U(z,T), R0(z,T), P_cut,sys(z,T) with S4/S6/S7 operating points) |
| AN3 | Deterministic scenarios (§5.2) | `run_scenarios.py` | `scenarios.json`: key `"{S}\|Tamb={T}\|soh={SOH}\|s0={s0}"` → `{tte_h, event, s_off_pct, soe_off_pct, t_disp0_h, P_mean_W, T_end_C, Ts_end_C, rho_perf}` | `fig_trajectories.pdf` (s, V, T, T_s for S1, S5a, S7, S9); `fig_tte_s0.pdf` (S15); `fig_throttle.pdf` (S7 at 15/25/35 °C: T_s(t), P(t), ρ_perf) |
| AN4 | Analytic checks | `run_analytic.py` | `analytic.json`: `F2_vs_L2_pct{P}`, `a_coeff{ohm, cut, total}`, `F4_check`, `F5_cap{Tamb}`, `Pcut_table` | table inside §4 of the paper |
| AN5 | CTMC (S13, S14) | `run_ctmc.py` | `ctmc.json`: `{profile}.{Pbar_W, sigma2_W2h, mu_fluid_h, mu_dynkin_h, sd_clt_h, ig_q05_q50_q95_h, mc_L0{mean_h, sd_h, q05_q50_q95_h, p_lt_24h}, mc_L2{...}, energy_share{}, component_share{}, variance_share{}, elasticities{}}` | `fig_ctmc_tte.pdf` (activity raster, SOC paths, TTE histogram with IG density and Dynkin mean); `fig_cv_scaling.pdf` (CV vs κ, exponential and lognormal) |
| AN6 | Validation V1–V7 | `run_validation.py` | `validation.json`: `V1.{dev}\|{case}.{model_h, obs_h, err_pct, P_W}`, `V1.additivity_residual_pct`, `V1.slope_test{D1, D2}`, `V2.{T}\|{rate}.{fQ_model, fQ_obs, rmse_V}`, `V3`, `V4.rmse_V`, `V5.{case}.dtte_pct`, `V6.{case}.{model_h, obs_h, err_pct, pass}`, `V7` | `fig_validation_bonn.pdf` (predicted vs observed with 90 % UQ bars; slope test inset); `fig_validation_hg2.pdf` (f_Q(T) at 0.5C and 2C, L2 and X9); `fig_validation_band.pdf` (V6) |
| AN7 | UQ | `run_uq.py` | `uq.json`: `"{S}\|s0\|Tamb".{median_h, q05_h, q95_h, cv, p_lt_target, s_off_q05_q50_q95}`, `variance_split.{S13,S14}.{usage_frac, param_frac}`, `delta_method_check`, `abc_bonn` | `fig_uq_violins.pdf`; `tab_tte_uq.tex` (T3) |
| AN8 | Local elasticities and OAT | `run_sensitivity.py` | `sensitivity.json`: `elasticities."{S}\|{T}\|{SOH}".{param}`, `oat."{S}".{factor}.{lo_h, hi_h}` | `fig_tornado.pdf` (S4 and S10a −10 °C aged; hardware and behaviour on one axis) |
| AN9 | Sobol | `run_sobol.py` | `sobol.json`: `{S}.{S1, ST, S1_conf, ST_conf}`, `S10a-G.{S2, interaction_share}`, `surrogate_check_maxerr_pct` | `fig_sobol.pdf` (S_i / S_Ti heatmap, factors × scenarios) |
| AN10 | Structural variants X1–X18 | `run_structural.py` | `structural.json`: `"{X}\|{S}".{dtte_pct, ds_off_pct}` | `tab_structural.tex` (T4); `fig_structural.pdf` (heatmap) |
| AN11 | Usage fluctuations (§10.5) | `run_usage_fluct.py` | `usage_fluct.json`: `rate_scaling`, `burstiness`, `order_effect`, `s0_curvature`, `heavy_tail`, `jitter`, `usage_params` | `fig_usage_fluct.pdf` |
| AN12 | Cold / aged stranded charge and O1 | `run_stranded.py` | `stranded.json`: `map."{Tamb}\|{SOH}".{s_off_pct, tte_h, event}`, `governor."{m_b}".{s_off_pct, tte_h}`, `S9` | `fig_stranded_map.pdf` (heatmap with and without the O1 governor) |
| AN13 | Aging A1–A7 | `run_aging.py` | `aging.json`: `calibration.{k_c, k_H, k_L, q_cal_1000}`, `{A}.{t_yr[], soh[], r_ratio[], years_to_80, tte_S13_h[], s_off_0C_pct[]}` | `fig_aging.pdf` (SOH vs years per strategy; TTE and s_off vs SOH) |
| AN14 | Recommendations | `run_recommendations.py` | `recommendations.json`: `levers.{id}.{s12_pct, ctmc_pct, ctmc_ci90_pct}`, `scenario_levers`, `os.{O1..O9}` | `fig_levers.pdf` (ranked bars with UQ error bars) |
| AN15 | Generalisation | `run_generalization.py` | `generalization.json`: `RP63.{Pi1..Pi5}`, `classes.{name}.{Pi1..Pi5, notes}` | `fig_generalization.pdf` |
| AN16 | Architecture diagram | `make_architecture.py` | – | `fig_architecture.pdf` (u(t) → loads → ECM ↔ thermal ↔ throttle; burst/brown-out; slow aging; feedback arrows T→R, T_s→P_SoC, SOH→Q,R) |

---

## 15. Implementation notes for `src/phonebatt`

### 15.1 Module map

All numerical functions accept NumPy arrays and broadcast over a trailing sample axis.

| Module | Contents |
|---|---|
| `params.py` | `Params` (frozen dataclass or mapping) with every §4 symbol; `RP63 = default()`; `UQ_SPEC: {name: (dist, lo, mode, hi)}`; device constructors `iphone14pro()`, `iphone13pro()`, `iphone14promax()`, `pixel8()`, `galaxy_s24()`, `pixel8a()`, `hg2_cell()`, `nasa_b0005()`; `sample(spec, n, method="sobol", seed=2026)` |
| `cell.py` | `a_R`, `a_deg`, `ocv(z, T, p)`, `dUdT(z)`, `R0(z, T, p)`, `R1(T, p)`, `tau1(T, p)`, `z_empty(p)`, `usable_energy(p, z0)`, `current(Ue, R0, P) -> (I, V, D)`, `peak_voltage(Ue, R0, Ppk)`, `display_soc(z, p, mode="D1"\|"D2")` and its inverse |
| `power.py` | `PRIMITIVES` (§5.1), `power_parts(u, p) -> (P_nonsoc, P_soc_dem, dPb, breakdown)`, `throttle(P_dem, Ts, p)`, `x_eff(P_soc, p)`, `steady_power(u, p, Tamb)` (F5) |
| `usage.py` | `UsageState` (§1.4 fields), `Schedule` (list of `(dt_s, UsageState)`, cyclic), `scenario(id, **opts)` for S1–S16, `CTMC(profile)` with `generator()`, `state_map`, `sample_paths(n, t_max, h, seed)` |
| `model.py` | scalar `Segment.rhs/state/gvals` for LSODA; vectorised `rhs_batch`, `event_batch`; `run(schedule, Tamb, s0, T0, p, solver="lsoda"\|"rk4", h) -> Result(tte_h, event, s_off, soe_off, t_disp0_h, traj)` |
| `integrate.py` | existing `rk4_batch` (unchanged), plus a thin segmented `solve_ivp` wrapper |
| `analytic.py` | `tte_L0`, `tte_L1` (F2), `a_coefficient` (F3), `brownout_root`, `safe_burst` (F4), `thermal_cap` (F5), `ctmc_moments` (F6–F7), `hub_closed_form` (F8), `ig_quantiles`, `semi_markov_sigma2` (F9) |
| `aging.py` | `g_stress(z, p)` (E16), `day_profile(strategy, Tamb, dsoc, **opts)`, `advance_day(state, profile, p)`, `calibrate(p) -> (k_H, k_L)`, `run_years(strategy, Tamb, years, p)` |
| `validation.py` | Bonn observations (read from `docs/research/data/validation/bonn_iphone_summary.json`), V6 table, HG2 and NASA loaders (read the CSV/JSON already in `docs/research/data/validation/`) |

### 15.2 Interfaces

```python
run(schedule, Tamb_C, s0=1.0, T0_C=None, p=RP63, solver="lsoda", h=2.0, t_max_h=2000.0, record=False)
    -> Result(tte_h: float, event: str, s_off: float, soe_off: float, t_disp0_h: float, traj: dict | None)
run_batch(schedules_or_ctmc, Tamb_C, theta: dict[str, ndarray], s0=1.0, h=2.0)
    -> dict(tte_h=ndarray, event=ndarray[str], s_off=ndarray)
ctmc_moments(L: ndarray, P: ndarray, Eu_Wh: float, x0: int = 0)
    -> dict(pi, Pbar, sigma2, h, mu_fluid, mu_dynkin, sd_clt)
```

### 15.3 Regression tests (pytest; values [PROTO]; relative tolerance unless stated)

| Quantity | Expected | Tolerance |
|---|---|---|
| `U(1, T_ref)`; mean OCV; `z_e` | 4.4496 V; 3.8699 V; 0.0263 | 1e-4 abs |
| `E_n`; `E_u` | 18.189 Wh; 17.788 Wh | 1e-3 |
| `R0(0.5, −20 °C)`; `R1(−20 °C)`; `τ1(−20 °C)` | 218.1 mΩ; 86.0 mΩ; 135.3 s | 1e-3 |
| Primitive powers (§5.1) | table | ±0.5 mW |
| Steady game power (F5) at 25 / 35 °C | 3.116 / 1.728 W | ±1 mW |
| L2 TTE and event (§5.4 main table) | table | 1 %; event identical |
| S9 at −20 °C, SOH 0.8 | 12.33 h, `s_off` 44.2 % | 2 %; ±1 pp |
| Stranded-charge map (§5.4) | table | ±0.3 pp; TTE 1 % |
| Bonn predictions (§8 V1) | table | 1 % |
| V6 band (§8) | table | 1 % |
| F2 vs isothermal L2 (T3) | ≤ 0.001 % | ≤ 0.05 % |
| `a` (F3) | 0.00818 W⁻¹ | 2 % |
| F4 check (−20 °C social) | `z*` 0.0881 | 1 % |
| CTMC typical: `P̄`, σ², μ_Dynkin, CLT SD | 0.4400 W, 0.2492 W²h, 41.25 h, 7.21 h | 1e-3 |
| Hub σ²: deviation matrix vs closed form; two-state formula | equal | 1e-10 |
| Monte Carlo L0 mean (20 000 paths, seed 1) vs Dynkin | – | ≤ 1 % |
| Monte Carlo L2 mean (400 paths, h = 10 s) typical / heavy | 39.53 / 21.22 h | 3 % (seed-dependent) |
| Aging calibration `k_H`, `k_L` | 2.697e-3, 1.288e-3 | 1 % |
| A1 SOH at 2 years (quasi-static driver) | 0.834 | ±0.005 |
| Levers (§12.1) | table | ±0.3 pp |

### 15.4 Compute budget (one core, NumPy)

* One L2 LSODA run: 0.1–0.5 s. All 17 §5.4 baseline runs took 9 s [PROTO].
* L2 CTMC Monte Carlo, 400 paths at h = 10 s: ≈ 20 s per profile [PROTO].
* Sobol on the surrogate (51,200 evaluations per scenario): vectorise F2 with 200-point Gauss–Legendre on
  `[z_c, z0]` and Newton iterations for `z_c`; well under a minute per scenario.
* Full-L2 Sobol (N = 256, k = 24: 12,800 runs) with `rk4_batch` at h = 2 s: ≈ 10–20 min per scenario.
* UQ, N = 2048 × ≈ 15 scenarios: ≈ 30 min. Total pipeline < 2 h.

---

## 16. Paper outline (English, LaTeX, `paper/main.tex`)

**Page budget.** The body is ≤ 21 pages, which keeps the total within COMAP's 25 pages:
summary sheet 1 + table of contents 1 + body ≤ 21 + references ≈ 1.5 = 24.5. The AI-use report is appended and
does not count.

| § | Section | Pages | Content | Figures / tables |
|---|---|---|---|---|
| – | Summary sheet | (1) | problem, model in one sentence, 6 headline numbers, top 3 user and OS recommendations | – |
| – | Contents | (1) | | – |
| 1 | Introduction | 0.75 | restatement; approach (hierarchy L0–L3); roadmap | – |
| 2 | Assumptions and notation | 1.0 | A1–A15 with justification; notation table | – |
| 3 | The continuous-time model | 4.5 | 3.1 L0 → L1 (0.75); 3.2 load model: display, SoC/DVFS, radio with nested timers, peripherals (1.25); 3.3 electro-thermal ECM, throttling, brown-out (1.5); 3.4 aging layer and charging (0.5); 3.5 SOC outputs and TTE events (0.5) | `fig_architecture.pdf`; T2 (condensed primitives) |
| 4 | Analytical insights | 2.0 | F1 (energy share = elasticity), F2–F3 (`a` = 0.0082 W⁻¹), F4 (safe burst budget), F5 (hot ambient lengthens gaming TTE), F6–F8 (Dynkin, CLT/IG, variance attribution) | `fig_cell_maps.pdf` |
| 5 | Parameterisation and data | 1.0 | parameter sources, the single calibration, data licences | T1 (≈ 35 rows) |
| 6 | Validation | 2.0 | Bonn hold-outs and slope test; HG2 temperature; V6 band; V3/V5 in brief; where the model fails (GNSS, Android floor) | `fig_validation_bonn.pdf`, `fig_validation_band.pdf` (HG2 as a panel) |
| 7 | Time-to-empty predictions | 2.5 | scenarios × s0 × T with UQ; trajectories; drivers of rapid drain; "surprisingly little" | T3; `fig_trajectories.pdf`; `fig_power_breakdown.pdf`; `fig_tte_s0.pdf` |
| 8 | Uncertainty and sensitivity | 2.5 | violins and variance split; CTMC distribution; tornado; Sobol; structural variants; usage fluctuations | `fig_uq_violins.pdf`, `fig_ctmc_tte.pdf`, `fig_tornado.pdf`, `fig_sobol.pdf`, T4 |
| 9 | Cold, aging and charging | 1.25 | stranded-charge map and O1; SOH under A1–A7; winter shutdowns of aged phones | `fig_stranded_map.pdf`, `fig_aging.pdf` |
| 10 | Recommendations | 1.5 | ranked user levers; OS strategies O1–O9; short user guide | `fig_levers.pdf` |
| 11 | Generalisation | 0.5 | scaling rules, Π groups, device classes | `fig_generalization.pdf` (small) |
| 12 | Strengths, limitations, extensions | 0.75 | | – |
| 13 | Conclusion | 0.25 | | – |
| | **Body total** | **20.5** | 0.5 page reserve | 13 figures, 4 tables |
| – | References | (≈ 1.5) | §17 entries cited in the text | – |
| – | AI use report | (not counted) | `paper/ai_report.tex` | – |

Optional figures (`fig_throttle.pdf`, `fig_cv_scaling.pdf`, `fig_usage_fluct.pdf`, `fig_validation_hg2.pdf`,
`fig_structural.pdf`) are produced but enter the paper only if space allows.

---

## 17. References to cite

Only sources recorded in `docs/research/*.md` are listed. Licence and access notes are in those files. Keys in
**bold** are the minimum set the paper must cite. Press and secondary sources are flagged and used only as cited
facts.

### 17.1 Device power

1. **[CH10]** A. Carroll, G. Heiser. "An Analysis of Power Consumption in a Smartphone." *USENIX ATC 2010*,
   pp. 271–284. https://www.usenix.org/legacy/event/usenix10/tech/full_papers/Carroll.pdf
2. **[CH13]** A. Carroll, G. Heiser. "The Systems Hacker's Guide to the Galaxy: Energy Usage in a Modern
   Smartphone." *APSys 2013*, Singapore. doi:10.1145/2500727.2500734
3. **[HU12]** J. Huang, F. Qian, A. Gerber, Z. M. Mao, S. Sen, O. Spatscheck. "A Close Examination of Performance
   and Power Characteristics of 4G LTE Networks." *MobiSys 2012*, pp. 225–238. doi:10.1145/2307636.2307658
4. **[NA21]** A. Narayanan, X. Zhang, R. Zhu, A. Hassan, S. Jin, et al. "A Variegated Look at 5G in the Wild:
   Performance, Power, and QoE Implications." *ACM SIGCOMM 2021*, pp. 610–625. doi:10.1145/3452296.3472923
5. **[XU20]** D. Xu, A. Zhou, X. Zhang, G. Wang, X. Liu, C. An, Y. Shi, L. Liu, H. Ma. "Understanding Operational
   5G: A First Measurement Study on Its Coverage, Performance and Energy Consumption." *ACM SIGCOMM 2020*.
   doi:10.1145/3387514.3405882
6. **[GU24]** A. Gupta, A. Heidari, A. Kalipattapu, I. K. Jain, D. Bharadia. "3 W's of Smartphone Power
   Consumption: Who, Where and How Much Is Draining My Battery?" *ACM MobiCom 2024* (poster).
   doi:10.1145/3636534.3695905
7. **[DI13]** N. Ding, D. Wagner, X. Chen, A. Pathak, Y. C. Hu, A. Rice. "Characterizing and Modeling the Impact of
   Wireless Signal Strength on Smartphone Battery Drain." *ACM SIGMETRICS 2013*, pp. 29–40.
   doi:10.1145/2494232.2466586 (DOI as given by a BibSonomy record)
8. [CH15] X. Chen, N. Ding, A. Jindal, Y. C. Hu, M. Gupta, R. Vannithamby. "Smartphone Energy Drain in the Wild:
   Analysis and Implications." *ACM SIGMETRICS 2015*, pp. 151–164. doi:10.1145/2745844.2745875 (numbers taken
   from M. Zhan's WPI CS528 presentation, 2016; secondary)
9. **[DZ11]** M. Dong, L. Zhong. "Chameleon: A Color-Adaptive Web Browser for Mobile OLED Displays."
   *MobiSys 2011*; arXiv:1101.1240
10. **[DH21]** P. Dash, Y. C. Hu. "How Much Battery Does Dark Mode Save? An Accurate OLED Display Power Profiler
    for Modern Smartphones." *MobiSys 2021* (values via the Purdue University news release, 28 July 2021;
    secondary)
11. **[CX13]** X. Chen, Y. Chen, Z. Ma, F. C. A. Fernandes. "How Is Energy Consumed in Smartphone Display
    Applications?" *HotMobile 2013*. https://hotmobile.org/2013/papers/full/17.pdf
12. [HE20] C. Herglotz, S. Coulombe, C. Vazquez, A. Vakili, A. Kaup, J.-C. Grenier. "Power Modeling for Video
    Streaming Applications on Mobile Devices." *IEEE Access* 8:70234–70244, 2020. doi:10.1109/ACCESS.2020.2986580
    (CC BY 4.0)
13. **[KW19]** B. Karki, M. Won. "Characterizing Power Consumption of Dual-Frequency GNSS of Smartphone."
    arXiv:1910.13041
14. [SU14] L. Sun, R. K. Sheshadri, W. Zheng, D. Koutsonikolas. "Modeling WiFi Active Power/Energy Consumption in
    Smartphones." *IEEE ICDCS 2014*, pp. 41–51. doi:10.1109/ICDCS.2014.13
15. [BH19] G. Bhat, S. Gumussoy, U. Y. Ogras. "Power and Thermal Analysis of Commercial Mobile Platforms:
    Experiments and Case Studies." *DATE 2019*; arXiv:1904.09814
16. **[CHI15]** V. Chiriac, S. Molloy, J. Anderson, K. Goodson. "A Figure of Merit for Smart Phone Thermal
    Management." *Electronics Cooling*, 30 Nov 2015 (reprinted 19 Sept 2023, CHI23).
    https://www.electronics-cooling.com/?p=19849
17. [AT20] A. Frumusanu. "Samsung Galaxy S20+ & Ultra (Snapdragon & Exynos) Battery Life Preview." AnandTech,
    23 Mar 2020 (figures via 9to5Google, 24 Mar 2020, and Ubergizmo, 25 Mar 2020; secondary, flag in the paper)
18. **[LX]** Linux kernel v6.6: `Documentation/power/energy-model.rst`; device trees `mt8183.dtsi`, `sdm845.dtsi`,
    `sm8250.dtsi`, `rk3399.dtsi`. https://github.com/torvalds/linux/tree/v6.6
19. [AOSP] Android Open Source Project. "Power profiles for Android." https://source.android.com/docs/core/power/values

### 17.2 Cell, ECM and battery models

20. [LIN26] Z. Lin, Y. Ying, H. Li, Y. Ren, T. Liu, P. Hou, H. Huang. "Decoding High-voltage LiCoO2: From
    Degradation to Stabilization Toward Durable Li-ion Batteries." *Advanced Materials* 38 (2026).
    doi:10.1002/adma.202523570
21. **[PIX9]** Google. "Pixel 9, Pixel 9 Pro, Pixel 9 Pro XL Spec Sheet" (2024), footnote 6.
22. [MR24] MacRumors. "iPhone 16 Battery Capacities Revealed." 19 Sep 2024 (press, citing Anatel filings).
23. **[CRM06]** M. Chen, G. A. Rincón-Mora. "Accurate electrical battery model capable of predicting runtime and
    I–V performance." *IEEE Trans. Energy Conversion* 21(2), 504–511 (2006). doi:10.1109/TEC.2006.874229
24. **[AI20]** W. Ai, L. Kraft, J. Sturm, A. Jossen, B. Wu. "Electrochemical Thermal-Mechanical Modelling of Stress
    Inhomogeneity in Lithium-Ion Pouch Cells." *J. Electrochem. Soc.* 167, 013512 (2020). doi:10.1149/2.0122001JES
25. [MQ19] S. G. Marquis, V. Sulzer, R. Timms, C. P. Please, S. J. Chapman. "An Asymptotic Derivation of a Single
    Particle Model with Electrolyte." *J. Electrochem. Soc.* 166, A3693–A3706 (2019). doi:10.1149/2.0341915jes
26. **[SULZER21]** V. Sulzer, S. G. Marquis, R. Timms, M. Robinson, S. J. Chapman. "Python Battery Mathematical
    Modelling (PyBaMM)." *J. Open Research Software* 9(1), 14 (2021). doi:10.5334/jors.309
27. [PLETT-N] G. L. Plett. *ECE4710/5710 Lecture Notes*, Ch. 2 "Equivalent-Circuit Cell Models," University of
    Colorado Colorado Springs. http://mocha-java.uccs.edu/ECE5710/ECE5710-Notes02.pdf
28. [BARC22] S. Barcellona, S. Colnago, G. Dotelli, S. Latorrata, L. Piegari. "Aging effect on the variation of
    Li-ion battery resistance as function of temperature and state of charge." *J. Energy Storage* 50, 104658
    (2022). doi:10.1016/j.est.2022.104658
29. [SUN23] B. Sun, X. Qi, D. Song, H. Ruan. "Review of Low-Temperature Performance, Modeling and Heating for
    Lithium-Ion Batteries." *Energies* 16(20), 7142 (2023). doi:10.3390/en16207142
30. [PINE] ShenZhen SaiJiaoYang Energy & Science Technology Co. "Lithium battery specification, model
    QZ01-396172, 2750 mAh" (PinePhone battery), rev. A0.
    https://xnux.eu/devices/feature/bat-pp/PinePhone_Battery_model_QZ01-396172-2750.pdf
31. [BU802a] Battery University. "BU-802a: How does Rising Internal Resistance affect Performance?"
    https://www.batteryuniversity.com/article/bu-802a-how-does-rising-internal-resistance-affect-performance/
32. [QC-FG] Qualcomm / AOSP kernel msm. `qpnp-fg.txt` device-tree binding (`qcom,fg-cutoff-voltage-mv`).
    https://android.googlesource.com/kernel/msm/+/973f4134d9deb396415846f902848f0a32cb4cfa/Documentation/devicetree/bindings/power/qpnp-fg.txt
33. [QC-SMB] Qualcomm / AOSP kernel msm. `smb1360-charger-fg.txt`.
    https://android.googlesource.com/kernel/msm/+/973f4134d9deb396415846f902848f0a32cb4cfa/Documentation/devicetree/bindings/power/smb1360-charger-fg.txt
34. [OVEJ19] V. J. Ovejas, A. Cuadras. "Effects of cycling on lithium-ion battery hysteresis and overvoltage."
    *Scientific Reports* 9, 14875 (2019). doi:10.1038/s41598-019-51474-5
35. **[BERN85]** D. Bernardi, E. Pawlikowski, J. Newman. "A General Energy Balance for Battery Systems."
    *J. Electrochem. Soc.* 132, 5–12 (1985). doi:10.1149/1.2113792
36. [VISW10] V. V. Viswanathan, D. Choi, D. Wang, W. Xu, S. Towne, R. E. Williford, et al. "Effect of entropy
    change of lithium intercalation in cathodes and anodes on Li-ion battery thermal management." *J. Power
    Sources* 195, 3720–3729 (2010). doi:10.1016/j.jpowsour.2009.11.103
37. **[JH09]** M. R. Jongerden, B. R. Haverkort. "Which battery model to use?" *IET Software* 3, 445–457 (2009).
    doi:10.1049/iet-sen.2009.0001
38. [MM93] J. F. Manwell, J. G. McGowan. "Lead acid battery storage model for hybrid energy systems." *Solar
    Energy* 50, 399–405 (1993). doi:10.1016/0038-092X(93)90060-2
39. [RV01] D. Rakhmatov, S. Vrudhula. "An analytical high-level battery model for use in energy management of
    portable electronic systems." *Proc. IEEE/ACM ICCAD 2001*, 488–493. doi:10.1109/ICCAD.2001.968687

### 17.3 Thermal behaviour, cold, aging and charging

40. **[KANG19]** S. Kang, H. Choi, S. Park, C. Park, J. Lee, U. Lee, S.-J. Lee. "Fire in Your Hands: Understanding
    Thermal Behavior of Smartphones." *ACM MobiCom 2019*. doi:10.1145/3300061.3300128
41. [PU23] H. Pu, X. Wu. "Investigating Skin Temperature-Based Overheating in mmWave Smartphones: Power and
    Thermal Models for Optimal Non-Throttling Performance." arXiv:2304.06512 (2023)
42. **[HE19]** L. He, Y. Lee, E. Kim, K. G. Shin. "Environment-Aware Estimation of Battery State-of-Charge for
    Mobile Devices." *ACM/IEEE ICCPS 2019*. https://rtcl.eecs.umich.edu/rtclweb/assets/publications/2019/ICCPS-liang19.pdf
43. **[LEE20]** Y. Lee, L. He, K. G. Shin. "Causes and Fixes of Unexpected Phone Shutoffs." *ACM MobiSys 2020*.
    doi:10.1145/3386901.3389024
44. [MAL99] H. Maleki, S. A. Hallaj, J. R. Selman, R. B. Dinwiddie, H. Wang. "Thermal Properties of Lithium-Ion
    Battery and Components." *J. Electrochem. Soc.* 146(3), 947–954 (1999). doi:10.1149/1.1391704
45. **[SCH14]** J. Schmalstieg, S. Käbitz, M. Ecker, D. U. Sauer. "A holistic aging model for Li(NiMnCo)O2 based
    18650 lithium-ion batteries." *J. Power Sources* 257, 325–334 (2014). doi:10.1016/j.jpowsour.2014.02.012
    (equations via CA21)
46. [CA21] M. Clarke, J. J. Alonso. "Lithium–Ion Battery Modeling for Aerospace Applications." *Journal of
    Aircraft* (2021). doi:10.2514/1.C036209
47. **[SCM18]** M. Schimpe, M. E. von Kuepach, M. Naumann, H. C. Hesse, K. Smith, A. Jossen. "Comprehensive
    Modeling of Temperature-Dependent Degradation Mechanisms in Lithium Iron Phosphate Batteries."
    *J. Electrochem. Soc.* 165(2), A181–A193 (2018). doi:10.1149/2.1181714jes
48. **[KEIL16]** P. Keil, S. F. Schuster, J. Wilhelm, J. Travi, A. Hauser, R. C. Karl, A. Jossen. "Calendar Aging of
    Lithium-Ion Batteries: I. Impact of the Graphite Anode on Capacity Fade." *J. Electrochem. Soc.* 163(9),
    A1872–A1880 (2016). doi:10.1149/2.0411609jes
49. [WU17] Y. Wu, P. Keil, S. F. Schuster, A. Jossen. "Impact of Temperature and Discharge Rate on the Aging of a
    LiCoO2/LiNi0.8Co0.15Al0.05O2 Lithium-Ion Pouch Cell." *J. Electrochem. Soc.* 164(7), A1438 (2017).
    doi:10.1149/2.0401707jes
50. [NEL18] K. J. Nelson, J. E. Harlow, J. R. Dahn. "A Comparison of NMC/Graphite Pouch Cells and Commercially
    Available LiCoO2/Graphite Pouch Cells Tested to High Potential." *J. Electrochem. Soc.* 165(3), A456 (2018).
    doi:10.1149/2.0041803jes
51. [WAL14] T. Waldmann et al. "Temperature dependent ageing mechanisms in Lithium-ion batteries – A Post-Mortem
    study." *J. Power Sources* 262, 129–135 (2014) (abstract only)
52. [ZIL19] I. Zilberman, S. Ludwig, A. Jossen. "Cell-to-cell variation of calendar aging and reversible
    self-discharge in 18650 nickel-rich, silicon–graphite lithium-ion cells." *J. Energy Storage* 26, 100900
    (2019). doi:10.1016/j.est.2019.100900
53. **[AP-BATT]** Apple Support. "iPhone battery and performance." Article 101575.
    https://support.apple.com/en-us/101575
54. [AP-TEMP] Apple Support. "If your iPhone or iPad gets too hot or too cold." Article 118431.
    https://support.apple.com/en-us/118431
55. [AP-OBC] Apple Support. "About Charge Limit and Optimized Battery Charging on iPhone." Article 108055.
    https://support.apple.com/108055
56. [AP-WHY] Apple. "Batteries – Why Lithium-ion?" https://www.apple.com/batteries/why-lithium-ion/
57. **[GOOG]** Google Pixel Phone Help. "Understand your Pixel battery."
    https://support.google.com/pixelphone/answer/15738128
58. [SAM] Samsung Singapore Support. "Galaxy Battery protection feature in One UI 6.1." Updated 11 Apr 2024.
    https://www.samsung.com/sg/support/mobile-devices/galaxy-battery-protection-feature-in-one-ui-6-1/
59. **[TI-JEITA]** J. Qian. "Li-ion battery-charger solutions for JEITA compliance." *Texas Instruments Analog
    Applications Journal*, 1Q 2010, SLYT365. https://www.ti.com/lit/pdf/slyt365
60. [UL-TS] F. Chiang (UL Taiwan). "IEC 62368-1 Technical Brief: Touch Temperature Limits." 22 Oct 2010.
    https://japan.ul.com/wp-content/uploads/sites/27/2014/06/1_techbrief_touchtemp.pdf
61. [AOSP-TH] Android Open Source Project. "Thermal mitigation."
    https://source.android.com/docs/core/power/thermal-mitigation

### 17.4 Validation data and published runtimes

62. **[BONN26]** F. Schweizer, J. Yu, E. Mille, L. M. Reimer, M. Kapsecker, J. Klinker, S. M. Jonas (2026). "Battery
    Impact of Individual Hardware Sensors in Modern iPhone Devices" [dataset]. bonndata. doi:10.60507/FK2/LSJGWW.
    CC BY 4.0
63. **[SCH26]** F. Schweizer, J. Yu, E. Mille, L. M. Reimer, M. Kapsecker, J. Klinker, S. Jonas (2026). "Smartphone
    Sensor Battery Consumption: A Standardized and Reproducible Test Protocol." *Sensors* 26(10), 2923.
    doi:10.3390/s26102923
64. **[KOL20]** P. Kollmeyer, C. Vidal, M. Naguib, M. Skells (2020). "LG 18650HG2 Li-ion Battery Data and Example
    Deep Neural Network xEV SOC Estimator Script." Mendeley Data v3. doi:10.17632/cp3473x7xv.3. CC BY 4.0
65. **[NASA07]** B. Saha, K. Goebel (2007). "Battery Data Set." NASA Prognostics Data Repository, NASA Ames Research
    Center. https://phm-datasets.s3.amazonaws.com/NASA/5.+Battery+Data+Set.zip (licence not specified; public
    access; acknowledgement requested; document this in the paper)
66. [OX17] C. Birkl, D. Howey (2017). "Oxford Battery Degradation Dataset 1." University of Oxford.
    doi:10.5287/bodleian:KO2kdmYGg. ODbL v1.0 / DbCL v1.0
67. **[APL14P]** Apple Inc. "iPhone 14 Pro — Technical Specifications." https://support.apple.com/en-us/111849
68. [APLBAT] Apple Inc. "iPhone — Battery Test Information." https://www.apple.com/iphone/battery.html
69. [N9M22] 9to5Mac (11 Sep 2022). Report on iPhone 14 battery sizes (Chemtrec filing).
    https://9to5mac.com/2022/09/11/iphone-14-battery-sizes/ (press)
70. [NBC22] Notebookcheck (2022). "Apple iPhone 14 Pro Max review: A gigantic brawny smartphone."
    https://www.notebookcheck.net/Apple-iPhone-14-Pro-Max-review-A-gigantic-brawny-smartphone.659750.0.html
71. [NBC23] Notebookcheck (2023). "Google Pixel 8 smartphone review: Compact and with 7 years of updates."
    https://www.notebookcheck.net/Google-Pixel-8-smartphone-review-Compact-and-with-7-years-of-updates.768436.0.html
72. [NBC24] Notebookcheck (2024). "Samsung Galaxy S24 review – The best small Android smartphone gets cheaper but
    not better." https://www.notebookcheck.net/Samsung-Galaxy-S24-review-The-best-small-Android-smartphone-gets-cheaper-but-not-better.809474.0.html
73. [LMAG22] Laptop Mag (2022). "iPhone 14 Pro battery life results are in — and we're not impressed."
    https://www.laptopmag.com/news/iphone-14-pro-battery-life-results-are-in-and-were-not-impressed
74. **[GPX]** Google. "Pixel phone hardware tech specs." https://support.google.com/pixelphone/answer/7158570

**Data licences for the paper's data statement.** Open-licensed data used: BONN26 (CC BY 4.0), KOL20 (CC BY 4.0),
OX17 (ODbL), Ai2020 parameters via PyBaMM (CC BY; BSD-3 code). NASA07 has no stated licence (public access,
acknowledgement requested); we document it and use it only for a secondary check. Published measurements and
specifications (papers, Apple, Google, Notebookcheck, Laptop Mag) are used as cited facts, as the problem allows.

---

## Appendix A. Derivations

**A.1 Current root (E2).** `R0 I² − Ũ I + P = 0` has roots `I± = (Ũ ± √D)/(2R0)`. The physical root is the smaller
current (higher terminal voltage), `I− = (Ũ − √D)/(2R0)`. Multiplying numerator and denominator by `(Ũ + √D)` gives
`I = 2P/(Ũ + √D)`, which has no cancellation as `P → 0`. Then `V = P/I = (Ũ + √D)/2`.

**A.2 L1 quadrature (F2).** With `R_Σ` and `v1 = R1 I` at quasi-equilibrium, `1/I = V_qs/P`, so
`dt = 3600 Q_max dz / I = 3600 Q_max V_qs dz / P`. Integrate from `z_c` to `z0`. At the cut-off `V = V_cut`, so
`I = P/V_cut` and `U(z_c) = V_cut + P R_Σ/V_cut`.

**A.3 Radio time-average with nested timers (E14).** Let bursts form a Poisson process `N` of rate λ, and let the
radio state at time t depend only on the age `A(t)` of the last burst. The radio is active if `A < T_act`, in the tail
if `T_act ≤ A < T_tail`, and idle otherwise. For a Poisson process, `Pr(A(t) ≥ s) = Pr(no arrival in (t − s, t]) =
e^(−λs)`. The time-average power is therefore
`P_idle + (P_tail − P_idle)(1 − e^(−λT_tail)) + (P_act − P_tail)(1 − e^(−λT_act))`. A promotion occurs when a burst
arrives while `A ≥ T_tail`. By PASTA its rate is `λ·e^(−λT_tail)`, which adds `λe^(−λT_tail)E_pro`. Bursts are
assumed short compared with `T_act`; their transfer energy is the per-bit term.

**A.4 First-order constant `a` (F3).** Expand `V_qs = U − R_Σ P/U + O(P²)` and the cut-off shift
`Δz_c = P R_Σ(z_e)/(V_cut U′(z_e))`. Then `∫_{z_c(P)} V_qs dz = ∫_{z_e} U dz − P∫R_Σ/U dz − U(z_e)Δz_c + O(P²)`, and
`U(z_e) = V_cut`.

**A.5 Brown-out threshold (F4).** In E10 set `V_pk = V_sys`. Then `Ũ − R0 I_pk = V_sys` with `I_pk = P_pk/V_sys`,
so `U(z*) = V_sys + v1 + R0(z*) P_pk/V_sys`. Use `v1 ≈ R1 P/U` (RC at equilibrium under the mean load).

**A.6 Dynkin identity (F6).** `h = Dm(P − P̄1)` solves the Poisson equation `−Λh = P − P̄1` with `πᵀh = 0`. Apply
Dynkin's formula to `h(X_t)` up to the bounded stopping time `TTE ∧ t` and let `t → ∞`:
`E[h(X_TTE)] − h(x0) = E∫_0^TTE (Λh)(X_s) ds = E∫_0^TTE (P̄ − P_{X_s}) ds = P̄·E[TTE] − E_u`. This is exact in L0,
because `∫_0^TTE P ds = E_u`. The depletion state is distributed approximately as `π_x P_x/P̄` (a level crossing of a
cumulative process happens in state x at a rate proportional to its power).

**A.7 Hub variance (F8).** One renewal cycle is an idle sojourn `D0` followed by a session `D_J`. Its centred reward
is `Y = (P0 − P̄)D0 + (P_J − P̄)D_J`, with `E[Y] = 0`, i.e. `Σ_j p_j(P_j − P̄)d_j = −(P0 − P̄)/λ_on`. With
`E[D²] = (1 + CV²)·E[D]²`:
`E[Y²] = (P0 − P̄)²(1 + CV0²)/λ_on² + Σ_j p_j(1 + CV_j²)d_j²(P_j − P̄)² − 2(P0 − P̄)²/λ_on²`. The renewal-reward CLT
gives `σ² = E[Y²]/E[L]`. For exponential sojourns (CV = 1) the idle terms cancel, which gives F8. The general form
is F9.

**A.8 Inverse Gaussian (F7).** The first passage of a Brownian motion with drift P̄ > 0 and variance rate σ² to the
level `E_u` is `IG(E_u/P̄, E_u²/σ²)`. Replacing the mean by μ_Dynkin corrects the initial-state and depletion-state
bias at O(1).

**A.9 Squared-loss aging form (E17–E19).** For `q = k·√τ`, where τ is a stress-weighted exposure,
`q² = k²τ ⇒ d(q²)/dt = k²·dτ/dt`. With time-varying stress, the rate-based integral `dq/dt = k(t)²/(2q)` is
equivalent to `d(q²)/dt = k(t)²`. The latter is linear, exact over piecewise-constant segments, and well defined at
`q = 0`.

**A.10 Energy closure (test T2).** `U·I = V·I + R0 I² + v1·I` and `I = v1/R1 + C1 dv1/dt`. Hence
`v1·I = v1²/R1 + d(½C1v1²)/dt` and `∫U I dt = ∫P dt + ∫R0 I² dt + ∫v1²/R1 dt + ½C1 v1(t_end)²`, where
`∫U I dt = −3600·Q_max·∫U dz` along the path. The entropic term exchanges heat with the chemistry, not electrical
energy, so it does not enter.

---

## Appendix B. Review of the two proposals

Scores from 1 (poor) to 10 (excellent).

| Criterion | PF (physics-first) | PS (parsimony + stochastic) |
|---|---|---|
| Physical grounding | 9 | 7 |
| Covers every problem requirement | 8 | 9 |
| Traceability of parameters to cited sources | 8 | 7 |
| Implementable in Python/SciPy within a day | 6 | 8 |
| Insight and explanatory power for the paper | 8 | 9 |
| Novelty | 7 | 8 |
| **Total (of 60)** | **46** | **48** |

**PF strengths.** Fully coupled ECM, a mechanistic load model, two voltage thresholds with a burst envelope that
explains cold and aged shutdowns, a stranded-charge map, a safe-burst budget for an OS governor, clean numerics with
energy-closure tests, and pre-registered blind Bonn predictions.

**PF errors and weaknesses** (all fixed above).

1. **Radio formula has no connected-active state.** At λ ≫ 1/T_act (gaming at 20 packets/s) the LTE radio sits at
   the 0.15 W tail level, so gaming over LTE costs about the same as over Wi-Fi. GU24 shows ≈ 0.73 W of cellular
   rail in a game. Fixed with the nested C-DRX timer (E14, Appendix A.3).
2. **Standby primitive SB routes traffic over cellular while Wi-Fi is associated** (`w = 1`, `ρ = cell`). This is
   inconsistent. Split into SB-W and SB-L.
3. **Brightness-law claim.** "γ_b = 2 fits DH21 at both ends" holds only near b = 0.3. On the S4 load the model
   saves ≈ 17 % at b = 0.5 against DH21's 3–9 % at 30–50 %. The power law also gives zero content power at minimum
   brightness, whereas CH13 has g(min) = 0.16. Kept γ_b = 2 as a compromise; the claim is corrected, with
   [PROTO] savings of 5–13 % at 30–50 % and 33 % at 100 %.
4. **Cycle aging (E14).** It sums the SEI and plating mechanisms inside one square, which creates cross terms. It
   drives both by total throughput with `c_ch = 0` on discharge, which leaves a spurious plating weight of
   `w_L·e^(−β_L) ≈ 2 %` during discharge. It also needs a `q(0) = 10⁻⁴` hack. SCM18 Eq. 23 uses separate `√Q_tot`
   and `√Q_ch` terms. Fixed by the squared-loss sub-states E17–E19.
5. **ΔP_cam = 1.4 W** from CH13 whole-device minus screen-on idle overstates the modern increment. CX13 subtraction
   gives 1.03–1.28 W, and the Bonn increments are 1.04–1.22 W. Set to 1.15 W.
6. **Internal inconsistency in F6.** The analytic CTMC mean (29.7 h) does not equal the stated formula
   `E_u/P̄ = 17.79/0.590 = 30.2 h`.
7. **"Typical day" S12 has 9.5 h of screen-on time**, including 1 h of gaming, 1.5 h of navigation and 2.5 h of
   video. That is a demanding day; it matched the Pixel 8a's 0.56 W only through low component powers. Renamed
   "scheduled full day". The CTMC typical day (6.5 h) is now the reference.
8. **Recommendation ranking** used mean demand while the headline TTE is a non-smooth L2 event. This is formalised
   in §12.1.
9. The entropic value at z = 0 (−0.449 mV/K) is cited to cell_ecm Appendix A2, which starts at z = 0.05. It does
   exist in the [RUN] JSON (`pybamm_ocv_fit.json`), which is now cited directly.

**PS strengths.** The model ladder (M0–M4) with one mechanism per level; energy share = elasticity; the exact
first-order constant `a`; the Dynkin mean; the CLT/IG law; hub-and-spoke variance attribution ("gaming is 5 % of
sessions but most of the variance"); usage elasticities; the semi-Markov extension; aging results for charging
strategies; and the widest validation band (10 published runtimes).

**PS errors and weaknesses** (all fixed above).

1. **The floor was tuned on the primary validation set.** `ΔP_awake = 0.08 W` was "informed by BONN26", and the
   idle run was then reported as a validation (+2.1 %); the blind error was not reported. It is now an explicit
   [CAL], with the blind prior reported (−5.1 %) and an optional ABC posterior.
2. **One threshold for everything.** `V_cut` is applied to the instantaneous burst voltage. There is no brown-out
   threshold and no sub-second burst mechanism, so the 42 % shutdown rests entirely on a synthetic 30 s / 4 W burst
   train. Now there are two thresholds plus an envelope, and explicit long bursts stay in `u(t)`.
3. **Sign clash.** `A_X(T)` (M3) grows in the cold, while `A(T;E)` (M4) grows in the heat, and `B` uses
   `1/A(T;E_L)` for cold acceleration. This is error-prone. Replaced by the named factors `a_R` and `a_deg`.
4. **Cycle aging.** One state with `½(A_H + plating)` inside the square, driven by discharge FCE, although plating
   depends on charge current. The prototype (`aging.py`) adds an ad-hoc "excess plating on charge" term that is not
   in the spec. `c_ref = 0.5C` is inconsistent with SCM18's 1C reference. Replaced by E17–E19.
5. **Lumped per-activity network and SoC constants** (`P_net,act` table, `P_SoC,dem`). They are mostly
   [ASSUMPTION]s and do not respond to throughput, burst rate, refresh or DVFS. That weakens the "network activity"
   and "processor load" requirements. Replaced by the mechanistic PF loads.
6. **Prototype/spec mismatch.** `ctmc.py` uses `λ_on = 3 h⁻¹` and Wi-Fi/4G means of 3 h / 2 h; the spec states
   4 h⁻¹ and 2.5 h (`ctmc2.py` matches the spec).
7. **A lower bound used as a point target.** Pixel 8a "at least 31 h" was scored as +24 %. It is now treated as a
   bound.
8. **Rated energy instead of measured capacity.** Validation used Apple's rated 12.38 Wh although SCH26 measured
   3.32 Ah (12.57 Wh usable on our OCV), a 1.5 % bias. Now uses the measured capacity.
9. **Affine brightness law.** It caps the dark-mode saving at ≈ 20 % at full brightness against DH21's 39–47 %.
   It was acknowledged only as a variant; γ_b is now in UQ.

**Shared omission (both).** Both used the coulombic map as "SOC". The problem glossary defines SOC as an energy
fraction. The final spec reports both: D1 for initial conditions and `s_off`, and D2 (energy fraction) as the
glossary SOC. The Bonn slope test V1.4 is expected to favour D2.
