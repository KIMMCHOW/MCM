# Model specification, proposal P-S: a parsimonious model hierarchy with stochastic usage

2026 MCM Problem A, *Modeling Smartphone Battery Drain*. Design proposal written 2026-11-17.

**Status.** This document is implementation-ready. Every parameter has a value, a unit, a source and an
uncertainty range. The numbers marked **[PROTO]** were computed in this session by the prototype scripts in
`docs/design/prototype_parsimony/*.py`, which implement the equations below. Use them as regression targets: the
production code in `src/phonebatt/` must reproduce them within the stated tolerance (§9.3).

**Source keys.** Keys such as `CH13`, `GU24` and `KEIL16` refer to the reference lists of the four research notes:

* `docs/research/device_power.md` (dp),
* `docs/research/cell_ecm.md` (ce),
* `docs/research/thermal_aging.md` (ta),
* `docs/research/validation_data.md` (vd).

`[RUN]` marks a value computed in an earlier session from open data or PyBaMM, as recorded in those notes.

**Confidence tags.**

* **H**: primary or official source.
* **M**: figure read-off, secondary source, or physics-model output.
* **L**: our assumption or derivation.

**DERIVED** marks our own arithmetic on cited numbers. **ASSUMPTION** marks a modelling choice with no direct
source; each one has a UQ range and appears in the sensitivity analysis.

**Units and signs.** SI units are used unless a table states otherwise; energy is in Wh where convenient.
Current I > 0 on discharge. Inside equations T is in kelvin. R_g = 8.314 J mol⁻¹ K⁻¹ and T_ref = 298.15 K.

---

## 0. Design philosophy and headline results

**Philosophy.** We use the simplest model that answers each question, and add physics one effect at a time.

1. **M0, the energy bucket.** `dE/dt = −P(u(t))` gives the closed form `TTE = E_u/P̄`. Every TTE elasticity with
   respect to a component power is then minus that component's energy share.
2. Each extension M1 → M4 adds **one** physical mechanism:
   * M1: OCV and series resistance, so a constant-power load draws more current as voltage falls;
   * M2: RC polarisation and recovery, or KiBaM bound charge as a variant;
   * M3: a lumped thermal node with Arrhenius resistance and thermal throttling;
   * M4: capacity fade and resistance growth (SOH).

   We quantify the extra TTE each mechanism produces. This is the *structural sensitivity*.
3. **Usage is a continuous-time Markov chain (CTMC)** over activity states and connectivity states. The battery
   state then follows a **piecewise-deterministic Markov process (PDMP)**. We derive closed forms for the fluid
   mean, the exact mean (via a Dynkin/Poisson-equation identity), the CLT variance and an inverse-Gaussian TTE
   law, and check them against Monte Carlo.
4. **No black-box fitting.** Parameters come from published measurements. The one calibration is the
   screen-on idle floor against the Bonn iPhone idle runs; it is flagged as a calibration, and everything else is
   validated out of sample (§6).

**Headline results [PROTO]** (reference phone RP-63, §3; 25 °C, new battery, from 100 %).

| Question | Result |
|---|---|
| Usable energy | E_n = 18.19 Wh (0→1 OCV integral); E_u = 17.79 Wh to the 3.40 V cut-off at vanishing load; 0.978 of nameplate |
| Constant-use TTE | standby (4G) 340 h; audio over Wi-Fi ≈ 96 h (M0); video Wi-Fi 28.9 h; video 5G 19.6 h; social Wi-Fi 15.6 h; social 4G 12.8 h; navigation 11.2 h; video call 4G 6.3 h; 3D game 5.9 h (throttled) |
| Typical-day CTMC | P̄ = 0.440 W; mean TTE 41.2 h (exact Dynkin formula; Monte Carlo 41.2 h); SD 6.6–6.9 h (CV 0.17) with M0; M3 Monte Carlo mean 39.9 h, 90 % interval 29.1–51.0 h |
| Heavy-day CTMC | P̄ = 0.809 W; M3 Monte Carlo mean 21.6 h, 90 % interval 14.5–29.7 h |
| Largest drains | 3D gaming (5 % of sessions; 26 % of typical-day energy and **72 % of TTE variance**); cellular video calls; weak signal (−14 % TTE when social browsing at −110 dBm); 5G instead of Wi-Fi for streaming (−32 %) |
| Surprisingly small | battery internal loss and rate-capacity effects at 25 °C (TTE −0.5 % to −3 %); KiBaM bound charge (< 0.5 %); GNSS (1.9 % of typical-day energy); Bluetooth (1.3 %); dark mode at 40 % brightness (+3.4 % TTE); halving background syncs (+0.6 %); a steady 1.75 W load at −10 °C with self-heating (−2 %) |
| Counter-intuitive | in gaming, hot ambient *lengthens* TTE because of throttling: 5.9 h at 25 °C vs 11.1 h at 35 °C, at −43 % performance |
| Cold, aged and bursty compound | 4 W, 30 s bursts every 10 min on an 80 %-SOH battery at −20 °C shut the phone down at **42 % displayed SOC**; on a new battery at −10 °C, at 7 % |

---

## 1. Model structure

### 1.1 State variables

| Symbol | Meaning | Unit | Range | Level | Time scale |
|---|---|---|---|---|---|
| E | remaining usable energy (M0 only) | Wh | [0, E_u] | M0 | hours |
| z | thermodynamic state of charge: fraction of the 0.2C capacity Q between OCV = 3.0 V (z = 0) and 4.45 V (z = 1) | – | [0, 1] | M1+ | hours |
| v₁ | RC polarisation voltage | V | ≥ 0 | M2+ | τ₁ = 30 s (25 °C) to 135 s (−20 °C) |
| δ | KiBaM height difference (variant M2-K, replaces v₁) | Ah | ≥ 0 | M2-K | 1/k′ ≈ 18 s |
| T | lumped phone/battery temperature | K | – | M3+ | τ_th = R_th·C_th = 450 s |
| q_cal, q_cyc | calendar and cycle capacity-loss fractions | – | [0, 0.5] | M4 | days to years |
| X(t) | usage state (activity × connectivity), a CTMC on a finite set 𝒮 | – | 𝒮 | all | minutes |

Derived: SOH = 1 − q_cal − q_cyc; Q_eff = Q₀·SOH; resistance growth r = γ_R·(q_cal + q_cyc).

### 1.2 Inputs: the usage vector u(t)

An **activity** is a preset of the continuous inputs below (§2.7). The CTMC switches between presets. Scenarios
may also override any component directly.

| Input | Symbol | Domain | Default (RP-63) |
|---|---|---|---|
| screen on/off | s_on | {0, 1} | per activity |
| brightness slider fraction | b | [0, 1] | 0.4 indoor (DH21: auto-brightness typically 30–40 %); 1.0 outdoor (non-HBM max) |
| linear-light average picture level | APL | [0, 1] | per activity (§2.2) |
| refresh rate | f_r | {60, 120} Hz | 120 while scrolling (social), else 60 (LTPO policy, ASSUMPTION) |
| SoC demand power | P_SoC,dem | W | per activity (§2.3) |
| radio access technology | c | {air, wifi, 4g, 5g} | from the connectivity CTMC |
| radio activity power at reference signal | P_net,act | W | per activity × RAT (§2.4) |
| cellular signal | RSRP | dBm | −85 (reference) |
| background wake-up rate | λ_bg | s⁻¹ | 1/60 |
| GNSS fix rate | ν_fix | Hz | 1 in navigation, else 0 |
| camera pipeline on | s_cam | {0, rear, front} | per activity |
| ambient temperature | T_amb | °C | 25 |
| charger | P_chg | W | 0 in TTE runs; daily profile in M4 |

### 1.3 Outputs

* SOC(t) in three forms (§1.6): z(t), the displayed SOC s(t) and the state of energy SOE(t).
* Terminal voltage V(t), current I(t), temperature T(t) and skin temperature T_skin(t).
* Throttled power P(t) and performance ratio ρ(t).
* **TTE** (§5.1); the displayed SOC at shutdown, s_off; the TTE distribution for stochastic usage.
* SOH trajectory over months (M4).

### 1.4 The hierarchy

All levels share the load model `P(t) = P_tot(u(t), T)` of §2. Each level adds the terms marked ▶.

**M0, energy bucket.**

```
dE/dt = −P(t) / 3600                     [E in Wh, P in W, t in s]
E(0)  = SOE₀ · E_u,       E_u = Q₀ ∫_{z_c0}^{1} U(z) dz          (Wh, Q₀ in Ah)
TTE   = inf{t : E(t) = 0}
```

* `z_c0 = U⁻¹(V_cut) = 0.0263`, the cut-off SOC at vanishing load [PROTO].
* E_u = 17.79 Wh [PROTO].
* M0 makes the battery an ideal energy store, so the usable energy does not depend on load. This is supported
  at phone loads by the Bonn additivity test: residual ≤ 0.15 % over a 6.5× power range (vd, [RUN]).

**M1, OCV + series resistance, constant-power load.**

```
dz/dt = − I / (3600 Q_eff)
I     = [U(z) − sqrt(U(z)² − 4 R_s P)] / (2 R_s),     R_s = R₀(z,T,r) + R₁(T,r)   (DC, quasi-steady)
V     = U(z) − I R_s
```

* ▶ **New physics:** at low SOC the voltage sag raises the current, and the cut-off V ≥ V_cut strands charge.
* R₁ is added to R_s because, for loads that are constant over many τ₁ (τ₁ = 30 s ≪ hours), the RC branch is
  at steady state.

**M2, RC polarisation and recovery (Thevenin 1RC).**

```
dv₁/dt = − v₁/τ₁(T) + I·R₁(T,r)/τ₁(T)
I      = [U_e − sqrt(U_e² − 4 R₀ P)] / (2 R₀),      U_e = U(z,T) − v₁
V      = U_e − I R₀
```

▶ **New physics:** load transients. A burst shorter than τ₁ sees only R₀. After a shutdown, V recovers by v₁
within about 3τ₁ (the "phone comes back on at 1–3 %" effect).

*Variant M2-K (kinetic battery model, KiBaM).* This replaces v₁ by the KiBaM height difference δ:

```
dδ/dt = I/(3600 c) − k′ δ,          z_surf = z − (1 − c) δ / Q_eff,      U evaluated at z_surf
```

* c = 0.166 (JH09).
* The combination κ_K = c·k′/(1 − c) is the only identifiable one. κ_K = 40 h⁻¹ (33–47), DERIVED by attributing
  *all* of the PyBaMM DFN rate loss, Q(1C)/Q(0.1C) = 0.978 (ce, [RUN] Ai2020; Marquis 0.981), to bound charge,
  so it is an upper bound on the effect. Then k′ = 0.0558 s⁻¹ (1/k′ = 18 s).
* Arrhenius: κ_K(T) = κ_K·A(T; E_a,K) with E_a,K = 35 kJ/mol (25–55; PyBaMM diffusion and exchange-current E_a,
  ce). Used only in the structural-sensitivity and HG2 cold-validation runs.

**M3, lumped thermal node, Arrhenius kinetics, throttling.**

```
C_th dT/dt = P + I² R₀ + v₁²/R₁ − I·T·∂U/∂T(z) − (T − T_amb)/R_th
U(z,T)     = U(z) + (T − T_ref)·∂U/∂T(z)
R₀(z,T,r)  = R₀,ref · [1 + A_R e^{−z/z₀}]/[1 + A_R e^{−0.5/z₀}] · A_R0(T) · (1 + r)
R₁(T,r)    = R₁,ref · A_R1(T) · (1 + r),        τ₁(T) = τ₁,ref · A_τ(T)
A_X(T)     = exp[(E_a,X/R_g)(1/T − 1/T_ref)]
T_skin     = T_amb + (T − T_amb)·R_skin/R_th
P_SoC      = P_SoC,fl + (P_SoC,dem − P_SoC,fl)·φ,     φ = 1/(1 + exp(−(T_thr − T_skin)/w))   (if P_SoC,dem > P_SoC,fl)
```

* ▶ **New physics:** cold raises resistance; self-heating recovers part of it; sustained power is thermally
  capped.
* The heat source is the whole device power P, because nearly all of it dissipates in the phone (ta §2.1), plus
  the battery's irreversible and reversible (Bernardi) heat.

**M4, aging (slow states; quasi-static within one discharge).**

```
Q_eff = Q₀ (1 − q),     q = q_cal + q_cyc,     r = γ_R q
dq_cal/dt_d = [k_cal · A(T; E_cal) · g(z)]² / (2 q_cal)                 (t_d in days)
dq_cyc/dt_d = [k_cyc · B(T, c_ch)]² · (dN_FCE/dt_d) / (2 q_cyc)          (N_FCE = discharge Ah / Q₀)
A(T;E) = exp[−(E/R_g)(1/T − 1/T_ref)]
B(T,c) = ½·[A(T;E_H) + exp(β_L(c − c_ref)) / A(T;E_L)]                   (B = 1 at 25 °C and c = c_ref)
```

* ▶ **New physics:** capacity fade and power fade.
* The rate form is from ta §2.4(b)–(c), with the SCM18 mechanism split: heat-driven SEI growth (E_H) and
  cold/fast-charge plating (E_L, an inverse Arrhenius term).
* Under constant conditions it integrates to q = k·A·g·√t.
* Initial value q(0) = 10⁻⁶ to remove the singularity.

### 1.5 Load–battery coupling, feasibility and shutdown

The PMIC makes the phone a **constant-power load** at the battery terminals (dp §2.1). All powers in §2 are
battery-terminal values, with the PMIC efficiency η ≈ 0.83 (CH13) already folded in.

* **Current.** I is the smaller root of `R I² − U_e I + P = 0`, with (R, U_e) = (R_s, U) in M1 and (R₀, U − v₁)
  in M2+.
* **Feasibility.** `P ≤ P_max,sol = U_e²/(4R)`. In practice the cut-off binds first:
  `P ≤ P_max,cut = V_cut·(U_e − V_cut)/R` (ta T3).
* **Shutdown** is the first time that either `V(t) ≤ V_cut` or `U_e² − 4RP ≤ 0`. Below V_cut, the phone's fuel
  gauge reports 0 and the device powers off (QC-FG, QC-SMB, LEE20).

### 1.6 Three SOC definitions

| Symbol | Definition | Use |
|---|---|---|
| z | coulombic, thermodynamic (state variable) | physics |
| s (displayed SOC) | s = (z − z_ref)/(1 − z_ref), with z_ref = z_c(P_ref = 0.5 W, 25 °C, new) = 0.028 [PROTO] | what the user sees; s_off = s at shutdown ("unexpected shutoff SOC") |
| SOE | E_rem(z; P_ref)/E_u(1; P_ref), with E_rem(z; P) = Q ∫_{z_c(P)}^{z} V_t dz | energy-based gauge |

* The initial condition from a displayed SOC s₀ is `z₀ = z_ref + s₀(1 − z_ref)`.
* **Testable consequence.** Under constant P, s(t) falls about 10 % faster per percent near empty than near full,
  because dz/dt ∝ 1/V. The mean V_t in the 99–50 / 50–20 / 20–1 % windows is in the ratio 1 : 0.92 : 0.90
  [PROTO]. SOE(t) falls linearly.
* The Bonn idle traces fall linearly to within 6 % (31.1, 30.8 and 32.4 min per % in the 99–50, 50–20 and
  20–1 % windows; vd [RUN]). V1.4 uses this to decide which gauge convention iOS follows.
* iOS also holds "100 %" for 6–7 % of the run. Model this as `s_disp = min(1, s/(1 − h))` with h = 0.065
  (vd b.4), for the Bonn comparison only.

### 1.7 Modelling assumptions

| # | Assumption | Justification | Tested by |
|---|---|---|---|
| A1 | Constant-power load at the terminals | PMIC regulation (dp §2.1); Bonn runs have near-constant mean power | M1 vs a constant-current variant (structural) |
| A2 | Component powers add | Bonn additivity residual ≤ 0.15 % (vd) | V1.3 hold-out |
| A3 | One thermal node; skin from a static R_skin/R_th ratio | three independent C_th estimates agree (ta T16–T22) | structural: isothermal vs M3 |
| A4 | Throttling is a smooth cap on SoC power | AOSP graded thermal levels (AOSP-TH); w and T_thr are ASSUMPTIONS | OAT on w, T_lim |
| A5 | Usage is Markov (exponential sojourns) | parsimony; analytic tractability | semi-Markov lognormal variant (§5.4.6) |
| A6 | Aging is frozen within one discharge | time-scale separation (days vs hours) | – |
| A7 | No OCV hysteresis | ≤ ~20 mV for LCO/graphite (ce a.2, unverified) | ±20 mV OCV offset in UQ |
| A8 | 4.45 V OCV = Ai2020 4.2 V curve + 0.242·z^3.86 | constructed (ce a.2, L) | OCV offset and V_cut in UQ |
| A9 | Pack R scales as R·Q = constant across phones | DFN R₀·Q ≈ 85–95 mΩ·Ah (ce a.3) | validation device configs |
| A10 | Background traffic is Poisson; radio tails overlap | closed-form tail fraction 1 − e^{−λT_tail} (dp §2.4) | explicit RRC-state variant |
| A11 | Battery-side powers include PMIC loss | measurements in dp §1.7 are at the battery terminals | – |
| A12 | Self-discharge ignored in TTE runs | 0.1 %/month vs about 0.3 %/h standby drain (ta S6) | M4 storage only |

---

## 2. Component power sub-models (battery-terminal, watts)

```
P_tot = s_on·[P_floor,on + P_disp] + (1 − s_on)·[P_sleep + P_bg(c, λ_bg)·m_sig]
        + P_SoC + P_net,act(a, c)·m_sig(RSRP) + P_GNSS + P_cam + P_aux
```

### 2.1 Floor, sleep, paging and background (Poisson tail formula)

* **Awake floor:** `P_floor,on = P_sleep + ΔP_awake = 0.015 + 0.080 = 0.095 W`. This also applies to screen-off
  activities that hold a wakelock (audio, voice).
* **Screen-off standby:**

  ```
  P_idle = P_sleep + P_bg,    P_bg = P_page,c + (P_tail,c − P_page,c)·(1 − e^{−λ_bg T_tail,c}) + λ_bg·E_cpu
  ```

  This is the overlapping-tail formula of dp §2.4, labelled ASSUMPTION A10. E_cpu = 0.1 J per wake-up (DERIVED:
  ΔP_awake × about 1 s). Coalescing N wake-ups into one burst replaces λ_bg by λ_bg/N in the tail term only
  (§7.5).

| Quantity | Value | Range | Source |
|---|---|---|---|
| P_sleep (SoC suspend, DRAM self-refresh, PMIC) | 15 mW | 8–30 | CH15 via slides: SoC suspend 3.8–5.1 mA; CH13: whole-phone suspend 24–27 mW |
| ΔP_awake (CPU awake, no display) | 80 mW | 50–250 | CH10 (Freerunner 269 vs 69 mW); AOSP example 50 mA (illustrative). **Lower-middle value informed by BONN26 idle, so a calibration (see V1)** |
| P_page: airplane / Wi-Fi / 4G / 5G NSA | 0 / 5 / 20 / 30 mW | Wi-Fi 2–10; 4G 5–35; 5G 20–60 | HU12 (Wi-Fi beacon 1.6 mW; LTE DRX 19.7 mW, DERIVED); CH15 (1.1 mA Wi-Fi; 2.3–8.3 mA paging); 5G = 1.5× 4G is an ASSUMPTION (GU24 reports about 2× 4G in screen-on standby) |
| P_tail, T_tail: 4G | 0.120 W, 10.2 s | 0.066–0.178 W; 5–11.6 s | NA21 Tables 2 and 7; HU12 |
| P_tail, T_tail: 5G NSA | 0.255 W, 10.3 s | 0.249–0.593 W; 10.2–21.4 s | NA21; XU20 (chained tails) |
| P_tail, T_tail: Wi-Fi | 0.12 W, 0.24 s | – | HU12 (negligible) |
| λ_bg | 1/60 s⁻¹ | 1/300–1/20 (log-uniform) | ASSUMPTION; the result is checked against the 0.15–0.5 %/h standby band (dp BG9) |

Resulting standby power [PROTO]: Wi-Fi 0.022 W, 4G 0.052 W, 5G 0.082 W. On 4G that is 0.29 %/h, or 14.2 days.

### 2.2 Display (OLED)

```
P_disp = P_d0 + p_A · A · g(b) · APL + κ_RR · max(0, f_r − 60 Hz)
g(b)   = 0.16 + 1.68 b          (b ≤ 0.5)
       = 1.0 + 2.0 (b − 0.5)    (b > 0.5)        through CH13's points g(min) = 0.16, g(0.5) = 1, g(max) = 2.0
```

**Parameters.**

* **P_d0 = 0.12 W** (0.08–0.25). CH13: 223 mW constant term on the 2012 S3, used as the upper value. A
  lower-middle default is consistent with the Pixel 7a display rail (GU24, about 0.30 W in total).
* **A = 97.4 cm²** for a 6.3-inch, 19.5:9 panel (geometry). It is 91.3 cm² at 6.1" and 110.2 cm² at 6.7".
* **p_A = 4.8 mW/cm²** at b = 0.5 and APL = 1 (range 3.0–7.1). DERIVED so that P_disp(b = 0.5, APL = 0.6) =
  0.40 W, the dp default for a 6.1–6.7" OLED (0.25–0.9; GU24: Pixel 7a about 0.30 W, S23+ about 0.45 W). The
  upper end, 7.1 mW/cm², is the 2012 S3 value: CH13 1124 mW at max white, less 223 mW, over 63.5 cm², divided
  by 2.
* **κ_RR = 2.8 mW/Hz** (1.0–3.0). AT20 via secondary reports: about 160–180 mW at 120 Hz vs 60 Hz on a
  fixed-rate panel (M/L).
* **APL by activity** (ASSUMPTIONS; ranges in §3):

  | Activity | APL |
  |---|---|
  | social, light UI | 0.50 |
  | web, light | 0.60 (DZ11: 65 % white pixels) |
  | dark-mode UI | 0.15 |
  | video | 0.25 |
  | game | 0.35 |
  | navigation map | 0.50 |
  | viewfinder and video call | 0.40 |

**Checks.**

* Max-white at max brightness: 1.06 W (dp gives 1.0–2.0, so we are at the low end).
* Dark-mode saving at b = 0.4 in social: 0.14 W, or 10 % of state power. That is +3.4 % TTE on the typical day,
  consistent with DH21's 3–9 % at 30–50 % brightness.

**Structural alternative S-g.** DH21 reports 39–47 % saving at 100 % brightness, which suggests a convex slider
law: `g(b) = 0.16 + 1.84 b²`. Run it as a structural variant.

### 2.3 SoC (CPU + GPU + memory), DVFS law, throttling

* **Demand per activity (W):**

  | Activity | P_SoC,dem (W) | Source |
  |---|---|---|
  | social | 0.45 | GU24: compute + memory 0.4–0.5 W in Instagram |
  | video, hardware decode | 0.20 | CH13: HW decode saves > 1 W vs software; HE20 |
  | game, 3D | **3.0** | CH13: 3D game 2.43 W total; NBC max-load totals 2.4–5.8 W |
  | navigation | 0.60 | dp: DERIVED 0.3–0.8 W map rendering |
  | video call | 0.60 | GU24: compute ≈ 0.50 W |
  | audio | 0.03 | ASSUMPTION |
  | camera | 0 (inside P_cam) | – |

* **DVFS law,** used for OS strategies and performance accounting (LX):

  ```
  P_dyn ∝ f·V(f)²,    V(f) = 0.47 V + 0.29 V/GHz · f
  ```

  This gives P ∝ f^1.9 along a real DVFS table (MT8183 A73). The performance ratio under a cap is therefore
  `ρ = (P_SoC/P_SoC,dem)^{1/1.9}`.
* **Throttling** (M3; §1.4):
  * T_lim = 45 °C (CHI15) and T_thr = T_lim − 2 K = 43 °C;
  * w = 1 K (ASSUMPTION, 0.5–2);
  * P_SoC,fl = 0.3 W (ASSUMPTION, 0.2–0.5);
  * the sharp-cap limit is `P_sust(T_amb) ≈ (T_thr − T_amb)/R_skin` (ta §2.2), i.e. 3.3 W at 25 °C and 1.45 W
    at 35 °C, within CHI15's 2.2–3.8 W.

### 2.4 Network radio

```
P_net = P_net,act(a, c) · m_sig(RSRP),      m_sig = 1 + γ_sig · max(0, (−85 − RSRP)/10 dB)
```

γ_sig = 0.26 per 10 dB (0.025–1.0) is DERIVED so that m_sig(−105) = 1.52. That reproduces DI13's +52 % energy
from −85 to −105 dBm. The upper range covers the Tx/Rx current ratios of 1.2–3.5×. m_sig applies to cellular
activity power and to the whole cellular background term P_bg (paging and tails). It does not apply to Wi-Fi.
Cell reselection and search under weak signal are not modelled, which is conservative.

**P_net,act at reference signal (W).**

| Activity | Wi-Fi | 4G | 5G NSA | Basis |
|---|---|---|---|---|
| social (feed, frequent small requests) | 0.10 (0.05–0.30) | 0.35 (0.20–0.95) | 0.50 (0.30–1.0) | GU24: Pixel 7a cellular rail ≈ 0.95 W in Instagram (Exynos modem, upper end); CH13: Wi-Fi ≈ 4× more efficient than cellular for web |
| video streaming | 0.10 (0.06–0.30) | 0.28 | 0.39 | Wi-Fi: Apple streamed − local video = 0.06–0.08 W (vd, DERIVED); 4G − Wi-Fi = +0.18 W (0.12–0.24) and 5G − 4G = +0.11 W (0.07–0.15), from GU24 Pixel 7a and S23+ YouTube totals |
| game | 0.03 | 0.12 | 0.25 | tail-level background traffic (NA21 tail powers) |
| navigation (tiles, traffic) | 0.05 | 0.20 (0.12–0.35) | 0.30 | ASSUMPTION between tail and streaming |
| video call (uplink-heavy) | 0.35 (0.2–0.8) | 1.50 (0.8–2.6) | 1.80 | GU24 Meet cellular rail 2.58 W (upper); HU12 LTE linear model ≈ 2.0 W at 1.5 Mbps up/down (2012 radio, upper); ASSUMPTION mid |
| camera (sync) | 0.02 | 0.05 | 0.08 | tail-level |
| audio streaming | 0.02 | 0.06 | 0.09 | tail-level |

Modern 4G tails (0.066–0.178 W) are 5–35× cheaper than the 2012 LTE constants (1.06 W for 11.6 s; dp finding
6). Using HU12 tails would overstate background drain, so HU12 is used only as an upper bound.

### 2.5 GNSS

```
P_GNSS = e_fix · ν_fix · D
```

* e_fix = 0.27 J per fix (0.23–0.32; KW19: 232 mJ single-frequency L1, 318 mJ dual-frequency L1+L5).
* ν_fix = 1 Hz in navigation, so P_GNSS = 0.27 W (CH13 S3: 0.39–0.43 W above idle, including the app).
* D is the duty factor of the receiver. D = 1 by default. D ∈ [0.1, 1] in the Bonn validation (§6), because
  iOS may duty-cycle a stationary receiver (unverified).
* A cold fix costs 154 J (CH13), which is 0.24 % of E_u: one-off and negligible.

### 2.6 Camera, audio, Bluetooth and negligible items

| Item | Value | Basis |
|---|---|---|
| P_cam, rear pipeline (sensor + ISP + CPU), above the awake floor | 1.15 W (0.9–2.0) | DERIVED from CX13 (Galaxy S2 recording 1.40–1.65 W total, AMOLED 0.17–0.18 W, minus an assumed 0.2 W floor); CH13 S3 2.26–2.61 W total gives the upper bound |
| Front camera in a video call | 0.30 W (0.2–0.5) | GU24 camera rail ≈ 0.26 W |
| Audio DSP/codec | 0.03 W | ASSUMPTION (CH13: max volume adds 5 mW) |
| Bluetooth audio link | 0.04 W (0.036–0.045) | CH10: +36–45 mW |
| Sensors (accelerometer 5 mW, gyroscope 30 mW at 50 Hz), flash storage < 8 mW, haptics | omitted from presets; included in the "surprisingly small" table | CH13; haptics unverified |
| Camera HDR/flash burst (cold scenario only) | burst total 4.0 W for 30 s | ASSUMPTION; inside CH13's > 7 W peak and NBC's 2.4–5.8 W max loads |

### 2.7 Activity presets and their powers [PROTO]

Conditions: b = 0.4 indoor, reference signal, 25 °C, before throttling.

| Activity a | s_on, b, APL, f_r | Components (4G) | P_a Wi-Fi / 4G / 5G (W) |
|---|---|---|---|
| idle (standby) | 0, –, –, – | sleep + paging + Poisson tails | 0.022 / 0.052 / 0.082 |
| social | 1, 0.4, 0.50, 120 | floor 0.095, display 0.483, SoC 0.45, net 0.35 | 1.128 / 1.378 / 1.528 |
| video (streaming) | 1, 0.4, 0.25, 60 | floor 0.095, display 0.217, SoC 0.20, net 0.28 | 0.612 / 0.792 / 0.902 |
| game (3D) | 1, 0.5, 0.35, 60 | floor 0.095, display 0.284, SoC 3.0, net 0.12 | 3.409 / 3.499 / 3.629 (demand); steady throttled ≈ 2.85 / 2.95 |
| nav | 1, 0.6, 0.50, 60 | floor 0.095, display 0.401, SoC 0.6, GNSS 0.27, net 0.20 | 1.416 / 1.566 / 1.666 |
| vcall | 1, 0.4, 0.40, 60 | floor 0.095, display 0.276, SoC 0.6, cam 0.30, net 1.50 | 1.621 / 2.771 / 3.071 |
| camera | 1, 0.6, 0.40, 60 | floor 0.095, display 0.344, cam 1.15, net 0.05 | 1.609 / 1.639 / 1.669 |
| audio (BT, screen off) | 0, –, –, – | floor 0.095, DSP 0.03, BT 0.04, net 0.06 | 0.185 / 0.225 / 0.255 |

---

## 3. Reference phone RP-63: full parameter table

RP-63 is a generic 2024-class phone with a 6.3-inch OLED (LTPO, 60/120 Hz), a 5G-capable modem and a
4.70 Ah / 4.45 V LCO-graphite pouch cell (18.2 Wh). Its capacity is that of current 6.3" flagships, for example
the Pixel 9 at 4.70 Ah typical and 4.558 Ah minimum (PIX9). It is not a model of any one product.

**Distributions for UQ.**

* LU: log-uniform over the range (used for positive quantities spanning at least 2×).
* U: uniform.
* "×(a–b)": a multiplicative factor, log-uniform.

| Group | Symbol | Value | Unit | UQ range (dist.) | Source / basis | Conf. |
|---|---|---|---|---|---|---|
| Cell | Q₀ | 4.70 | Ah | 4.56–4.80 (U) | PIX9 typical / minimum; ce a.1 | H |
| | U(z), 4.2 V base | −0.60117e^{−36.426z} + 3.6193 + 0.44311z − 0.48605z² + 0.63127z³ | V | δU ±0.020 V (U) | ce a.2: [RUN] fit to the Ai2020 OCV (CC BY), RMSE 7 mV | M |
| | high-voltage stretch | + 0.242·z^3.86 | V | (with δU) | ce a.2 construction (gives V_max 4.45 V, mean 3.87 V) | L |
| | V_max / mean OCV | 4.450 / 3.870 | V | – | ce a.1; [PROTO] mean 3.8699 | M |
| | V_cut (system cut-off) | 3.40 | V | 3.20–3.50 (U) | QC-FG, QC-SMB examples; LEE20 "≈ 3.4 V" | M/H |
| | R₀,ref (pack: cell + PCM + flex, z = 0.5, 25 °C) | 0.060 | Ω | 0.035–0.130 (LU) | ce synthesis: PINE pack ≤ 100 mΩ at 1 kHz, PCM 20–60 mΩ; BU802a 130 mΩ; [RUN] R₀·Q | L |
| | A_R, z₀ (R₀ rise at low SOC) | 0.80, 0.139 | – | A_R 0.5–2.0 (U) | ce a.3 [RUN] | M |
| | R₁,ref, τ₁,ref | 0.010 Ω, 30 s | Ω, s | R₁ 0.003–0.025 (LU); τ₁ 10–60 (LU) | ce a.3 [RUN] (Ai2020 35 s, Marquis 29 s) | M |
| | E_a,R0 / E_a,R1 / E_a,τ | 18 / 30 / 21 | kJ/mol | 15–25 / 25–55 / 15–30 (U) | ce a.3 [RUN], BARC22, SUN23; ta C6; vd HG2 apparent 26.4 | M |
| | ∂U/∂T(z) | ce Appendix A2 table (−0.25 mV/K at z = 0.5) | mV/K | ×(0.5–1.5) | ce [RUN] on Ai2020 | M–L |
| | KiBaM (M2-K): c, κ_K, E_a,K | 0.166, 40 h⁻¹, 35 kJ/mol | –, h⁻¹, kJ/mol | κ_K 33–47; E_a,K 25–55 | JH09; DERIVED from [RUN] DFN rate data; PyBaMM E_a | M/L |
| | η_charge | 0.99 | – | 0.98–1.00 | PLETT-N | M |
| Derived | E_n = Q₀·mean U | 18.19 | Wh | – | [PROTO] | – |
| | z_c0, E_u (P → 0) | 0.0263; 17.79 | –; Wh | – | [PROTO] | – |
| | z_ref (P_ref = 0.5 W) | 0.028 | – | – | [PROTO] | – |
| Thermal | C_th | 150 | J/K | 100–230 (U) | ta T16–T22 (three estimates) | M |
| | R_th (bulk) | 3.0 | K/W | 2–4 (U) | ta T14, T26 (HE19) | M/L |
| | R_skin | 5.5 | K/W | 4–9 (U) | ta T1 (PU23 Pixel 5: 5 K/W), T6 (KANG19) | M |
| | T_lim; T_thr = T_lim − 2 K | 45; 43 | °C | T_lim 41–48 (U) | CHI15; KANG19; IEC 62368-1 (UL-TS) | H |
| | w (throttle softness) | 1.0 | K | 0.5–2 (U) | ASSUMPTION (AOSP-TH graded levels) | L |
| | P_SoC,fl | 0.30 | W | 0.2–0.5 (U) | ASSUMPTION | L |
| Aging | k_cal | 1.44×10⁻³ | d^−½ | so q_cal(1 yr, 25 °C, 50 % SOC) = 2–4 % (U) | SCH14 via CA21 (2.4 % at 3.7 V); KEIL16 (2–3 % in 10 months) | M |
| | E_cal | 58 | kJ/mol | 20–60 (U) | SCH14 (6976 K); SCM18 LFP 20.6 | M |
| | g(z), SOC stress | piecewise-linear through (0, 0.67), (0.5, 1.0), (1, 1.83) | – | step variant (ta) | KEIL16 NCA at 25 °C (figure read-off) | M |
| | k_cyc | 4.9×10⁻³ | FCE^−½ | 3.5–7.2×10⁻³ (U) | DERIVED: q_cal + q_cyc = 0.20 at N₈₀ = 1000 FCE, 1 FCE/day, 25 °C (AP-BATT iPhone 15+; GOOG Pixel 8a+); 500-cycle rating gives 7.2×10⁻³ | M/L |
| | E_H, E_L, β_L, c_ref | 33 kJ/mol, 55 kJ/mol, 2.64 h, 0.5 h⁻¹ | – | E_H 25–45; E_L 40–70 | SCM18 (LFP; structure only); c_ref is an ASSUMPTION | M/L |
| | γ_R (r = γ_R q) | 2.5 | – | 1.5–6.6 (U) | SCH14 G8 (+41–65 % at 20 % fade); KEIL16 NCA storage | L |
| | self-discharge | 0.1 | %/month | 0.05–2 | ZIL19; BU802b (storage runs only) | M/L |
| Display | A | 97.4 | cm² | – | geometry, 6.3" 19.5:9 | – |
| | P_d0 | 0.12 | W | 0.08–0.25 (U) | CH13 (0.223 upper); GU24 | M |
| | p_A | 4.8 | mW/cm² | 3.0–7.1 (U) | DERIVED (§2.2) from GU24 and the dp default; CH13 upper | L/M |
| | g(b) | §2.2 | – | structural S-g | CH13 | H |
| | κ_RR | 2.8 | mW/Hz | 1.0–3.0 (U) | AT20 (secondary) | M/L |
| | APL multiplier | 1 | – | ×(0.7–1.3) | ASSUMPTION | L |
| | b indoor / outdoor | 0.4 / 1.0 | – | ±0.1 (U) | DH21 (30–40 % auto indoor) | M |
| SoC | P_SoC,dem per activity | §2.3 | W | ×(0.7–1.4) global; game 2.0–4.5 (U) | GU24, CH13, NBC, dp | M |
| Base | P_sleep; ΔP_awake | 0.015; 0.080 | W | 0.008–0.030; 0.05–0.25 (LU) | §2.1 | M/L |
| | E_cpu per wake-up | 0.10 | J | 0.05–0.30 (LU) | DERIVED | L |
| Network | P_page, P_tail, T_tail | §2.1 | – | §2.1 ranges (U) | HU12, NA21, XU20, CH15 | H/M |
| | P_net,act table | §2.4 | W | global ×(0.5–2.0) (LU) plus the per-cell ranges | GU24, NA21, Apple, HU12 | M/L |
| | γ_sig | 0.26 | per 10 dB | 0.025–1.0 (LU) | DI13 (DERIVED) | M |
| | λ_bg | 1/60 | s⁻¹ | 1/300–1/20 (LU) | ASSUMPTION | L |
| GNSS | e_fix; ν_fix | 0.27 J; 1 Hz | – | 0.23–0.32 (U) | KW19 | H |
| Camera | P_cam rear / front | 1.15 / 0.30 | W | 0.9–2.0 / 0.2–0.5 (U) | CX13, CH13, GU24 | M |
| Aux | BT; audio DSP | 0.04; 0.03 | W | 0.036–0.045; 0.01–0.06 | CH10; ASSUMPTION | H/L |
| Usage | CTMC rates | §4.1 | – | λ_on ×(0.7–1.3); d_j ×(0.7–1.3) each (LU); p ~ Dirichlet(50 p) | ASSUMPTION | L |

**Validation device configurations** (same physics, different constants).

* iPhone 14 Pro: E = 12.38 Wh; iPhone 13 Pro: 11.97 Wh (N9M22 via vd).
* iPhone 14 Pro Max: 16.68 Wh. Pixel 8: 4.575 Ah × 3.87 V. Pixel 8a: 4.492 Ah. S24: 4.0 Ah (vd).
* Q = E/3.87 V. R scaled by R·Q = constant (A9). Display area set by the diagonal.

---

## 4. Usage scenarios

Common settings: start from s₀ = 100 % (z₀ = 1). T(0) = T_amb (the phone is soaked at ambient). New battery
unless stated. RAT, b and RSRP as listed. The TTE column gives the [PROTO] M3 result, with M0 and M1 as noted.

| ID | Scenario | u(t) definition | TTE [PROTO] |
|---|---|---|---|
| S1 | Standby, overnight-type | constant: a = idle, c = 4G, λ_bg = 1/60 s⁻¹, 25 °C | 340 h (14.2 d; M0 340.1, M1 340.0) |
| S2 | Light/social, Wi-Fi | constant: a = social, c = Wi-Fi, b = 0.4, f_r = 120 | 15.6 h (M0 15.78) |
| S3a | Video streaming, Wi-Fi | constant: a = video, c = Wi-Fi | 28.9 h (M0 29.05) |
| S3b | Video streaming, 5G NSA | constant: a = video, c = 5G | 19.6 h (M0 19.72): **−32 %** vs S3a |
| S4 | Navigation (GNSS 1 Hz, 4G, car mount) | constant: a = nav, c = 4G, b = 0.6, 25 °C. S4-out: b = 1.0 (sun) | 11.2 h; S4-out 10.0 h |
| S5 | 3D gaming, sustained | constant: a = game, c = 4G, b = 0.5, P_SoC,dem = 3.0 W. Variants at 15 / 25 / 35 °C and without throttling | 25 °C: 5.90 h (steady P = 2.84 W); 15 °C: 4.93 h; **35 °C: 11.09 h** (P = 1.53 W, ρ = 0.57); no throttle: 4.96 h |
| S6 | Video call over 4G | constant: a = vcall, c = 4G. S6-w: over Wi-Fi | 6.29 h; S6-w ≈ 11.0 h (M0) |
| S7 | Cold outdoor with bursts ("winter hike") | periodic, 600 s cycle: 570 s at 0.60 W (screen-off tracking: floor + GNSS 1 Hz + 4G tail + …), then 30 s at 4.0 W (camera + HDR + flash + upload). S7a T_amb = −10 °C; S7b −20 °C; S7c 25 °C control | S7a: 21.3 h, s_off = 6.8 %. S7b: 19.3 h, s_off = 15.5 %. Aged (S8c below): −10 °C 14.7 h, s_off 20 %; −20 °C 10.7 h, **s_off 42 %** |
| S7′ | Cold steady navigation | constant: nav, b = 1.0, T_amb = −10 °C, phone cold-soaked | 9.81 h vs 10.01 h at 25 °C (−2 %; self-heating keeps T ≈ −4 °C) |
| S8a | Aged battery, typical day | S10 profile, SOH = 0.80, r = γ_R·0.2 = 0.5, 25 °C | ≈ 32.2 h (quasi-static, −20 %) |
| S8b | Aged and cold, typical day | as S8a at T_amb = 0 °C | ≈ 32.0 h (quasi-static) |
| S8c | Aged and cold, bursty | S7 profile, SOH = 0.80 | above |
| S9 | Weak signal | constant: a = social, c = 4G, RSRP = −110 dBm (m_sig = 1.65) vs −85 dBm | 10.94 h vs 12.77 h (**−14 %**) |
| S10 | Typical day (stochastic) | product CTMC, §4.1 "typical", connectivity Wi-Fi ⇄ 4G (mean 2.5 h each), 25 °C | M3 Monte Carlo: mean 39.9 h, SD 6.6 h, 5/50/95 % = 29.1 / 39.8 / 51.0 h |
| S11 | Heavy day (stochastic) | §4.1 "heavy" | M3 Monte Carlo: mean 21.6 h, 5/50/95 % = 14.5 / 21.4 / 29.7 h |
| S12 | Hot climate, typical day and aging | S10 at T_amb = 32 °C plus the M4 3-year run (§7.3) | SOH after 2 years: 0.767 vs 0.820 at 25 °C (charge to 100 %) |
| S13 | Initial-SOC sweep | S2, S3a, S5 and S10 from s₀ ∈ {100, 80, 50, 20, 10} % | figure F4 |

S8b uses quasi-static M1 at T = T_amb + 3 K/W × P̄. Implementation runs full M3 Monte Carlo.

### 4.1 CTMC usage profiles

**Activity chain (hub and spoke).** From idle, sessions start at rate λ_on. A session of activity j is chosen
with probability p_j, lasts Exp(mean d_j), and then returns to idle. The generator Λ_A has entries
Λ_A[idle, j] = λ_on·p_j and Λ_A[j, idle] = 1/d_j.

**Connectivity chain.** Wi-Fi ⇄ 4G with rates 1/2.5 h⁻¹ each way. The product generator is the Kronecker sum
`Λ = Λ_A ⊗ I + I ⊗ Λ_C` (16 states). A 5G variant swaps 4G for 5G.

| Profile | λ_on (h⁻¹) | (p_j; d_j in min) for social, video, game, nav, vcall, camera, audio |
|---|---|---|
| typical | 4 | (0.55; 5), (0.12; 15), (0.05; 20), (0.04; 20), (0.04; 10), (0.08; 2), (0.12; 30) |
| heavy | 4 | (0.45; 8), (0.20; 25), (0.15; 30), (0.05; 25), (0.05; 15), (0.05; 3), (0.05; 30) |

These rates are ASSUMPTIONS, chosen for a plausible screen-on fraction: 27 % (6.5 h/day) typical and 48 %
heavy. They are checked against Google's median-user claim for the Pixel 8a (V2). Uncertainty is carried in §7.1.

**Typical-day results [PROTO].**

* Stationary activity distribution π: idle 0.588, social 0.108, video 0.071, game 0.039, nav 0.031, vcall 0.016,
  camera 0.006, audio 0.141.
* Energy shares: social 0.307, game 0.258, video 0.113, nav 0.106, vcall 0.078, audio 0.066, idle 0.050,
  camera 0.023.
* Component shares, which equal −(TTE elasticity) in M0: SoC 0.434, display 0.221, network 0.147, awake floor
  0.089, standby 0.050, camera 0.027, GNSS 0.019, BT 0.013.

**Time-inhomogeneous variant (Monte Carlo only, structural check).** Use λ_on(t) = 0 from 00:00 to 07:00
(night standby) and 5.5 h⁻¹ otherwise, which keeps about the same daily screen time.

---

## 5. Time to empty: definition, numerics and analytic results

### 5.1 Definition

```
TTE(s₀, u, θ) = inf{ t > 0 : V(t) ≤ V_cut  or  U_e(t)² − 4 R P(t) ≤ 0 }       (M1–M4)
TTE_M0        = inf{ t > 0 : E(t) ≤ 0 }
```

* TTE is measured from the moment the stated initial SOC applies.
* Report alongside it the displayed SOC at shutdown, s_off = (z(TTE) − z_ref)/(1 − z_ref). s_off > 2 % means
  an "unexpected shutoff".
* Also report the time to displayed 0 %, t(s = 0). At 25 °C on a new battery it coincides with TTE to within
  1.5 % for P ≤ 3.5 W: s(TTE) = −0.2 % at 0.05 W and +1.4 % at 3.4 W [PROTO].
* For stochastic usage, TTE is a random variable. Report its mean, SD, 5/50/95 % quantiles and P(TTE < 24 h).

### 5.2 Numerical method

**Stiffness.**

* Time scales: τ₁ = 30 s at 25 °C, rising to 135 s at −20 °C (2RC variant: τ = 6.7 s).
* τ_th = 450 s; SOC: 5–340 h; aging: days.
* The stiffness ratio is about 10³–10⁵. It is mild for explicit RK4 if h ≤ τ_min/3.
* The current is algebraic but explicit (quadratic root), so there is no differential-algebraic system.

**Integrator.** Use the existing vectorised `rk4_batch` (`src/phonebatt/integrate.py`). Step sizes:

| Case | h |
|---|---|
| bursty or throttled runs (S5, S6, S7) | 2 s |
| CTMC runs | 5 s |
| constant light loads (S1–S4) | 10 s |

Optional splitting update for v₁, exact for piecewise-constant I over a step:
`v₁ ← v₁ e^{−h/τ₁} + R₁ I (1 − e^{−h/τ₁})`. It makes h up to 30 s safe in standby.

**Events.** Use two event functions, g₁ = V − V_cut and g₂ = U_e² − 4RP. The crossing time is found by linear
interpolation inside the step (already implemented: second-order accurate, NaN treated as an event).

**Inputs.**

* CTMC paths are pre-sampled as jump-time lists and mapped to the step grid; a jump falls on the nearest step,
  with error at most h per sojourn and no bias.
* Deterministic periodic inputs are evaluated analytically at each stage.

**Reference solver.** For deterministic scenarios use `scipy.integrate.solve_ivp(method="LSODA", rtol=1e-8,
atol=1e-10)` with terminal events, segment by segment at input discontinuities.

**Acceptance.** |TTE(h) − TTE(h/2)| / TTE < 0.1 %, and the RK4 vs LSODA difference < 0.2 %, on S2, S5 and S7a.

**Long horizons.** For S1 standby and M4, use the M1 quadrature (§5.3.2) or h = 30 s with the exponential v₁
update. Aging uses daily operator splitting: integrate one representative day of the fast model, then advance
(q_cal, q_cyc) one day with the day's T and SOC histograms.

### 5.3 Analytic and approximate TTE formulas

**5.3.1 M0.**

* Constant P: `TTE = SOE₀·E_u/P`.
* Piecewise-constant deterministic u: TTE solves `Σ_k P_k Δt_k = SOE₀ E_u`.
* **Elasticity and drain decomposition.** With `P̄ = Σ_i P̄_i` (time-averaged component powers),

  ```
  ∂ln TTE/∂ln P_i = −P̄_i/P̄ = −s_i
  TTE_without_i / TTE = 1/(1 − s_i)
  ```

  So the energy-share pie chart *is* the elasticity chart.
* A small change δP gives `ΔTTE/TTE = −δ + δ² − …` with δ = δP/P̄.

**5.3.2 M1, exact quadrature (constant P, quasi-steady R_s).** Use dt = 3600·Q·dz/I and
1/I = (U + √(U² − 4R_sP))/(2P):

```
TTE_M1(z₀; P) = (3600 Q_eff / P) ∫_{z_c(P)}^{z₀} V_t(z; P) dz,     V_t = ½ [U + sqrt(U² − 4 R_s(z) P)]
z_c(P) solves  U(z_c) = V_cut + P·R_s(z_c)/V_cut                   (scalar root; brentq on [1e−5, z₀])
```

**First-order perturbation in P.** Use V_t = U − P R_s/U + O(P²) and z_c(P) − z_c0 ≈ P R_s(z_c0)/(V_cut U′(z_c0)):

```
E_u(P) ≈ E_u0 (1 − a P),     a = [ ∫_{z_c0}^{1} R_s/U dz + R_s(z_c0)/U′(z_c0) ] / ∫_{z_c0}^{1} U dz
TTE_M1 ≈ (E_u0/P)(1 − a P),   ∂ln TTE/∂ln P = −1 − aP/(1 − aP)
```

* The two terms are **ohmic loss** (about 0.0050 W⁻¹) and **cut-off stranding** (about 0.0032 W⁻¹; U′(z_c0)
  ≈ 8.8 V per unit SOC).
* RP-63 at 25 °C: **a ≈ 0.0082 W⁻¹**, matching the [PROTO] quadrature (M1/M0 − 1 = −0.95 % at 1.13 W and
  −3.05 % at 3.4 W).
* In the cold or with age, a scales roughly with the Arrhenius factor and with (1 + r). TTE then becomes
  noticeably superlinear in P (§5.3.4).

**5.3.3 M2 bursts and recovery.**

* For a burst of length t_b, the effective resistance is R_eff(t_b) = R₀ + R₁(1 − e^{−t_b/τ₁}).
* For bursts ≪ τ₁ only R₀ matters. The cold-shutdown map (§5.3.4) is evaluated with R_eff(t_b).
* After a shutdown at current I, the open-circuit relaxation lifts V by v₁ ≈ R₁ I within about 3τ₁.
* **M2-K.** The stranded fraction at constant C-rate C is `s_K = C/κ_K`: 0.13 % at 0.05C and 0.50 % at 0.2C
  [PROTO]. Phone loads are 0.01–0.2C, so the rate-capacity effect is negligible. This is consistent with the
  Bonn runs and with DFN (0.2–1 % at 0.2–0.5C, ce).

**5.3.4 Cold, aged and low-SOC shutdown (closed form).** Quasi-steady, for bursts ≥ 3τ₁:

```
P_max,cut(z, T, r) = V_cut [U(z,T) − V_cut] / R_s(z, T, r)
z_off(P_pk, T, r) solves  U(z_off, T) = V_cut + P_pk R_s(z_off, T, r)/V_cut
s_off = (z_off − z_ref)/(1 − z_ref)
```

P_max,cut values [PROTO]:

| Condition | z = 0.5 | z = 0.2 | z = 0.1 |
|---|---|---|---|
| 25 °C, new | 20.2 W | 12.5 W | 9.1 W |
| −10 °C, new | 6.8 W | 4.3 W | 3.1 W |
| −10 °C, SOH 0.8 | 4.5 W | 2.9 W | 2.1 W |
| −20 °C, SOH 0.8 | 3.1 W | 2.0 W | 1.4 W |

A 3–4 W burst is therefore infeasible on a cold, aged phone at mid-to-low SOC. This is the mechanism of the
unexpected shutdowns described by Apple (AP-BATT; LEE20 measured 15–35 % displayed SOC at −5 °C on 2015–17
phones).

The fraction of charge stranded relative to warm conditions is `Δs_off ≈ z_off(P_pk, T, r) − z_c0`. In the
cold-soaked limit (no self-heating), the relative TTE loss is about z_off/(1 − z_c0) for z_off ≫ z_c0.

**5.3.5 M3 thermal.**

* Step response: `T(t) = T_amb + R_th P_heat + (T₀ − T_amb − R_th P_heat) e^{−t/τ_th}`.
* Steady-state skin temperature: T_skin = T_amb + R_skin·P_heat.
* **Throttled gaming**, sharp-cap limit w → 0:

  ```
  TTE_game ≈ E_u(P*)/P*,     P* = min(P_dem, (T_thr − T_amb)/R_skin − P_batt,heat)
  dTTE/dT_amb = E_u/(R_skin P*²) > 0    whenever the cap binds
  ```

  The performance ratio is `ρ = ((P* − P_non-SoC)/P_SoC,dem)^{1/1.9}`.
  [PROTO], smooth cap: 25 °C P* = 2.84 W, TTE 5.90 h, ρ ≈ 0.88; 35 °C P* = 1.53 W, 11.1 h, ρ ≈ 0.57.
* **Cold self-heating.** The battery runs at T ≈ T_amb + R_th P. At P = 1.75 W and −10 °C that is −4.8 °C, so the
  Arrhenius factor on R₀ falls from 2.6 (cold-soaked) to 2.2.

**5.3.6 M4 aging.** At 25 °C, `TTE(q) ≈ TTE_new·(1 − q)·(1 − a P γ_R q)`. The second factor is small at 25 °C
(< 1 % at P = 1.5 W, q = 0.2) and large in the cold. [PROTO] quasi-static TTE:

| P | 25 °C, SOH 1.0 / 0.9 / 0.8 / 0.7 | −10 °C, SOH 1.0 / 0.9 / 0.8 / 0.7 |
|---|---|---|
| 0.44 W | 40.3 / 36.2 / 32.2 / 28.1 h | 40.0 / 35.9 / 31.8 / 27.8 h |
| 3.0 W | 5.81 / 5.20 / 4.59 / 3.99 h | 5.54 / 4.85 / 4.15 / 3.43 h |

At 3.0 W and −10 °C the cut-off rises from z = 0.063 to z = 0.157 between SOH 1.0 and 0.7.

### 5.4 Stochastic usage: PDMP results

Let X(t) be the CTMC on 𝒮 with generator Λ, stationary distribution π, and state powers P_x. In M3 the powers
are the quasi-static (throttled) values. The battery state follows `dE/dt = −P_{X(t)}` (M0) or the M1–M3 flow
with P = P_{X(t)}: a PDMP.

**5.4.1 Fluid mean.** `P̄ = πᵀP` and `TTE_fluid = E_u/P̄`.

**5.4.2 Exact mean (Dynkin / Poisson-equation identity).**

* Let D = (1πᵀ − Λ)⁻¹ − 1πᵀ be the deviation matrix and h = D P (so −Λh = P − P̄1 and πᵀh = 0).
* Dynkin's formula applied to h(X_t) up to the stopping time gives, exactly for M0:

  ```
  E[TTE | X₀ = x₀] = ( E_u − h(x₀) + E[h(X_TTE)] ) / P̄
  E[h(X_TTE)] ≈ Σ_x π_x P_x h_x / P̄          (depletion happens preferentially in high-power states)
  ```

* The second line is the inspection-paradox correction.
* [PROTO] typical day: fluid 40.39 h; corrected 41.19 h; Monte Carlo (20 000 paths) 41.18 h.
  Heavy day: 21.98 / 22.79 / 22.62 h.

**5.4.3 CLT variance.**

* The asymptotic variance rate of the consumed energy is

  ```
  σ² = 2 Σ_x π_x (P_x − P̄)(D(P − P̄1))_x           [W² h]
  ```

* Inverting the drift-diffusion gives `Var(TTE) ≈ σ² E_u/P̄³`, so

  ```
  CV_TTE = sqrt( 2 τ_c / μ ) · CV_P,      CV_P² = Var_π(P)/P̄²,   τ_c = σ²/(2 Var_π P)  (integrated autocorrelation time)
  ```

* **Interpretation:** relative TTE uncertainty shrinks like 1/√(TTE/τ_c). Long, power-contrasting sessions (large
  τ_c and CV_P) widen it.
* [PROTO] typical: σ² = 0.2265 W²h, CLT SD 6.87 h vs Monte Carlo 6.56 h (CV 0.17). Heavy: SD 4.83 vs 4.68 h
  (CV 0.22). The CLT is about 5 % conservative.

**5.4.4 Hub-and-spoke closed form.** For the activity chain alone (renewal-reward over idle → session cycles of
mean length E[L] = 1/λ_on + Σ_j p_j d_j):

```
P̄  = (P_idle/λ_on + Σ_j p_j d_j P_j) / E[L]
σ² = 2 Σ_j p_j d_j² (P_j − P̄)² / E[L]                  (the idle–cross terms cancel exactly)
```

* Each activity's variance contribution is `2 p_j d_j²(P_j − P̄)²/(σ² E[L])`.
* [PROTO] typical: game 0.72, nav 0.11, vcall 0.07, social 0.05. Heavy: game 0.91.
* So **a few 20–30 min gaming sessions dominate TTE uncertainty** while taking 4 % of the time.
* The closed form matches the deviation-matrix σ² to machine precision [PROTO, typical, 4G only:
  0.238185 W²h both ways, difference 8×10⁻¹⁶].

**5.4.5 Inverse-Gaussian law.** Treat consumption as Brownian motion with drift P̄ and variance σ². Then
`TTE ~ IG(mean μ_D, shape λ_IG = E_u²/σ²)`, with μ_D the Dynkin mean.
[PROTO] typical 5/50/95 %: IG 30.7 / 40.6 / 53.8 h vs Monte Carlo 30.5 / 41.2 / 52.0 h.

**5.4.6 Semi-Markov generalisation.** For non-exponential sojourns with squared coefficient of variation CV_j²
(CV₀² for idle):

```
σ² = [ Σ_j p_j (1 + CV_j²) d_j² (P_j − P̄)² + (P_idle − P̄)²(CV₀² − 1)/λ_on² ] / E[L]
```

It reduces to §5.4.4 at CV = 1. Heavy-tailed sessions (lognormal with CV_j = 1.5) multiply the variance by
(1 + 2.25)/2 = 1.6.

**5.4.7 Convexity (Jensen) penalty in M1.** With the internal loss a·P² per unit time,

```
E[TTE]_M1 ≈ E[TTE]_M0 · (1 − a P̄ − a Var_π(P)/P̄)
```

* The last term, −a P̄ CV_P², is the cost of burstiness at equal mean power. For the typical day it is about −1 %.
* The cut-off is set by the power of the state the phone is in as z falls through z_c(P_x), so bursts near empty
  raise the effective cut-off: [PROTO] z_off is 0.033 on average and at most 0.041, vs z_c0 = 0.026.
* [PROTO] M3 Monte Carlo / M0 Monte Carlo mean ratio: **0.967 (typical), 0.954 (heavy)**.

**5.4.8 Usage-parameter elasticities (hub chain, M0).**

```
∂ln E[TTE]/∂ln P_x = −π_x P_x/P̄
∂ln TTE/∂ln d_j   = −π_j (P_j − P̄)/P̄
∂ln TTE/∂ln λ_on  = −π_idle (1 − π_idle)(P̄_session − P_idle)/P̄
```

[PROTO] typical: d_game −0.22, d_social −0.20, d_audio +0.075 (longer music sessions *extend* TTE because they
displace higher-power use), λ_on −0.54.

---

## 6. Validation plan

Principle: parameters come from the literature. The data below are **not** used for fitting, except where marked
CAL. Report the relative error ε = (model − observed)/observed and whether the observation lies in the model's
90 % parameter-UQ interval (coverage).

**V1. Bonn iPhone dataset** (BONN26, CC BY 4.0; summaries in `docs/research/data/validation/bonn_iphone_*`).

* Configuration: E = 12.38 Wh (14 Pro) or 11.97 Wh (13 Pro). OLED on and black, so the display term is P_d0 only.
  Airplane mode (P_page = P_net = 0). 25 °C. Floor = P_sleep + ΔP_awake.
* V1.1 **CAL**: black-screen idle. The floor default (0.215 W in total) was chosen within the literature ranges to
  be consistent with this run. Report it as a calibration check, not a validation.
* V1.2 **Hold-out**: rear camera, GPS "best", 13 Pro idle. Targets are the medians and min–max ranges.
* V1.3 **Additivity hold-out**: camera + GPS, predicted from the components.
* V1.4 SOC shape: the slope ratio over the 99–50 / 50–20 / 20–1 % windows for s(t) (charge gauge) vs SOE(t)
  (energy gauge), against 31.1 / 30.8 / 32.4 min/%.

| Target | Observed (median) | Model a priori [PROTO, M1] | ε | Comment |
|---|---|---|---|---|
| 14 Pro idle (CAL) | 55.0 h | 56.2 h | +2.1 % | calibration-informed |
| 14 Pro camera | 8.43 h | 8.73 h | +3.6 % | P_cam prior 1.15 W |
| 14 Pro GPS best | 48.5 h | 24.8 h (D = 1) – 49.9 h (D = 0.1) | −49 % / +3 % | **fails at D = 1.** The +0.03 W observed implies receiver duty cycling; reported as a "surprisingly little" finding and carried as D ∈ [0.1, 1] |
| 14 Pro camera + GPS (D = 1) | 8.38 h | 7.27 h | −13 % | additivity itself holds in the data (residual 0.15 %) |
| 13 Pro idle | 50.8 h | 54.3 h | +6.9 % | pre-used units (−3 % capacity, SCH26) |
| 13 Pro camera | 9.53 h | 8.44 h | −11 % | – |
| 13 Pro GPS best | 29.8 h | 24.0 h (D = 1) | −19 % | – |

Pass criteria: |ε| ≤ 15 % on hold-outs, or the observation inside the 90 % UQ interval.

**V2. Published specifications band** (cited facts, not open data). Tolerance ±25 %. Brightness for "150 nits"
tests: b = 150/L_max with L_max = 1000 cd/m² (14 Pro typical max per Apple; the same value is an ASSUMPTION for
Android phones). Apple tests use b = 0.4 (brightness undisclosed).

| Benchmark (source) | Observed | Model [PROTO] | ε |
|---|---|---|---|
| iPhone 14 Pro local video (APL14P) | 23 h | 23.8 h | +3 % |
| iPhone 14 Pro streamed video (APL14P) | 20 h | 19.8 h | −1 % |
| iPhone 14 Pro audio, BT (APL14P) | 75 h | 65.3 h | −13 % |
| iPhone 14 Pro Max Wi-Fi web, 150 nits (NBC22) | 15.8 h | 15.2 h | −4 % |
| iPhone 14 Pro Max local video, 150 nits (NBC22) | 28.1 h | 34.6 h | +23 % |
| iPhone 14 Pro Max idle at minimum brightness (NBC22) | 57.9 h | 61.2 h | +6 % |
| iPhone 14 Pro cellular web, 150 nits (LMAG22) | 10.4 h | 9.2 h | −11 % |
| Pixel 8 Wi-Fi web (NBC23) | 13.6 h | 16.4 h | +21 % |
| Galaxy S24 Wi-Fi web (NBC24) | 14.7 h | 14.4 h | −2 % |
| Pixel 8a typical mixed use (GPX) | ≥ 31 h | 38.5 h (S10 profile) | +24 % |
| Max-load power (NBC22–24) | 2.4–5.8 W | game unthrottled at b = 1: 3.7 W | inside |
| Standby drain band (dp BG9, DERIVED) | 0.15–0.5 %/h | 0.29 %/h | inside |

**V3. Cell-level temperature: LG HG2** (KOL20, CC BY 4.0; 18 CSVs in `docs/research/data/validation/lg_hg2_raw`).

* CAL at 25 °C only: the OCV from the C/20 curve (table), and R₀ + R₁ = 0.038 Ω split 6 : 1.
* Predict, using our activation energies and no refit, the 0.5C and 2C capacity ratios at −20, −10, 0, 10 and
  40 °C. Observed at 0.5C: 0.73 / 0.86 / 0.91 / 0.93 / 0.97. At 2C: 0.49 / 0.73 / 0.88 / 0.93 / 0.98.
* Expected outcome [PROTO estimate]: M3 alone gives R_s(−20 °C) ≈ 0.165 Ω vs the apparent 0.25 Ω, so it will
  **under-predict** cold loss. The structural variant M3 + M2-K (κ_K Arrhenius, E_a,K = 35 kJ/mol) is the test of
  whether diffusion-limited availability closes the gap.
* Targets: |Δf_Q| ≤ 0.05 for T ≥ 0 °C; ±0.10 at −20 °C; RMSE_V ≤ 50 mV for T ≥ 0 °C.
* Also report the apparent E_a of R₀ + R₁ against the 26.4 kJ/mol fitted to HG2.

**V4. Aging: NASA PCoE B0005** (no licence stated, public; vd).

* Inputs: the measured fade 1.862 → 1.331 Ah, and the measured R_e/R_ct growth (EIS).
* Check that M4's two effects explain the fresh → aged TTE at 2 A, 55.8 → 39.7 min (−28.9 %): (i) Q fade alone
  gives −28.5 %; (ii) the remainder must come from R growth through the cut-off.
* Also run the 4 °C (B0055) vs 24 °C capacity ratio, 0.71, as a harsh temperature check.

**V5. RC dynamics: Oxford Artemis 1 Hz drive cycle** (ODbL; optional). Drive the 1RC model with the logged
current at 40 °C. Take the OCV from the low-current segments (CAL) and report the V(t) RMSE. MANUAL_REQUIRED if
the pseudo-OCV file (254 MB) cannot be downloaded.

**V6. Physics reference: PyBaMM DFN, Ai2020** (BSD-3 code, CC BY data).

* Run constant-power phone loads (0.5, 1.5 and 3 W scaled to 2.28 Ah, i.e. the same C-rates) at 25, 0 and −10 °C,
  plus a 30 s / 4C-equivalent burst train.
* Compare with M1, M2, M2-K and M3 (isothermal) TTE and V(t).
* Target: TTE difference ≤ 3 % at 25 °C. Report the difference in the cold as structural error.

**V7. Published phone-level cold result** (LEE20; published measurement, cited).

* Build a Nexus-5X-like configuration: Q = 2.7 Ah, R_int(20 °C) = 0.25 Ω total, V_cut = 3.4 V, with a bursty UI
  load (I_max 1–2 A).
* Check the reported −15 % TTE at 0 °C (159 → 135 min) and shutdown at displayed SOC up to 33 % at −5 °C.
  Pass if inside the UQ band.

---

## 7. Analysis plans

### 7.1 Uncertainty quantification

**Sources of uncertainty.**

1. **Parameters (epistemic):** the §3 distributions, d ≈ 30.
2. **Usage (aleatory):** CTMC randomness at fixed θ.
3. **Usage-model parameters:** λ_on, p and d.
4. **Initial-SOC reading:** ±2 % displayed (ASSUMPTION) and the s ↔ z map.
5. **Model form:** the structural ensemble (§7.2.4).

**Propagation.**

* **Deterministic scenarios (S1–S9):** Monte Carlo over θ. N = 4000 Latin-hypercube samples on M3 via
  `rk4_batch`. Report median, 5–95 % interval and CV.
* **First-order check (M0 delta method):** `CV_TTE² ≈ CV_E² + Σ_i s_i² CV_{P_i}²`. Compare with the Monte Carlo
  result.
* **Stochastic scenarios (S10, S11):** nested sampling, 400 θ × 200 usage paths with M3. Apply the law of total
  variance:

  ```
  Var(TTE) = E_θ[Var(TTE | θ)] + Var_θ(E[TTE | θ])
  ```

  The first term is usage and the second is parameters. Report the fraction from usage. Use the §5.4 closed
  forms for the inner moments as a cross-check (inner Monte Carlo with 200 paths vs analytic: target ≤ 3 % on
  mean, ≤ 10 % on SD).
* **Bonn-constrained posterior (optional, ABC rejection):** keep θ samples whose predicted Bonn 14 Pro idle TTE
  is within ±5 % of 55.0 h. Report how much the reference-phone intervals narrow. This makes the single
  calibration explicit and quantitative.

### 7.2 Sensitivity analysis

**7.2.1 Local elasticities.**

* Compute `e_θ = ∂ln TTE/∂ln θ` by central differences (±1 %) for every θ in every scenario, using M3.
* Check against the analytic values: −s_i in M0 (§5.3.1), −(1 + aP) for P in M1 (§5.3.2), and the CTMC
  elasticities (§5.4.8).
* Present as a heat map (F8).

**7.2.2 OAT tornado.** Move each θ to its range endpoints, with the others nominal, for S2, S5, S7a and S10
(mean). Rank by |ΔTTE| (F9).

**7.2.3 Global Sobol.**

* Use SALib Saltelli sampling for first-order S₁ and total S_T indices, with 95 % bootstrap CIs (1000 resamples).
* Constant scenarios: N = 1024 base samples × (2d + 2) on the **M3 quasi-static surrogate**. The surrogate is the
  M1 quadrature at the steady thermal state T = T_amb + R_th P_heat (with throttled P* for S5). It is accurate
  when τ_th, τ₁ ≪ TTE, so check it against full M3 on 200 points (≤ 1 % error required).
* Dynamic scenarios (S7a, S10 mean): N = 256 on full M3.
* Expected leaders:
  * S2: the network and SoC multipliers, p_A, P_d0, b;
  * S5: R_skin, T_lim and the SoC demand (because of throttling);
  * S7: R₀,ref, E_a,R0, V_cut, A_R and γ_R;
  * S10: usage rates (λ_on, d_game) beside the power multipliers.

**7.2.4 Structural sensitivity.** Run every scenario under each variant and report ΔTTE (%) against full M3 (F5).

| Variant | Mechanism changed | Expected size [PROTO where available] |
|---|---|---|
| M0 | energy bucket | +0.5 to +3 % (25 °C); +3.4 % (S10 mean); larger in the cold and aged |
| M1 | no RC, no thermal | ≤ 0.2 % for constant loads; matters for S7 bursts |
| M2 (no thermal) vs M3 | isothermal at T_amb | gaming −16 % (no throttle); cold: more pessimistic (no self-heating) |
| M2-K vs M2 | bound charge instead of RC | < 0.5 % at 25 °C; tested on HG2 in the cold |
| 2RC (2.5 mΩ / 6.7 s + 14.1 mΩ / 47 s, ce [RUN]) | faster and slower polarisation | expected < 0.5 % except bursts |
| Constant current instead of constant power | load model | +1–3 % (overstates TTE near empty) |
| Brightness law S-g | convex slider law | dark-mode benefit at b = 1 vs 0.4 |
| Explicit RRC tail state r(t) | network tail dynamics | standby ±10 % depending on burst clustering |
| Exponential vs lognormal sessions | usage law | variance ×1.6, mean unchanged |
| Homogeneous vs diurnal CTMC | usage law | mean and quantile shifts |
| Entropic OCV off; R₀(z) flat (A_R = 0) | small terms | < 0.5 % at 25 °C; low-SOC cold cut-off |

**7.2.5 Usage fluctuations.**

1. TTE distribution at fixed θ (F6), and the CV law of §5.4.3 vs Monte Carlo when all d_j are scaled by
   κ ∈ [0.25, 4] (F7). The predicted CV ∝ √κ.
2. Rate-parameter uncertainty (±30 %) propagated through the closed forms.
3. Fluid vs PDMP mean bias (Dynkin correction; < 2 %).
4. Session heavy tails (§5.4.6).

### 7.3 Aging plan

* **Model:** M4 (§1.4) with daily operator splitting, over 3 years, for the S10 user.
  * Daily discharge: E_day = 24 h × P̄, which is 0.58 FCE/day.
  * Charging overnight at 0.7C, with T = T_amb + 8 K during charging (ASSUMPTION, 5–12 K; KANG19 measured
    ≈ 37 °C while fast-charging an idle phone).
  * Daytime T = T_amb + 3 K.
* **Factors:** charging strategy × climate.
  * Strategies: charge to 100 % and hold overnight; OBC (hold at 80 % until 1 h before unplugging); 80 % limit.
  * Climates: T_amb = 25 or 32 °C.

[PROTO] SOH:

| Strategy, climate | 1 yr | 2 yr | 3 yr |
|---|---|---|---|
| 100 %, 25 °C | 0.871 | 0.820 | 0.780 |
| OBC, 25 °C | 0.873 | 0.823 | 0.784 |
| 80 % limit, 25 °C | 0.881 | 0.834 | 0.797 |
| 100 %, 32 °C | 0.834 | 0.767 | 0.715 |
| 80 % limit, 32 °C | 0.852 | 0.791 | 0.744 |

**Finding:** heat costs 5.3 percentage points of SOH at 2 years, about four times the 80 % limit's 1.4–2.4 pp
gain. OBC gains only 0.3 pp in this model.

* **Outputs:**
  * TTE(SOH) for S2, S5, S7 and S10;
  * the cold-shutdown map s_off(T_amb, P_pk) at SOH = 1.0, 0.9 and 0.8 (F10);
  * the year in which the typical user's 5 % TTE quantile drops below 16 h ("no longer lasts a waking day").
* **Validation:** V4 (NASA). Aging shape: the √t law is the baseline; the t^0.75 calendar law (SCH14) is a
  structural variant.

### 7.4 Recommendations to users (ranking by TTE gain)

Compute each lever as ΔTTE on S10 (M0 closed form, then confirmed with M3 Monte Carlo), and on the scenario where
it matters most.

[PROTO] M0 on the typical day (base 40.4 h):

| Lever | ΔTTE (h) | ΔTTE (%) |
|---|---|---|
| Skip 3D-gaming sessions (5 % of sessions) | +10.5 | +26 |
| Prefer Wi-Fi: share 50 → 80 % | +2.5 | +6.2 |
| Cap refresh at 60 Hz | +1.7 | +4.3 |
| Dark mode at b = 0.4 | +1.4 | +3.4 |
| Brightness 0.4 → 0.3 | +0.6 | +1.4 |
| Halve background syncs | +0.25 | +0.6 |
| Avoided: weak signal (−110 dBm) | −3.2 | −8.0 |
| Avoided: 5G instead of 4G | −2.3 | −5.8 |
| Avoided: brightness 0.6 | −1.2 | −3.0 |

Scenario-specific levers:

* Video on Wi-Fi rather than 5G: +48 % (S3).
* Video call on Wi-Fi rather than 4G: +71 % (M0).
* Keep the phone warm in the cold (inner pocket) when bursts are expected: prevents shutdowns at 7–42 % displayed
  SOC (S7).
* Do not game in a hot car, which throttles performance by −43 %.
* Aging: avoid heat while charging; use an 80 % limit if a day's use is < 60 % of capacity.

**Presentation:** a ranked bar chart (F14) with error bars from θ-UQ, and a one-page user guide.

### 7.5 OS power-management strategies (model-grounded)

1. **Tail-aware coalescing of background traffic.** Batching N wake-ups changes the tail fraction from
   1 − e^{−λT} to 1 − e^{−λT/N}. With N = 5 on 4G, standby power falls from 0.052 to about 0.040 W (+30 % standby
   TTE).
2. **RAT steering.** Keep low-throughput states on 4G or Wi-Fi rather than 5G (XU20 dynamic switching −24.8 %).
   Value: the 4G vs 5G entries of §2.7; −5.8 % TTE avoided on the typical day.
3. **Signal-aware deferral.** Defer bulk transfers while RSRP < −100 dBm. The value is the integral of m_sig − 1
   over the deferred bytes.
4. **Adaptive refresh (LTPO policy).** 120 Hz only during touch-scroll: +4.3 % on the typical day.
5. **Thermal-aware frame capping instead of reactive throttling.** Because P ∝ f^1.9, a 20 % frequency cap cuts
   SoC power to 0.65× while losing 20 % performance. This keeps T_skin below T_thr and avoids throttle
   oscillation. Evaluate on S5: TTE and the mean ρ.
6. **Brown-out guard (peak-power governor).** Cap instantaneous demand at (1 − m)·P_max,cut(z, T, SOH), using the
   §5.3.4 closed form evaluated online with m = 0.1. Evaluate on S7b/S8c: s_off and TTE with and without the guard
   (LEE20 reported +38 % delivered capacity in the cold with current regulation; Apple performance management,
   AP-BATT).
7. **Model-based remaining-time display with uncertainty.** Show `TTE_hat = (E_rem(z, T, SOH) − h(x) + h̄_P)/P̄_EWMA`
   with an inverse-Gaussian 80 % band (§5.4.5) driven by the user's learned CTMC rates.
8. **Charging policy.** Use OBC or an 80 % limit plus thermal-aware charge current (JEITA zones, TI-JEITA). The
   value is the §7.3 table.

### 7.6 Generalisation to other portable devices

The structure carries over unchanged. Only constants and presets change. These scaling rules are DERIVED from the
physics; no device data were collected for other devices.

| Quantity | Scaling rule |
|---|---|
| Energy | E_u = Q·Ū·n_s for n_s cells in series (laptops); TTE = E_u/P̄ still holds |
| Pack resistance | R ∝ n_s/(n_p Q_cell) (R·Q constant per cell) |
| Display | P_disp ∝ A·g(b)·APL (OLED) or c₀ + c₁b + c₂b² (LCD, dp §2.2); tablets have 3–4× the area |
| Thermal | C_th ∝ mass, R_th ∝ 1/(h·A_surf), so P_sust ∝ (T_lim − T_amb)·h·A_surf: larger devices sustain more power |
| Watches | tiny C_th (τ about a minute), standby- and BLE-dominated, P̄ in mW |
| Usage | a device-specific activity set; the same PDMP formulas |

Examples:

* **Tablet:** display-dominated; Wi-Fi share higher; gaming less throttled.
* **Laptop:** CPU-dominated; multi-cell; constant power via DC-DC; cut-off per cell.
* **Smartwatch:** standby floor dominates (§5.3.1 elasticity → s_sleep ≈ 1).
* **E-reader:** bistable display, so energy per page turn replaces P_disp; the CTMC becomes a Poisson page-turn
  process and TTE ≈ E_u/(P_sleep + ν_page E_page).

---

## 8. Figures and tables for the paper

1. **F1, model schematic.** The M0 → M4 hierarchy: states, inputs, couplings, and the one mechanism each level
   adds.
2. **T1, parameter table.** §3, condensed to about 35 rows, with sources and UQ ranges.
3. **F2, power breakdown by activity.** Stacked bars of the §2.7 components for Wi-Fi / 4G / 5G, with V2
   benchmark markers.
4. **F3, SOC(t) for scenarios S1–S9.** M3, with throttled-power and temperature insets for S5 and S7.
5. **F4, TTE vs initial SOC.** s₀ = 10–100 % for S2, S3a, S5 and S10, comparing s and SOE.
6. **T2, structural sensitivity table.** ΔTTE of each model variant per scenario (§7.2.4).
7. **F6, stochastic TTE.** CTMC sample paths (activity raster + SOC), and the TTE histogram from Monte Carlo
   (M0 and M3) with the inverse-Gaussian density and the Dynkin mean.
8. **F7, CV of TTE vs session-length scale κ.** Analytic √(2τ_c/μ)·CV_P against Monte Carlo, for exponential and
   lognormal sessions.
9. **F8/F9, local sensitivity.** Elasticity heat map (parameters × scenarios) and OAT tornado (S10, S7a).
10. **F10, cold/aged shutdown map.** s_off over (T_amb, P_pk) at SOH 1.0 / 0.9 / 0.8, with LEE20 points overlaid.
11. **F11, thermal throttling.** S5 at 15 / 25 / 35 °C: T_skin(t), P(t), TTE and performance ρ.
12. **F12, validation.** Bonn predicted vs observed with 90 % intervals; the V2 benchmark band (±25 %); HG2
    f_Q(T) model vs data (M3 and M3 + M2-K).
13. **F13, aging.** SOH over 3 years for strategy × climate, and TTE vs SOH at 25 °C and −10 °C.
14. **F14 / T3, recommendations.** Ranked hours gained for user levers and OS strategies, with UQ error bars.
15. **F15, Sobol indices.** S₁/S_T for S2, S5, S7a and S10, with bootstrap CIs.

---

## 9. Implementation notes

### 9.1 Module map (`src/phonebatt/`)

| Module | Contents |
|---|---|
| `params.py` | Frozen dataclasses: `Param(value, lo, hi, dist, unit, source, conf)`. `RP63`, `IPHONE14P`, `IPHONE13P`, `HG2CELL`, … Function `sample(theta_spec, n, method="lhs")` |
| `ocv.py` | `U(z)`, `dUdz(z)`, `dUdT(z)` (table interpolation), `stretch(dV, p)` |
| `cell.py` | `R0(z, T, r)`, `R1(T, r)`, `tau1(T)`, `current(P, Ue, R)` (quadratic; returns NaN if infeasible), `Pmax_cut(z, T, r)`, `z_cut(P, T, r)` |
| `power.py` | Component functions (§2): `P_disp`, `P_bg`, `m_sig`, `P_gnss`, `state_power(activity, rat, **overrides) -> (P, components)`, `PRESETS` |
| `thermal.py` | `T_skin`, `throttle(P_soc_dem, T_skin)`, `perf_ratio` |
| `models.py` | RHS builders `rhs_M0`, `rhs_M1`, `rhs_M2`, `rhs_M2K`, `rhs_M3(theta, schedule)`. All are vectorised over samples, with state shape (n_state, n) |
| `usage.py` | `hub_generator`, `kron_sum`, `stationary`, `deviation_matrix`, `moments(L, P) -> (pi, Pbar, sigma2, h)`, `hub_closed_form`, `sample_paths(L, n, t_max, h)`, `dynkin_mean`, `ig_quantiles` |
| `tte.py` | `tte_M0`, `tte_M1_quadrature`, `tte_perturbation_a`, `simulate(model, theta, schedule, T_amb, z0, h)` wrapping `rk4_batch`, `reference_solve_ivp` |
| `aging.py` | `day_profile(strategy, T_amb, E_day)`, `advance_day(q, profile)`, `run_years` |
| `analysis/` | Scripts `run_scenarios.py`, `run_validation.py`, `run_uq.py`, `run_sensitivity.py`, `run_structural.py`, `run_aging.py`, `run_recommendations.py` → `results/*.json` and `paper/figures/*.pdf` |

### 9.2 Interfaces

```python
schedule(t: float) -> dict      # keys: activity or explicit components, rat, b, apl, fr, rsrp,
                                #       nu_fix, lam_bg, T_amb
simulate(model="M3", theta=RP63, schedule=..., z0=1.0, T0=None, h=2.0, t_max=400*3600)
    -> {"tte_h", "z_off", "s_off", "trace"}
moments(L, P) -> {"pi", "Pbar", "sigma2", "h", "mu_fluid", "mu_dynkin", "sd_clt"}
```

### 9.3 Regression tests (pytest; tolerances relative unless stated)

| Quantity | Expected [PROTO] | Tolerance |
|---|---|---|
| mean OCV over [0, 1] | 3.8699 V | 1e-4 |
| E_n | 18.189 Wh | 1e-4 |
| z_c0 | 0.0263 | 1e-3 abs |
| E_u0 | 17.788 Wh | 1e-3 |
| display area (6.3") | 97.4 cm² | 1e-3 |
| state powers (§2.7) | table values | 1e-3 |
| M1 TTE S2 | 15.63 h | 0.5 % |
| M1 TTE S3b | 19.57 h | 0.5 % |
| M1 TTE S6 | 6.26 h | 0.5 % |
| M3 TTE S5 at 25 °C | 5.90 h | 1 % |
| M3 TTE S5 at 35 °C | 11.09 h | 1 % |
| M3 S7b aged | TTE 10.67 h, s_off 42.4 % | 2 % |
| CTMC typical | P̄ = 0.4404 W, σ² = 0.2265 W²h, Dynkin mean 41.19 h | 1e-3 |
| hub closed form vs deviation matrix | equal | 1e-10 |
| Monte Carlo M0 mean vs Dynkin | – | ≤ 1 % at 20 000 paths |

### 9.4 Compute budget (numpy, one core)

* M3 deterministic run: < 1 s vectorised.
* 2000-path CTMC Monte Carlo on M3 at h = 5 s: about 2 min [PROTO].
* Sobol on the quasi-static surrogate (N = 1024, d ≈ 30, about 63 k evaluations): vectorise the quadrature with
  200-point Gauss–Legendre on [z_c, 1] and Newton iteration for z_c, for under 1 min per scenario.
* Full-M3 Sobol (N = 256, about 16 k runs): about 10 min per scenario.

### 9.5 Prototype provenance

`docs/design/prototype_parsimony/`. Run each script from inside that directory with `python3 -I <script>`; the
scripts add `.` to `sys.path`.

| Script | Produces |
|---|---|
| `proto.py` | OCV, energies, state powers, M0/M1 |
| `m3.py` | M3 deterministic scenarios and cold bursts |
| `ctmc.py`, `ctmc2.py` | CTMC analytics, Monte Carlo M0/M3 |
| `recs.py` | shares and levers |
| `aging.py` | SOH table |

Prototype simplifications to remove in production:

* In `ctmc*.py` and `recs.py` the game state uses fixed throttled powers in M0.
* `aging.py` uses an hourly day profile.

---

## Appendix A: derivations

**A.1 M1 quadrature.** From R_s I² − U I + P = 0, the smaller root is I = (U − √Δ)/(2R_s) with Δ = U² − 4R_sP.
Hence 1/I = 2R_s/(U − √Δ) = (U + √Δ)/(2P) = V_t/P, and dt = 3600 Q dz/I = 3600 Q V_t dz/P. Integrate from z_c
to z₀. The cut-off condition V = V_cut, with V = P/I, gives I = P/V_cut, so U(z_c) = V_cut + P R_s/V_cut.

**A.2 Perturbation constant a.** Expand V_t = U − R_sP/U + O(P²) and the cut-off shift
Δz_c = P R_s/(V_cut U′). Then ∫_{z_c(P)} V_t dz = ∫_{z_c0} U dz − P ∫ R_s/U dz − U(z_c0)·Δz_c + O(P²). With
U(z_c0) = V_cut, the last term is P R_s(z_c0)/U′(z_c0).

**A.3 Dynkin identity.** h = D P solves the Poisson equation −Λh = P − P̄1. For f = h, Dynkin's formula with the
bounded stopping time TTE ∧ t, letting t → ∞, gives E[h(X_TTE)] − h(x₀) = E∫₀^TTE (Λh)(X_s) ds =
E∫₀^TTE (P̄ − P_{X_s}) ds = P̄ E[TTE] − E_u. This holds in M0 because ∫₀^TTE P ds = E_u exactly.

**A.4 Hub variance.** Over one cycle, Y = (P₀ − P̄)D₀ + (P_J − P̄)D_J has E[Y] = 0, which gives
Σ_j p_j (P_j − P̄) d_j = −(P₀ − P̄)/λ_on. Then

```
E[Y²] = (P₀ − P̄)²(1 + CV₀²)/λ_on² + Σ_j p_j (1 + CV_j²) d_j² (P_j − P̄)² − 2(P₀ − P̄)²/λ_on²
```

The renewal-reward CLT gives σ² = E[Y²]/E[L]. With exponential sojourns (CV = 1) the idle terms cancel, which
gives §5.4.4.

**A.5 Inverse Gaussian.** The first passage of a Brownian motion with drift P̄ > 0 and variance σ² to level E_u
is IG(E_u/P̄, E_u²/σ²). Replace the mean by the Dynkin mean to correct the initial-state and depletion-state bias.
