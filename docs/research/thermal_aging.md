# Phone thermal behaviour, cold effects, aging and charging history: research note

Scope: parameters and model forms for the thermal, cold-temperature, aging (SOH) and charging-history parts of our continuous-time smartphone battery model (2026 MCM Problem A). The companion note `device_power.md` covers component power draw. The cell OCV/ECM and PyBaMM capacity-vs-temperature work is in `docs/research/data/cell_ecm_pybamm/`.

```
Fast states (hours):  SOC(t), V_RC(t), T(t)
  C_th dT/dt  = P_heat(t) - (T - T_amb)/R_th              P_heat = P_load + I^2 R_int (+ charger loss)
  R_int       = R_int(T, SOC, SOH)                        rises steeply in the cold, rises with age
  shutdown    : V_term = OCV(SOC) - I R_int - V_RC < V_cut  (V_cut ~ 3.4 V)
  throttling  : P_SoC capped once T_skin -> T_lim (~45 C)

Slow states (weeks to years):  q_loss(t), r_growth(t)
  dq_cal/dt   = f_cal(T, SOC)          (SEI growth, sqrt(t)-type, Arrhenius in T, rises with SOC)
  dq_cyc/dt   = f_cyc(T, I, SOC)       (charge-throughput driven; extra low-T / high-rate term = plating)
  Q_max = Q_0 (1 - q_loss),  R_int = R_0 (1 + r_growth)
```

Conventions:

* **Confidence.** H = read directly from the text or a table of a primary or peer-reviewed source, or an official manufacturer page. M = read off a figure (±10–20 %), taken from an abstract or a secondary restatement of a primary source, or from press coverage of an official statement. L = derived/extrapolated by us, background knowledge not verified in this session, or a non-peer-reviewed industry source.
* **DERIVED** rows are our own arithmetic on cited numbers. The arithmetic is shown so it can be re-checked; the script used is reproduced in §2.8.
* **Source keys** (e.g. `KANG19`) resolve to the numbered reference list in §5.
* Temperatures are in °C unless a formula needs kelvin. R_g = 8.314 J/(mol·K).

---

## 0. Executive summary: recommended parameter set

For a 2020–2025 phone (170–230 g, single-cell graphite/LCO-family battery, 4000–5000 mAh). Defaults are for model runs; ranges are for uncertainty quantification (UQ).

| Quantity | Symbol | Default | Range for UQ | Unit | Basis (details in §1) | Conf. |
|---|---|---|---|---|---|---|
| Lumped heat capacity of the phone | C_th | 150 | 100–230 | J/K | Three independent DERIVED estimates agree: mass × c_p; τ·R from HE19; Kang Fig. 7 energy/ΔT | M |
| Bulk (battery node) thermal resistance to ambient, per W of total device heat | R_th | 3 | 2–4 | K/W | DERIVED from HE19 Fig. 10; convection + radiation estimate | M/L |
| Skin hot-spot thermal resistance to ambient | R_skin | 5.5 | 4–9 | K/W | PU23 Pixel 5: 5 K/W; KANG19 Fig. 6: ~5–7 K/W; CHI15 (2.2–3.8 W at 45 °C skin) | M |
| Bulk thermal time constant | τ_th = R_th·C_th | 7.5 | 5–15 | min | HE19 Fig. 11 (exponential fit, Nexus 5X battery) | M |
| Skin hot-spot time constant | τ_skin | 3–6 | 2–6 | min | KANG19 restart fits 137–222 s; DERIVED 5.4–6.1 min from time-to-45 °C | M |
| Skin temperature limit (design/throttling target) | T_lim | 45 | 41–48 | °C | CHI15 (45 °C design limit); KANG19 (pain at 42–45 °C); IEC 62368-1 TS1 48 °C for contact > 1 min (UL-TS) | H |
| Sustainable total device power at T_amb = 25 °C | P_sust | 3.0 | 2.2–4.0 | W | CHI15 2.2→3.5 W (3.8 W ideal); DERIVED (45−25)/R_skin | M |
| Sustainable power vs ambient | P_sust(T_amb) | (T_lim − T_amb)/R_skin | — | W | DERIVED: ≈ 1.7–2 W at 35 °C ambient | L |
| Apple ambient operating range | — | 0 to 35 | — | °C | AP-TEMP | H |
| Charging temperature window (JEITA, handheld) | — | full rate 10–45; reduced 0–10 and 45–60; none < 0 or > 60 | — | °C | TI-JEITA | H |
| Phone-level battery resistance at 25 °C (incl. path) | R_int,25 | 0.15 | 0.10–0.33 | Ω | HE19 0.1 Ω (Galaxy S3 pulse); LEE20 Fig. 10: 0.13–0.33 Ω at shutoff | M |
| Resistance ratio R(−5 °C)/R(25 °C) at low SOC | — | 2.5 | 2.0–3.3 | × | LEE20 Fig. 10 (4 phones) | M |
| Effective Arrhenius energy of phone R_int | E_R | 18 | 12–27 | kJ/mol | DERIVED from LEE20 ratios | L |
| Phone cut-off voltage | V_cut | 3.4 | 3.0–3.5 | V | LEE20 ("approximately 3.4 V"; "typically below 3.4 V") | H |
| Displayed SOC at an unexpected cold shutoff | — | up to 33–35 % | 15–35 % (−5 °C); 30 % (−15 °C) | % | LEE20; HE19 | H |
| TTE penalty of 0 °C vs 25 °C (Nexus 5X, UI workload, no mitigation) | — | −15 % | — | — | LEE20: 159 → 135 min | H |
| Rated cycle life to 80 % capacity (Apple iPhone ≤ 14 / iPhone 15+) | N_80 | 500 / 1000 | — | full cycles | AP-BATT | H |
| Rated cycle life to 80 % capacity (Pixel 3–8 Pro / Pixel 8a+) | N_80 | 800 / 1000 | — | full cycles | GOOG | H |
| Calendar fade, Arrhenius energy (NMC 18650) | E_cal | 58 | 20.6 (LFP) – 58 | kJ/mol | SCH14 via CA21 (6976 K × R_g, DERIVED); SCM18 (LFP 20.6 kJ/mol) | M |
| Calendar fade acceleration per +10 K near 25 °C | — | ×2.1 | ×1.3 (LFP) – ×2.1 | × | DERIVED from E_cal | M |
| Calendar fade at 25 °C, mid SOC | q_cal(1 yr) | 2.5–3 % | 2–4 % | — | SCH14 model, 3.7–3.8 V (DERIVED); KEIL16 2–3 % after 10 months | M |
| Calendar fade at 25 °C, full charge | q_cal(1 yr) | ~4–5 % | 4–6 % | — | SCH14 at 4.1–4.2 V (DERIVED); KEIL16 NCA 5.5 % after 10 months at 100 % SOC | M |
| Resistance growth at 80 % SOH | R/R_new | 1.5 | 1.4–2.3 | × | DERIVED from SCH14 (+41–65 %). KEIL16 storage read-offs: NMC +48 % at 23 % fade; NCA +83 % at only 12.5 % fade (50 °C, 100 % SOC) | L |
| Temperature of minimum cycle aging | T_opt | 25 | ~20–30 | °C | WAL14 (Arrhenius plot with two branches crossing near 25 °C); WU17 (25 °C mildest) | M |
| Cell reversible self-discharge | k_sd | 0.1 | 0.05–2 | %/month | ZIL19 4.2 µA on 3.5 Ah (DERIVED 0.09 %/month); BU802b 1–2 %/month + 3 %/month protection circuit | M/L |

**Bottom line for the model.** (i) The phone heats with τ ≈ 5–15 min to a steady rise of ≈ 3–6 K per watt, so thermal throttling caps sustained power at ≈ 3 W. (ii) Cold mainly acts through R_int, which roughly doubles to triples between 25 °C and −5 °C. It causes shutoffs at a displayed SOC of 15–35 % under bursty load; most of the "lost" charge comes back once the phone is warm. (iii) Aging cuts Q_max by 20 % over 500–1000 full cycles and raises R_int by ≈ 40–100 %. Heat speeds up aging (≈ ×2 per +10 K), and so do long stays at high SOC. (iv) Self-discharge is negligible next to standby drain.

---

## 1. Detailed parameter tables

### 1.1 Lumped thermal parameters of the phone

| # | Name | Symbol | Value | Range | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|---|---|
| T1 | Skin-to-ambient thermal resistance (steady state) | R4 = (T_skin − T_amb)/P_tot | 0.005 °C/mW = 5 K/W | — | K/W | Pixel 5. Linear regression of skin temperature against total power (mmWave 5G + CPU loads) | PU23 | H |
| T2 | Total power in mmWave tests that drive overheating | P_tot | 5.7 (2 Gbps DL only); 6.3–8.8 (CPU + transceiver) | — | W | Pixel 5 | PU23 | H |
| T3 | Hot-spot steady temperature, Nexus 5X (hottest region, IR camera) | T_ss | Skype 58.1; Hangouts 62.3; Abyssrium 45.9; PUBG 45.2; video recording 50.8; voice call 33.0 | — | °C | Ambient 24–28 °C; 30-min runs; Table 1 | KANG19 | H |
| T4 | Same, across 9 phones (Nexus 5/6/5X, Pixel, Pixel 2, Galaxy S7, P20 lite, iPhone 7/7+/8) | T_ss | video chat 39.7–62.3; 3D game 33.6–46.2; video rec. 40.1–51.6; voice 31.9–38.8 | — | °C | Same | KANG19 Table 1 | H |
| T5 | Whole-device power during these runs (Nexus 5X) | P | Hangouts ≈ 6; Skype ≈ 3.6; Abyssrium ≈ 3.5; PUBG ≈ 2.8 | ±0.5 | W | Read off Fig. 6 (battery V·I via /proc) | KANG19 | M |
| T6 | Hot-spot R from T3/T5 | R_skin | 5.2 (Hangouts), 6.7 (Skype), 5.1 (Abyssrium), 6.8 (PUBG) | 5–7 | K/W | DERIVED: (T at 1800 s − 26 °C)/P, Fig. 6 | KANG19 | M/L |
| T7 | Time to steady state | — | within 20 min (all phones); within 30 min (all apps) | — | min | — | KANG19 | H |
| T8 | Time for Nexus 5X / Galaxy S7 to reach steady state, 3D game | — | 352 s / 606 s | — | s | User-study game (Yokai Saga) | KANG19 | H |
| T9 | Time for skin to reach 45 °C, Skype, Nexus 5X, without / with 1.5 mm TPU+PC case | — | 418 (SD 22) / 643 (SD 46) | — | s | Steady states 54 °C / 48 °C | KANG19 | H |
| T10 | First-order τ implied by T9 | τ_skin | 368 (no case) / 323 (case) | — | s | DERIVED: τ = −t / ln(1 − (45−26)/(T_ss−26)) | KANG19 | L |
| T11 | Exponential rate after app restart (fit T(t) = T∞ − (T∞−T0)e^(−t/τ)) | 1/τ | video 45→62 ×10⁻⁴ s⁻¹; PUBG 61→73 ×10⁻⁴ s⁻¹ | — | s⁻¹ | Nexus 5X; i.e. τ ≈ 137–222 s (DERIVED) | KANG19 | H (rates) |
| T12 | Battery temperature settling (exponential fit, Nexus 5X in a thermal chamber) | τ_bulk | ≈ 7.5 | 7–8 | min | DERIVED from Fig. 11: 17 K → 9 K at 5 min, → 4.5 K at 10 min. 34 traces fitted with adj. R² close to 1 (Fig. 12) | HE19 | M |
| T13 | Battery temperature vs current in a chamber (−13 to 30 °C, 27 runs, ≥ 1 h each) | ΔT_bat | ≈ 5–15 K above ambient at I² ≈ 0.6–1.25 A² | — | K | Read off Fig. 10. T_bat is linear in I², and T_amb adds an offset | HE19 | M |
| T14 | Bulk R_th per W of whole-device power | R_th | ≈ 2–4 | — | K/W | DERIVED: ΔT_bat/(3.8 V · I) from T13; assumes all device power ends up as heat in the phone | HE19 | L |
| T15 | Cold soak: idle Nexus 5X placed at −12 °C | — | battery → −11.8 °C; transient ≈ 60 min; R → ≈ 0.66 Ω | — | — | Fig. 6 and text | HE19 | H |
| T16 | Apparent heat capacity from energy vs surface temperature | E/ΔT | ≈ 150–170 | — | J/K | DERIVED from KANG19 Fig. 7: ≈ 2000 J raises hot spot from ≈ 28 to ≈ 40–41 °C. Losses and hot-spot non-uniformity bias this in opposite directions | KANG19 | L |
| T17 | Li-ion cell specific heat (Sony US-18650, LCO) | c_p,cell | 0.96 ± 0.02 (2.75 V OCV); 1.04 ± 0.02 (3.75 V) | — | J/(g·K) | Calorimetry | MAL99 | H |
| T18 | Aluminium specific heat (literature value used by modellers) | c_p,Al | 900 | — | J/(kg·K) | Cited as a literature value | SS18 | M |
| T19 | Effective specific heat of a phone (glass, Al/steel frame, PCB, battery) | c_eff | 0.9 | 0.75–1.0 | J/(g·K) | Glass/steel/plastic values not verified this session | — | L |
| T20 | Phone mass | m | 190 | 150–230 | g | Typical modern spec sheets (not fetched). iPhone 6s teardown: 142 g total, 25 g battery (search snippet only) | — | L |
| T21 | Lumped heat capacity m·c_eff | C_th | 171 | 112–230 | J/K | DERIVED: 0.15–0.23 kg × 750–1000 J/(kg·K) | — | L |
| T22 | C_th from τ·R | C_th | 150 | 112–240 | J/K | DERIVED: τ = 450 s (T12) with R_th = 3 K/W (2.5–4) | HE19 | L |
| T23 | Radiative heat-transfer coefficient (ε = 0.9, 35 °C surface) | h_rad = 4εσT³ | 6.0 | 5.5–6.5 | W/(m²·K) | DERIVED (physics; ε assumed) | — | L |
| T24 | Natural-convection coefficient | h_conv | 5–10 | 5–25 | W/(m²·K) | General range quoted in a smartphone thermal-simulation context (search snippet; source not opened) | — | L |
| T25 | Exposed area of a 150 × 72 × 8 mm slab | A | 0.025 | 0.02–0.03 | m² | DERIVED geometry | — | L |
| T26 | Isothermal-phone R_th = 1/((h_conv + h_rad)A) | R_th | 2.5–3.6 | — | K/W | DERIVED: consistent with T14. Skin hot spots sit higher (T1, T6) because the phone is not isothermal | — | L |
| T27 | Steady-state wait used in thermal design testing | — | 40 | — | min | Qualcomm/Stanford | CHI15 | H |

Notes on thermal coupling:

* **Battery self-heating is a small term.** Battery Joule heat I²R ≈ (1 A)² × 0.15 Ω ≈ 0.15 W, while the device dissipates 2–6 W (T5). For the battery temperature, the driving term is whole-device power through the shared chassis. HE19 attribute heat to I² only; their I² regressor is a proxy for device power.
* **Hands, pockets and cases add a term.** They change R_th (a case slowed the rise to 45 °C by 54 %, T9). Model them as a multiplicative factor on R_th and C_th for scenario studies. We have no sourced multiplier; treat it as a UQ range of 0.8–1.5× (L).

### 1.2 Thermal limits, throttling and safety thresholds

| # | Name | Value | Unit | Conditions / quote | Source | Conf. |
|---|---|---|---|---|---|---|
| L1 | Warm-sensation threshold | 33–35 | °C | Human perception studies cited in KANG19 | KANG19 | H |
| L2 | Pain onset | 42–45 | °C | "starts to feel pain at around 42–45 °C". Basis of EN563 limits (cited) | KANG19 | H |
| L3 | Skin-temperature design limit used in phone thermal design | 45 | °C | — | CHI15 | H |
| L4 | IEC 62368-1 TS1 max touch temperature, contact > 1 min | 48 (metal and plastic) | °C | TS2 = TS1 + 10 K. 10 s–1 min: metal 51, plastic 60. Based on ISO 13732-1:2006; tested at 25 °C ambient | UL-TS | H |
| L5 | Sustainable heat dissipation at the 45 °C skin limit | 2.2 → 3.5 (larger spreader); 3.8 (ideal spreading); +1.2 W from good vs poor spreading | W | Simulation and test of phone designs | CHI15 | H |
| L6 | Skin-temperature sensitivity of dissipation | 0.16 | W per °C | Same article (2023 reprint) | CHI15 | M |
| L7 | Android battery-temperature shutdown example | 68 | °C | "a device is shut down when the battery temperature exceeds the threshold (e.g., 68 °C)" | KANG19 | M |
| L8 | Android thermal status levels | NONE, LIGHT, MODERATE, SEVERE, CRITICAL, EMERGENCY, SHUTDOWN | — | Skin, battery, GPU, CPU and USB sensors; "device skin temperature is the most important system to track"; example critical threshold 65–72 °C (illustrative) | AOSP-TH | H |
| L9 | Apple behaviour when hot | Charging (incl. wireless) slows or stops; display dims; reduced frame rates or longer processing | — | — | AP-TEMP | H |
| L10 | Apple: limit charging above 80 % when warm | "Software may limit charging above 80% when the recommended battery temperatures are exceeded" | — | — | AP-MAX | H |
| L11 | Apple: ambient > 35 °C | "can permanently damage battery capacity"; ideal comfort zone 16–22 °C | °C | — | AP-MAX | H |
| L12 | Google: avoid prolonged heat | Avoid "6 or more hours" at "more than 35 °C" ambient; charge at about 25 °C | — | — | GOOG | H |
| L13 | Throttling cost on performance | DVFS frequency reduction / core shut-off. Nexus 5X 1.8 GHz cores switched off above threshold; AnTuTu score falls with ambient temperature | — | Figure-level only | KANG19 | M |
| L14 | Throttling performance loss up to 34 % | See device_power.md (BH19) | % | Cross-reference, not re-verified here | BH19 (device_power.md) | M |

### 1.3 Cold behaviour: resistance, shutoffs, recoverable capacity

| # | Name | Value | Range | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|---|
| C1 | Apple statement on cold | Outside 0–35 °C cold "might temporarily shorten battery life and could cause your device to turn off"; "battery life will return to normal" when warmed | — | — | — | AP-TEMP | H |
| C2 | Apple: impedance and cold | "A battery's impedance will temporarily increase at a low state of charge and in a cold temperature environment"; aged batteries have higher impedance and lower peak power | — | — | — | AP-BATT, AP-2017 | H |
| C3 | Phone battery resistance vs temperature (Nexus 5X) | ≈ 0.22 (40 °C), ≈ 0.25 (20 °C), ≈ 0.3 (0 °C), 0.3–0.36 (−10 to −15 °C) | — | Ω | Fig. 7, during use. Fitted 2-exponential law r(T) = a₁e^(b₁T) + c₁e^(d₁T) (RMSE < 0.03 Ω) | LEE20, HE19 | M |
| C4 | R0 vs temperature (Nexus 5X, 100 cycles of data) | ≈ 0.7 (−20 °C), ≈ 0.35 (0 °C), ≈ 0.25 (20 °C), ≈ 0.2 (40 °C) | — | Ω | Fig. 15a | LEE20 | M |
| C5 | Resistance at shutoff, 25 °C → −5 °C | Nexus 5X 0.13 → 0.43; Nexus 6P 0.33 → 0.67; iPhone 5S 0.26 → 0.55; iPhone SE 0.20 → 0.57 | — | Ω | Fig. 10a | LEE20 | M |
| C6 | Effective Arrhenius energy of R_int | 15 (Nexus 5X R0 −20/25 °C ratio 2.8 → 14.4); 12–27 (shutoff ratios 1.7–3.3 at −5 °C) | 12–27 | kJ/mol | DERIVED: E = R_g ln(ratio)/(1/T − 1/298.15). Below ≈ 5 °C the rise is steeper than a single Arrhenius law (HE19 text) | LEE20, HE19 | L |
| C7 | Typical instantaneous phone battery resistance (example) | 0.1 | — | Ω | Galaxy S3: 1 A → 0.5 A step, ΔV = 0.05 V | HE19 | H |
| C8 | Phone cut-off voltage | ≈ 3.4 | 3.0–3.5 | V | "phone shuts off when the voltage drops to about 3.4V"; cut-off of 8 phones "typically below 3.4V" | LEE20 | H |
| C9 | Max SOC at an unexpected shutoff, Nexus 5X | up to 33 | 0–33 | % | −5 to 25 °C ambient, I_max 1–2 A (Fig. 9b) | LEE20 | H |
| C10 | SOC at shutoff at −5 °C, iPhone 5S vs iPhone SE (with iOS 10.2.1) | 35 vs 15 | — | % | — | LEE20 | H |
| C11 | Voltage drop before shutoff vs after warming | 0.49 → 0.14 | — | V | Nexus 5X, −5 °C then room temperature | LEE20 | H |
| C12 | Capacity recovered after warming the shut-off phone | +330 mAh (+18 min of use) | — | mAh | Nexus 5X, no recharge | LEE20 | H |
| C13 | Delivered capacity in cold before shutoff (no mitigation) | 1897 mAh; +730 mAh (+38.4 %) with current regulation (BPM) | — | mAh | Nexus 5X, cold | LEE20 | H |
| C14 | Operating time 25 °C → 0 °C (UI exerciser, 10 runs) | 159 → 135 min (−15.1 %) without mitigation; 163 → 154 min (−6.1 %) with BPM | — | min | Nexus 5X | LEE20 | H |
| C15 | Shutoff at −15 °C | Xperia Z shut off showing 30 % SOC; Nexus 5X fuel gauge "could be over 90%" at shutoff vs ≈ 3 % with a temperature-aware estimate | — | % | 35 case studies | HE19 | H |
| C16 | User prevalence | 60 % of 121 surveyed users had experienced unexpected shutoffs; 43 % of iPhone users still affected after iOS 10.2.1 | — | % | Survey | LEE20 | H |
| C17 | iOS 10.2.1 effect | 80 % fewer unexpected shutdowns on iPhone 6s, 70 % on iPhone 6 | — | % | Apple statement, Feb 2017, via press | PRESS-17 | M |
| C18 | iPhone 6s shutdown root cause (2016) | Battery component of Sept–Oct 2015 units over-exposed to air → faster degradation; iPhone "designed to shut down automatically… such as extremely cold temperature" | — | — | Apple China statement, Dec 2016, via press | PRESS-16 | M |
| C19 | Apple performance management | iPhone 11 and later: iOS "dynamically manages performance peaks to prevent the device from unexpectedly shutting down" (longer app launches, lower frame rates, dimmer backlight) | — | — | — | AP-BATT | H |

### 1.4 Aging: manufacturer statements

| # | Statement | Value | Source | Conf. |
|---|---|---|---|---|
| A1 | iPhone 14 and earlier: "retain 80 percent of their original capacity at 500 complete charge cycles" under ideal conditions | 500 cycles → 80 % | AP-BATT | H |
| A2 | iPhone 15 and later: "retain 80 percent of their original capacity at 1000 complete charge cycles" under ideal conditions | 1000 cycles → 80 % | AP-BATT | H |
| A3 | Charge-cycle definition: one cycle = discharging an amount equal to 100 % of capacity, "but not necessarily all from one charge" | — | AP-WHY | H |
| A4 | Pixel 3 to Pixel 8 Pro and Pixel Fold: "up to 80% capacity for about 800 charge cycles". Pixel 8a and later: about 1000 | 800 / 1000 cycles → 80 % | GOOG | H |
| A5 | Pixel Battery Health Assistance "adjusts the battery's maximum voltage in stages". Press reports say the stages run from 200 to 1000 cycles (not on the Google page we read) | — | GOOG (press for stage limits) | H / M |
| A6 | Apple storage advice: store at ≈ 50 %. Fully charged storage "may lose some capacity"; fully discharged storage risks "deep discharge state" | — | AP-MAX | H |

### 1.5 Aging: laboratory measurements and model parameters

| # | Name | Value | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|
| G1 | Calendar fade law (NMC/graphite 18650, Sanyo UR18650E, 2.05 Ah) | q_cal = α_cap(V,T)·t^0.75, α_cap = (7.543·V − 23.75)·10⁶·exp(−6976/T) | t in days, V in volts (storage voltage), T in K | > 60 cells; 3 temperatures, various SOC | SCH14 via CA21 | M |
| G2 | Calendar resistance growth | r_cal = α_res·t^0.75, α_res = (5.270·V − 16.32)·10⁵·exp(−5986/T) | same | same | SCH14 via CA21 | M |
| G3 | Cycle fade law | q_cyc = β_cap·√Q, β_cap = 7.348·10⁻³(V̄ − 3.667)² + 7.600·10⁻⁴ + 4.081·10⁻³·ΔDOD | Q = charge throughput in Ah; V̄ = quadratic-mean voltage; ΔDOD as a fraction (our reading, see caveats) | Various DOD and mean SOC | SCH14 via CA21 | M |
| G4 | Cycle resistance growth | r_cyc = β_res·Q, β_res = 2.153·10⁻⁴(V̄ − 3.725)² − 1.521·10⁻⁵ + 2.798·10⁻⁴·ΔDOD | same | — | SCH14 via CA21 | M |
| G5 | Activation energies implied by G1/G2 | E_cal = 6976·R_g = 58.0; E_res = 5986·R_g = 49.8 | kJ/mol | DERIVED | SCH14 | M |
| G6 | Calendar fade after 1 yr (G1, DERIVED) | 25 °C: 2.4 % (3.7 V), 3.7 % (4.0 V), 4.6 % (4.2 V). 35 °C: 5.1 / 7.9 / 9.8 %. 45 °C: 10.4 / 16.1 / 19.9 % | % | Model extrapolation; 4.2 V is above the cell's tested window (see caveats) | SCH14 | L |
| G7 | Cycles to 80 % (G3 only, V̄ = 3.7 V, DERIVED, Q counted as discharge Ah) | 830 (100 % DOD), 1200 (80 %), 2470 (50 %), 7780 (20 %) | FCE | No calendar term; factor-2 ambiguity in Q definition | SCH14 | L |
| G8 | Resistance growth when 20 % capacity is lost (DERIVED) | cycling: +45 % (100 % DOD) to +65 % (20 % DOD); calendar: +41–42 % | % | From G1–G4 | SCH14 | L |
| G9 | Calendar fade after 10 months vs storage SOC: NCA (Panasonic NCR18650PD 2.8 Ah) | 25 °C: ≈ 2 % (0 % SOC), ≈ 3 % (50 %), ≈ 5.5 % (100 %). 40 °C: ≈ 3 / 5 / 9.5 %. 50 °C: ≈ 3.5 / 7.5 / 12.5 % | % | Read off Fig. 2a. Plateaus over SOC ranges of 20–30 %; step at ≈ 60 % SOC | KEIL16 | M |
| G10 | Same, NMC (Sanyo UR18650E) | 25 °C: ≈ 1.5 / 2.5 / 6.5 %. 50 °C: ≈ 3 / 12 / 23 % (0 / 50 / 100 % SOC) | % | Fig. 2b | KEIL16 | M |
| G11 | Resistance rise after 10 months, NCA | 25 °C: +15 % (0 % SOC) → +35 % (100 %). 50 °C: +17 % → +83 % | % | Fig. 2d. No direct correlation between capacity fade and resistance rise | KEIL16 | M |
| G12 | LFP calendar fade rate after 9 months | ≈ 0.2 (25 °C), ≈ 0.5 (50 °C) | %-points/month | Text | KEIL16 | H |
| G13 | Mechanism of SOC dependence | Low graphite potential (lithiation > 50 %) accelerates SEI growth and lithium loss. Cathode high-voltage effects "seem to play only a minor role". "High storage SoCs… should be avoided" | — | 16 SOCs, 3 temperatures | KEIL16 | H |
| G14 | Calendar model with SOC via anode potential (LFP, Sony US26650FTC1, 3 Ah) | Q_L,cal = k_cal(T,SOC)·√t, k_cal,ref = 3.69·10⁻⁴ h⁻⁰·⁵ (25 °C, 50 % SOC), E_a,cal = 20.6 kJ/mol, α = 0.384, k₀ = 0.142, U_a,ref = 0.123 V | — | Storage and cycle tests | SCM18 | H |
| G15 | High-temperature cycling term (LFP) | k_cyc,HighT,ref = 1.46·10⁻⁴ Ah⁻⁰·⁵, E_a = 32.7 kJ/mol (applied to √Q_tot) | — | 1C, 25 °C reference | SCM18 | H |
| G16 | Low-temperature (plating-type) cycling term (LFP) | k_cyc,LowT,ref = 4.01·10⁻⁴ Ah⁻⁰·⁵, E_a = 55.5 kJ/mol with **inverse** Arrhenius (faster when colder), β_LowT = 2.64 h on (I_ch − I_ref)/C₀, applied to √Q_ch | — | Charge-current dependent | SCM18 | H |
| G17 | Low-T, high-SOC term (LFP) | k_ref = 2.03·10⁻⁶ Ah⁻¹, E_a = 2.33·10⁵ J/mol (inverse), β = 7.84 h, active above SOC_ref = 82 %, linear in Q_ch | — | CC-CV tests | SCM18 | H |
| G18 | Low-T term acceleration (DERIVED from G16) | ×3.3 at 10 °C, ×7.8 at 0 °C (vs 25 °C, 1C charge) | × | — | SCM18 | L |
| G19 | Temperature optimum of cycle aging | Arrhenius plot with two regimes. Below 25 °C aging rises as T falls (Li plating); above 25 °C it rises with T (SEI growth, cathode degradation) | — | Cells cycled −20 to 70 °C | WAL14 | M |
| G20 | LCO/NCA blend pouch cell (Kokam 5 Ah) | 25 °C mildest; 10 °C and 40 °C both accelerate fade. 20 % loss at ≈ 3000–5000 EFC. Discharge rate 1C vs 5C: ≈ 4800 vs ≈ 3500 EFC. Ohmic R at EOL at 5C: 11.3 mΩ vs 2.7 mΩ new (> 3×) | — | 2.7–4.2 V | WU17 | H |
| G21 | LCO 18650 (Sony US18650S, 1.8 Ah) cycled at 25/45/50/55 °C | RT: 30 % fade after 800 cycles. 45 °C: 36 % after 800. 55 °C: ≈ 70 % after 490 | % | CC 1 A to 4.2 V, CV to 50 mA, 2.0–4.2 V | RAM02 (abstract, via search) + Popov et al. ECS abstract (protocol) | M |
| G22 | First-principles SEI side-reaction model (LCO/graphite) | J_s = −i₀s·a·exp(−α_c F η_s/(R T)), η_s = φ₁ − φ₂ − U_ref,s − J R_film/a. Film growth ∂δ/∂t = −J_s M_P/(a ρ_P F); R_film = R_SEI + δ/κ_P. U_ref,s = 0.4 V vs Li; i₀s = 1.5·10⁻⁶ A/m²; ρ_P = 2.1·10³ kg/m³; κ_P = 1 S/m; R_SEI,0 = 100 Ω·cm²; M_P printed as "7.3 × 10⁴ mol/kg" (units as printed) | — | Side reaction during charge only | RAM04 | H |
| G23 | Effect of end-of-charge voltage (simulated, 10 cycles) | 7.2 / 4.4 / 3.8 % fade for 4.2 / 4.0 / 3.9 V. "Overcharging the cell by 0.1 V… can reduce cycle life by up to 60%" (cited) | % | Simulation | RAM04 | H |
| G24 | LCO cycle-life drivers (900 mAh prismatic) | "Cycle-life is greatly influenced by the charge conditions, but is relatively insensitive to the discharge conditions". High charge cut-off voltage and long float at ≥ 4.2 V are most severe; > 1C charge damaging; lower DOD "does not improve the cycle-life" | — | Abstract | CHOI02 | M |
| G25 | Commercial LCO cells taken from phones | Upper cut-off 4.35 or 4.4 V is common; rated > 75 % retention after 300 cycles at 23–25 °C. At 40 °C with 24 h top-of-charge holds, cells at 4.35–4.4 V showed "severe capacity fade" | — | — | NEL18 | H |
| G26 | LFP cycle-life model | Q_loss = B·exp[(−31700 + 370.3·C_rate)/(R_g T)]·(Ah)^0.55, B depends on C-rate | % (Ah throughput) | A123 2.2 Ah, −30 to 60 °C, C/2–10C | WANG11 via FAR24 and search snippets | M |
| G27 | End-of-life definitions | "20% capacity loss or 30% internal resistance growth… whichever comes first" (vehicle practice) | — | — | YM12 | H |
| G28 | Field evidence (crowd-sourced Carat data) | 55 % of 9560 devices showed some full-charge-capacity loss; for some models median loss > 20 % with IQR > 20 pp | — | Android, crowd-sourced | HOQ16 | H (abstract) |
| G29 | Phone battery parameters vs cycle (Nexus 5X, low SOC 0–30 %) | R0 rises from ≈ 0.15 Ω (1st) to ≈ 0.2 (100th) and ≈ 0.3 Ω (200th); RC τ falls with cycling | Ω | Fig. 14 | LEE20 | M |
| G30 | Aged phones benefit more from current regulation | 300th-cycle battery: +42.8 % delivered capacity, +26 min with BPM. Operation time gains up to 2.03× (Nexus 6P, 414th cycle) | — | — | LEE20 | H |

### 1.6 Charging history and charging controls

| # | Name | Value | Source | Conf. |
|---|---|---|---|---|
| H1 | Apple Optimized Battery Charging | Learns the daily routine and holds charging at 80 % in some situations, "reducing the time your iPhone spends fully charged". On by default | AP-OBC | H |
| H2 | Apple Charge Limit (iPhone 15 and later) | 80–100 % in 5 % steps; OBC available only at 100 %. The phone "will occasionally charge to 100 percent to maintain accurate battery state-of-charge estimates" | AP-OBC | H |
| H3 | Apple charging profile | "fast charging to quickly reach 80% of its capacity, then switches to slower trickle charging" | AP-WHY | H |
| H4 | Apple fast-charge claims (page dated 22 Oct 2026) | 50 % in ≈ 30 min (20 W, iPhone 17e/Air); ≈ 20 min (40 W, iPhone 17 family); ≈ 15 min (60 W, iPhone 18 Pro). DERIVED average C-rate over the first half: 1C / 1.5C / 2C | AP-FAST | H (claims) / L (C-rate) |
| H5 | Apple: heat while charging | "Charging your device when it's inside certain styles of cases may generate excess heat, which can affect battery capacity" | AP-MAX | H |
| H6 | Samsung Battery protection (One UI 6.1+) | Basic: stop at 100 %, resume at 95 %. Maximum: stop at 80 %. Adaptive: Maximum while asleep, Basic before waking. Older One UI: single 85 % cap. "Keeping your battery at a full 100% charge for a long time can reduce its lifespan" | SAM | H |
| H7 | Google: use Adaptive Charging; charge at ≈ 25 °C with USB-C PD/PPS | — | GOOG | H |
| H8 | JEITA zones, single-cell handheld | T1 = 0, T2 = 10, T3 = 45, T4 = 50, T5 = 60 °C. 10–45 °C: up to 1C, 4.25 V max (4.2 V typ.). 45–50 °C: ≤ 4.15 V. 50–60 °C: ≤ 4.10 V. Below 0 °C or above 60 °C: no charge. Example charger IC halves current at 0–10 °C and drops to 4.06 V at 45–60 °C. Cold-charging rationale: lithium deposits on the anode → permanent loss | TI-JEITA | H |
| H9 | Measured surface temperature while (fast) charging, phone idle | ≈ 37 °C (2 % → 100 % in 100 min) | KANG19 | H |
| H10 | Charging while running video chat at 0–7 % SOC | Charge rate rises once SOC hits 0 % → fingerprint-area temperature 72 °C, 14 K above the non-charging case | KANG19 | H |
| H11 | User charging behaviour (> 80 000 Android phones, 4 weeks) | Start-of-charge SOC is spread evenly. Half of charge events end above 90 %. > 80 % of charge events last < 2 h (mean 1.2 h). Mean interval between charges 27.3 h. 70 % of users charge at least daily | GUO17 | H |
| H12 | Effect of fast charging on phone aging | Informal tests: ≈ 0.3–0.5 pp extra loss over 500 cycles (iPhone 12, Android 18 W vs 120 W-class). Not peer-reviewed. Lab evidence: higher charge current and lower temperature drive plating (G16, G24, TOM19) | press/YouTube; SCM18; CHOI02; TOM19 | L / H |

### 1.7 Self-discharge

| # | Name | Value | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|
| S1 | Reversible self-discharge current, 24 Ni-rich/Si-graphite 18650 cells (3.5 Ah) | 4.2 µA mean (relative variance 10 %) | µA | 70 % SOC, 25 °C, after 10 months of storage | ZIL19 (abstract) | H |
| S2 | Same as a rate | 0.088 | %/month | DERIVED: 4.2 µA × 730 h / 3.5 Ah | ZIL19 | L |
| S3 | Industry rule of thumb | ≈ 5 % in the first 24 h, then 1–2 %/month, "the protection circuit adds another 3 percent per month" | %/month | Cadex/Battery University | BU802b | L |
| S4 | Same source, table "per month" | Full charge: 6 % (0 °C), 20 % (25 °C), 35 % (60 °C). 40–60 % charge: 2 / 4 / 15 % | %/month | Inconsistent with S3; likely includes irreversible capacity loss. Do not use without checking | BU802b | L |
| S5 | Relative statement | Li-ion self-discharge "10 times lower than that of NiMH" | — | — | TI-JEITA | M |
| S6 | Comparison with phone standby drain | Standby 30 mW on a 19.3 Wh battery ≈ 3.7 %/day, vs self-discharge 0.003–0.07 %/day | %/day | DERIVED; P_sleep from device_power.md; 5000 mAh × 3.85 V nominal (unverified) | — | L |

---

## 2. Model forms and equations

### 2.1 Lumped thermal model (one node)

```
C_th · dT/dt = P_heat(t) − (T − T_amb)/R_th                                  (T1)
P_heat = P_load(t) + I(t)² R_int(T, SOC, SOH) + P_chg,loss(t) [− I T ∂U_oc/∂T]
```

* **P_load.** Power drawn at the battery terminals by the device (device_power.md). Nearly all of it ends up as heat inside the phone. Emitted light and RF are a small share; we did not quantify it.
* **I² R_int.** Battery internal loss. The Bernardi energy balance (BER85) adds the reversible entropic term −I·T·∂U/∂T. SS18 found the reversible heat had "low influence" for their cell, so we drop it by default.
* **P_chg,loss.** Charger and PMIC loss when plugged in. We did not find a sourced efficiency; use (1 − η)·P_in with η ∈ [0.85, 0.95] (L) for UQ.
* **Parameters.** C_th ≈ 150 J/K, R_th ≈ 3 K/W, so τ = R_th C_th ≈ 450 s. This matches HE19's 7.5-min exponential settling. The analytic step response is

  T(t) = T_amb + R_th P + (T₀ − T_amb − R_th P)·e^(−t/τ)   (HE19 Eq. 9).

* **Skin hot spot.** We use a quasi-steady offset rather than a second node:

  T_skin ≈ T_amb + (T − T_amb)·(R_skin/R_th), with R_skin/R_th ≈ 1.5–2.5 (DERIVED: 5–7 K/W ÷ 3 K/W).

  A two-node version (hot spot / bulk) can be fitted to KANG19's restart time constants (2–4 min). We do not have sourced node capacities, so it would be calibration, not data.

### 2.2 Thermal throttling (sustained power cap)

Steady state requires P_heat ≤ P_sust(T_amb) = (T_lim − T_amb)/R_skin.

Default: T_lim = 45 °C, R_skin = 5.5 K/W, so P_sust ≈ 3.6 W at 25 °C. This sits inside CHI15's 2.2–3.8 W.

A smooth, continuous-time cap on the SoC share of power (a modelling choice, not sourced):

```
P_SoC,eff = P_SoC,floor + (P_SoC,dem − P_SoC,floor) · σ((T_thr − T_skin)/w),   σ(x) = 1/(1+e^(−x))
```

With T_thr ≈ T_lim − 2 K and w ≈ 1 K this behaves like Android's graded LIGHT → SEVERE throttling (AOSP-TH). A consequence (DERIVED): at 35 °C ambient P_sust falls to ≈ 1.7–2 W, so in hot weather a capped game draws less power (TTE can lengthen) at the cost of performance. Leakage power also rises with temperature (KANG19), which works the other way. At 0 °C ambient P_sust ≈ 7–9 W, so throttling rarely binds in the cold.

### 2.3 Temperature-dependent resistance, cold shutdown and recoverable capacity

Recommended resistance law (SS18 form with a temperature-independent part):

```
R_int(T, SOC, SOH) = [R_∞ + R_a · exp( (E_R/R_g)(1/T − 1/T_ref) )] · f_SOC(SOC) · (1 + r_growth)      (T2)
```

* **Default.** One effective Arrhenius term with E_R ≈ 15–20 kJ/mol (DERIVED from LEE20 phone data; range 12–27). This reproduces R(−5 °C)/R(25 °C) ≈ 2–3 and R(−20 °C)/R(25 °C) ≈ 2.8.
* **Data-driven alternative.** HE19 fitted r(T) = a₁e^(b₁T) + c₁e^(d₁T) with RMSE < 0.03 Ω. They did not publish the coefficients.
* **f_SOC.** Resistance rises at low SOC (C2, LEE20 Fig. 14). Take its shape from the cell ECM note.

**Shutdown condition (constant-power load).** The terminal voltage is V_term = OCV(SOC) − I R_int − V_RC, with I = P_load / V_term. The phone shuts off when V_term < V_cut ≈ 3.4 V. The largest power deliverable without hitting V_cut is

```
P_max,cut = V_cut · (OCV − V_cut − V_RC) / R_int                                          (T3)
```

DERIVED examples with V_RC = 0:

| OCV | R_int | P_max,cut | Situation |
|---|---|---|---|
| 3.7 V | 0.15 Ω | 6.8 W | Warm |
| 3.7 V | 0.6 Ω | 1.7 W | Cold |
| 3.95 V | 0.6 Ω | 3.1 W | Cold, high SOC |

A cold phone at mid-to-low SOC therefore cannot serve a 2–3 W burst. This matches LEE20's "may shut off with an OCV as high as 3.95 V". Separately, a solution of the constant-power equation exists only while P ≤ OCV²/(4R_int); in practice V_cut binds first.

**Recoverable capacity.** Cold does not destroy charge. The charge left when V_cut is hit at low T can still be delivered once T rises (C12: +330 mAh). In the model this is automatic: SOC is not reset at shutoff, so restarting the ODE at a higher T yields further TTE. Self-heating also helps: at 3 W and R_th = 3 K/W the battery runs ≈ 9 K above ambient (DERIVED; HE19 Fig. 10 shows −3 °C battery at −13 °C ambient under ≈ 1 A).

### 2.4 Aging laws in continuous time

Separate time scales: SOH changes over weeks, a discharge lasts hours. Within one TTE simulation treat q_loss and r_growth as constants. Integrate the aging ODEs over the usage history (days to years).

**(a) Semi-empirical NMC model (SCH14) in rate form.** Use the "rate-based integral" approach recommended by SCM18 (after Thomas et al.) so that time-varying T and SOC are handled without an explicit t:

```
q_cal = α_cap t^0.75     ⇒  dq_cal/dt = 0.75 · α_cap(V,T)^(4/3) · q_cal^(−1/3)        [t in days]
r_cal = α_res t^0.75     ⇒  dr_cal/dt = 0.75 · α_res(V,T)^(4/3) · r_cal^(−1/3)
q_cyc = β_cap √Q         ⇒  dq_cyc/dt = β_cap(V̄,ΔDOD)² · |I| / (2 q_cyc)              [I in A, t in h]
r_cyc = β_res Q          ⇒  dr_cyc/dt = β_res(V̄,ΔDOD) · |I|
```

Coefficients are G1–G4. Start from q(0) = ε (e.g. 10⁻⁴) to avoid the singular initial rate, or start from the phone's actual age.

To transfer to a phone cell:

1. Map by SOC, not voltage. Use the NMC cell's OCV at the phone's SOC as V, because phone LCO-family cells reach 4.35–4.45 V (NEL18).
2. Scale throughput as Q_equiv = (phone Ah throughput) × 2.05 Ah / C_phone.

**(b) Mechanism-separated form (SCM18 structure, LFP-fitted).** This is the most useful structure for a phone model, because it separates heat-driven and cold/fast-charge-driven damage:

```
k_cal(T,SOC)  = k_cal,ref · exp[−(E_cal/R_g)(1/T − 1/T_ref)] · {exp[α F (U_a,ref − U_a(SOC))/(R_g T_ref)] + k₀}
k_H(T)        = k_H,ref  · exp[−(E_H/R_g)(1/T − 1/T_ref)]                           (SEI/cathode, hot)
k_L(T,I_ch)   = k_L,ref  · exp[+(E_L/R_g)(1/T − 1/T_ref)] · exp[β_L (I_ch − I_ref)/C₀]  (plating, cold/fast)
Q_loss        = k_cal √t + k_H √Q_tot + k_L √Q_ch + k_LH(T,I_ch,SOC>82%) · Q_ch          (SCM18 Eq. 23)
```

Rate form: dQ_loss/dt = k_cal²/(2 q_cal) + k_H² |I|/(2 q_H) + k_L² I_ch/(2 q_L) + k_LH I_ch, with one sub-state per mechanism.

* U_a(SOC) is the graphite potential. It is what makes calendar aging step up above ≈ 50–60 % SOC (KEIL16, G9/G13).
* With LFP parameters (G14–G17): calendar loss ≈ 3.5 % after 1 yr at 25 °C, 50 % SOC (DERIVED). The cold-charging term is ×3.3 at 10 °C and ×7.8 at 0 °C (G18). This explains the JEITA reduced-current zone (H8).

**(c) Recommended lumped law for our paper.** This is (b) with phone-calibrated constants:

```
dq_cal/dt = k_c² A(T; E_cal) g(SOC)² / (2 q_cal)
dq_cyc/dt = |I|/C_nom · k_y² [A(T; E_H) + A⁻¹(T; E_L) e^(β (I_ch/C_nom − c_ref))]² / (2 q_cyc)
r_growth  = γ · q_loss,   γ ≈ 2–3, UQ 2–7   (R +40–60 % at 20 % fade per SCH14 G8; NCA storage in KEIL16 G9/G11 gives γ ≈ 6.6)
A(T; E)   = exp[−(E/R_g)(1/T − 1/T_ref)],   T_ref = 298.15 K
```

| Constant | Value | Basis |
|---|---|---|
| E_cal | 58 kJ/mol (UQ: 20–60) | SCH14 / SCM18 |
| E_H | 33 kJ/mol | SCM18 |
| E_L | 55 kJ/mol | SCM18 |
| g(SOC) | piecewise: ≈ 0.7 (SOC < 20 %), 1.0 (20–60 %), ≈ 1.6–1.9 (> 60–70 %) | KEIL16 NCA at 25 °C, normalised to 50 % SOC; read-off, M |
| k_c | Calibrated so q_cal(1 yr, 25 °C, 50 % SOC) ≈ 2.5–3 % | SCH14, KEIL16 |
| k_y | Calibrated so that q_cal + q_cyc = 20 % at N_80 = 500 or 1000 full cycles at 25 °C and 1 cycle/day (A1, A2) | — |

Calibration example (DERIVED): over 500 days at 25 °C and 3.85 V average, SCH14 gives q_cal ≈ 3.9 %, so the cycle part is ≈ 16 % at 500 FCE and k_y ≈ 0.16/√500 ≈ 0.0072 per √FCE. For the 1000-cycle rating, q_cal(1000 d) ≈ 6.5 %, so the cycle part is 13.5 % and k_y ≈ 0.0043.

**(d) Mechanistic option (RAM04).** Tafel side reaction during charge, with film thickness δ and resistance R_film = R_SEI + δ/κ_P as states (G22). It is physically richest but needs a P2D-level electrode potential. Use it only to justify the √t / Arrhenius / high-SOC structure of (b)–(c).

### 2.5 Charging-history stress (what OBC and 80 % limits change)

* **Calendar stress is time-weighted over SOC.** For overnight charging, compare 8 h at 100 % with 8 h at 80 % plus 1 h at 100 % (OBC). With SCH14's voltage law at 25 °C, 4.1 V vs 3.85 V gives α_cap ratio 1.36, and 4.2 vs 3.85 V gives 1.50 (DERIVED). With KEIL16's NCA step, 100 % vs 50 % storage is ≈ 1.8× (G9). Effect on yearly fade: of order +0.5–1.5 %-points per year at 25 °C. Larger when warm: the same dwell at 35 °C costs ≈ 2.1× more (G5).
* **Heat while charging matters more than charge power per se.** KANG19 measured ≈ 37 °C while fast-charging an idle phone, but 72 °C when charging while video-chatting at low SOC (H10). Through A(T), a 10 K rise doubles calendar-type aging.
* **Cold charging.** Use the k_L term (inverse Arrhenius and charge-current exponential). JEITA halves current at 0–10 °C and forbids charging below 0 °C (H8), which in the model means I_ch → 0.5 I_ch or 0.
* **Charge voltage.** Lower end-of-charge voltage gives less SEI side reaction (RAM04 simulation: 7.2 → 3.8 % for 4.2 → 3.9 V). Phone LCO cells charged to 4.35–4.4 V age much faster when held at top of charge at 40 °C (NEL18). This is the physical basis of 80 % limits and Pixel Battery Health Assistance voltage steps.

### 2.6 Self-discharge term

dSOC/dt ⊃ −k_sd(T, SOC), with k_sd ≈ 0.1 %/month at 25 °C (S2). Use an Arrhenius multiplier with the calendar E (assumption) for UQ up to 1–2 %/month (S3). Effect on TTE: < 0.1 % of a day's budget (S6). It can be dropped from TTE runs and kept only for multi-week storage scenarios.

### 2.7 Coupling to TTE, and scaling to other devices

* **Aged battery.** Q_max = Q₀(1 − q_loss) scales TTE roughly linearly. R_int = R₀(1 + r_growth) adds an earlier V_cut crossing. The cut-off penalty is strongly nonlinear in the cold and at low SOC (T3). This is the mechanism behind Apple's 2016–2017 shutdowns and "performance management" (C2, C17–C19).
* **Other devices.** C_th ∝ mass; R_th ∝ 1/(h·A); τ ∝ m/A, i.e. ∝ thickness × density × c_p. Tablets and laptops (larger A, higher P_sust) and wearables (tiny C_th, τ of order a minute) follow from the same equations. This scaling is DERIVED from physics; no device data were collected for it here.

### 2.8 Arithmetic used for DERIVED rows (re-runnable)

```python
import math; R=8.314
acap=lambda V,T:(7.543*V-23.75)*1e6*math.exp(-6976/T)     # SCH14 via CA21, t in days
ares=lambda V,T:(5.270*V-16.32)*1e5*math.exp(-5986/T)
bcap=lambda V,D:7.348e-3*(V-3.667)**2+7.600e-4+4.081e-3*D
bres=lambda V,D:2.153e-4*(V-3.725)**2-1.521e-5+2.798e-4*D
acap(3.7,298.15)*365**0.75            # 0.0239  -> 2.4 %/yr
math.exp(6976*(1/298.15-1/308.15))    # 2.14   per +10 K
(0.2/bcap(3.7,1.0))**2/2.05           # 830 FCE to 80 % (cycle term only)
R*math.log(2.8)/(1/253.15-1/298.15)/1e3   # 14.4 kJ/mol (LEE20 R0 ratio -20/25 C)
3.4*(3.7-3.4)/0.6                     # 1.7 W  (T3, cold)
-5/math.log(9/17)                     # 7.9 min (HE19 Fig. 11)
4.2e-6*730/3.5*100                    # 0.088 %/month (ZIL19)
```

---

## 3. Key findings

1. **The phone is a ~150 J/K, ~3 K/W body with τ ≈ 5–15 min.** Mass × specific heat, HE19's measured exponential settling (≈ 7.5 min) and KANG19's energy-vs-temperature curve all give C_th ≈ 100–230 J/K. Skin hot spots run hotter than the bulk (R_skin ≈ 5–7 K/W vs ≈ 3 K/W). A one-node ODE with a hot-spot offset is adequate for TTE.
2. **Throttling caps sustained power near 3 W, and the cap depends on ambient.** With a 45 °C skin limit, P_sust ≈ 2.2–3.8 W (CHI15) at room temperature. Video chat or recording pushes hot spots past the 42–45 °C pain band within 6–8 min (KANG19). P_sust = (T_lim − T_amb)/R_skin implies roughly half the sustainable power at 35 °C ambient.
3. **Cold acts through resistance, not lost charge.** R_int rises 2–3.3× from 25 °C to −5 °C (≈ 15–20 kJ/mol effective; steeper below 5 °C). Phones shut off at V_cut ≈ 3.4 V while showing 15–35 % SOC (−5 °C), and 30 % or more at −15 °C. A cold, warmed phone gives back 330 mAh. The operating-time penalty at 0 °C was 15 % for a UI workload (LEE20). Apple confirms the effect is temporary.
4. **Aging, cold and low SOC compound.** Aging raises R_int (and cuts peak power) as well as cutting capacity. The same burst load that is safe on a new warm battery can trip V_cut on an aged cold one (T3). This is exactly Apple's 2017 explanation and the reason for iOS performance management.
5. **Rated life: 20 % capacity loss in 500–1000 full cycles.** iPhone ≤ 14: 500 cycles; iPhone 15+: 1000; Pixel 3–8 Pro: 800; Pixel 8a+: 1000. At 1 cycle/day that is ≈ 1.4–2.7 years.
6. **Heat is the dominant controllable aging stress.** Calendar aging grows ≈ ×2 per +10 K (E ≈ 58 kJ/mol, NMC). Storage at 50 °C for 10 months cost NCA/NMC cells 12–23 % at full charge vs 5.5–6.5 % at 25 °C (KEIL16). Apple and Google both warn that ambient > 35 °C permanently damages capacity.
7. **High-SOC dwell matters, but less than heat.** Calendar fade is a step function of SOC tied to graphite potential (step near 60 % SOC), not a smooth ramp. Keeping a phone at 100 % rather than ≈ 50–80 % costs ≈ 1.3–1.9× calendar fade. This is the physical rationale for Optimized/Adaptive charging and 80 % limits.
8. **There is a temperature sweet spot near 25 °C for cycling.** Below it, plating-type damage grows (inverse Arrhenius, worse with charge current and at high SOC). Above it, SEI growth grows (WAL14, WU17, SCM18). Fast charging in the cold is the worst combination; JEITA zones encode this.
9. **Surprisingly small effects.**
   * Self-discharge (≈ 0.1 %/month cell-level; ≤ 2 %/month with protection electronics) is two to three orders of magnitude below standby drain.
   * The battery's own I²R heat (~0.1–0.3 W) is minor next to device dissipation (2–6 W).
   * Discharge rate matters little for pouch-cell aging (1C vs 5C: 4800 vs 3500 EFC; WU17). For LCO, reducing DOD did not help (CHOI02).
   * Charger wattage alone: informal tests show < 1 pp extra wear over 500 cycles.
10. **Users keep batteries high and warm.** Half of charge events end above 90 % SOC (GUO17). Charging while using the phone at low SOC produced the hottest temperatures measured, up to 72 °C (KANG19). Recommendations should target heat during charge and long full-charge dwell, not "never fast-charge".

---

## 4. Caveats

* **Schmalstieg 2014 equations are second-hand.** They come from a restatement (CA21), corroborated by two other papers' coefficient tables seen in search snippets. The original paper was not accessible here. Units still to confirm: whether ΔDOD is a fraction or a percentage (we use a fraction, which gives plausible results), whether Q counts charge plus discharge Ah (a factor-2 ambiguity, i.e. √2 on β_cap), and the voltage range of validity. Phone cells at 4.35–4.45 V lie outside the tested window. One restatement writes the exponent as exp(1 − 6976/T), which differs; we use exp(−6976/T).
* **Chemistry mismatch.** No open, quantitative aging model for a modern phone LCO/graphite (4.4 V) cell was found. SCH14 is NMC, SCM18 and WANG11 are LFP, KEIL16 is NCA/NMC/LFP, WU17 is an LCO/NCA blend. The SCM18 activation energies (LFP) and the SCH14 E_cal (NMC) differ by 3×, so we carry E_cal = 20–60 kJ/mol in UQ.
* **Figure read-offs (M).** KEIL16 Fig. 2, LEE20 Figs. 7/10/14/15, HE19 Figs. 10/11 and KANG19 Figs. 6/7 carry ±10–20 % reading error. Phone-level "battery resistance" in LEE20/HE19 includes connector, protection-circuit and fuel-gauge path resistance and RC polarization at the measurement time scale. It is not pure cell R0.
* **Thermal parameters are device- and grip-specific.** R_th, C_th and τ come from Nexus 5X/Pixel 5-era phones and a few designs. Cases, hands, pockets and desks change R_th. Phone masses and the specific heats of glass, steel and plastic were not verified in this session.
* **Manufacturer pages are living documents.** Apple and Google pages carry 2025–2026 dates and were read on 2026-11-17; content may change. The 2017 Apple letter's original URL now returns 404, so its quotes are via a blog archive. The 2016–2017 shutdown statistics come via press (M).
* **Self-discharge sources disagree.** The only peer-reviewed number (ZIL19) is a reversible self-discharge current on a non-phone cell. Battery University's table conflicts with its own 1–2 %/month text.
* **Fast-charging evidence for phones is informal.** The controlled evidence concerns cells (plating at low T / high current), not phone charging systems with thermal management.
* **Access limits.** Several primary PDFs (Ramadass 2002, Wang 2011, Choi & Lim 2002, Waldmann 2014, Zimmerman 2004) were behind paywalls or returned errors. Their numbers are taken from abstracts, restatements or search snippets and graded M.
* **Data licences.** CALCE (LCO CS2/CX2 cycling data) is "open access" with a citation request but no explicit licence. NASA PCoE 18650 aging data: "License not specified", access level public. Both are usable for validation if the paper documents this. The CC BY papers (KEIL16, SCM18, WU17, NEL18, SS18, TOM19) are the safest sources for reproduced figures and data.

---

## 5. References

Access and licence notes describe what we observed. "©" means the publisher holds copyright and only cited numbers are used. Web pages were accessed 2026-11-17.

1. **[AP-TEMP]** Apple Support. "If your iPhone or iPad gets too hot or too cold." Article 118431 (formerly HT201678), page date 14 May 2026. https://support.apple.com/en-us/118431 (free to read, © Apple).
2. **[AP-BATT]** Apple Support. "iPhone battery and performance." Article 101575, page date 1 June 2026. https://support.apple.com/en-us/101575 (free, © Apple).
3. **[AP-OBC]** Apple Support. "About Charge Limit and Optimized Battery Charging on iPhone." Article 108055 (formerly HT210512), page date 8 Dec 2025. https://support.apple.com/108055 (free, © Apple).
4. **[AP-MAX]** Apple. "Batteries – Maximizing Performance." https://www.apple.com/batteries/maximizing-performance/ (free, © Apple).
5. **[AP-WHY]** Apple. "Batteries – Why Lithium-ion?" https://www.apple.com/batteries/why-lithium-ion/ (free, © Apple).
6. **[AP-FAST]** Apple Support. "Fast charge your iPhone." Article 102574, page date 22 Oct 2026. https://support.apple.com/en-us/102574 (free, © Apple).
7. **[AP-2017]** Apple. "A Message to Our Customers about iPhone Batteries and Performance." 28 Dec 2017. The original apple.com newsroom URL returned HTTP 404 in this session. Quotations taken from the archive in M. Tsai's blog, https://mjtsai.com/blog/2017/12/28/apple-iphone-battery-message/ (secondary).
8. **[PRESS-17]** M. Tsai blog, 24 Feb 2017, quoting MacRumors/TechCrunch on Apple's iOS 10.2.1 statement (80 % / 70 % fewer shutdowns). https://mjtsai.com/blog/2017/02/24/ios-10-2-1-update-reduces-unexpected-shutdowns/ (secondary/press).
9. **[PRESS-16]** iPhone in Canada, 2 Dec 2016, "Apple Explains iPhone 6s Battery Shutdown Problem to Chinese Customers" (quoting Apple's China statement via Business Insider). https://www.iphoneincanada.ca/2016/12/02/apple-iphone-6s-battery-shutdown-problem/ (secondary/press).
10. **[GOOG]** Google Pixel Phone Help. "Understand your Pixel battery." https://support.google.com/pixelphone/answer/15738128 (free, © Google). Stage limits for Battery Health Assistance (200–1000 cycles) are from 9to5Google press coverage, https://9to5google.com/2025/03/31/pixel-battery-health-assistance/ (secondary, not opened).
11. **[SAM]** Samsung Singapore Support. "Galaxy Battery protection feature in One UI 6.1." Last updated 11 Apr 2024. https://www.samsung.com/sg/support/mobile-devices/galaxy-battery-protection-feature-in-one-ui-6-1/ . Legacy 85 % "Protect battery": https://www.samsung.com/in/support/mobile-devices/battery-protection-feature-in-samsung-galaxy-phones (free, © Samsung). No official Samsung operating-temperature range was found.
12. **[AOSP-TH]** Android Open Source Project. "Thermal mitigation." https://source.android.com/docs/core/power/thermal-mitigation (AOSP documentation; licence not checked).
13. **[TI-JEITA]** J. Qian. "Li-ion battery-charger solutions for JEITA compliance." *Texas Instruments Analog Applications Journal*, 1Q 2010, literature no. SLYT365. https://www.ti.com/lit/pdf/slyt365 (free, © TI). JEITA/BAJ guideline released 20 Apr 2007, as stated there.
14. **[UL-TS]** F. Chiang (UL Taiwan). "IEC 62368-1 Technical Brief: Touch Temperature Limits." 22 Oct 2010. https://japan.ul.com/wp-content/uploads/sites/27/2014/06/1_techbrief_touchtemp.pdf (free, © UL; summarises IEC 62368-1 Table 42 / ISO 13732-1:2006).
15. **[KANG19]** S. Kang, H. Choi, S. Park, C. Park, J. Lee, U. Lee, S.-J. Lee. "Fire in Your Hands: Understanding Thermal Behavior of Smartphones." *ACM MobiCom 2019*, Los Cabos. DOI: 10.1145/3300061.3300128. Author PDF: https://ic.kaist.ac.kr/publications/papers/kang2019fire.pdf (free author copy, © owner/authors, publication rights ACM).
16. **[PU23]** H. Pu, X. Wu. "Investigating Skin Temperature-Based Overheating in mmWave Smartphones: Power and Thermal Models for Optimal Non-Throttling Performance." arXiv:2304.06512 (2023). https://arxiv.org/abs/2304.06512 (arXiv preprint).
17. **[CHI15]** V. Chiriac, S. Molloy, J. Anderson, K. Goodson. "A Figure of Merit for Smart Phone Thermal Management." *Electronics Cooling*, 30 Nov 2015. https://www.electronics-cooling.com/?p=19849 . Reprinted under V. Chiriac, 19 Sept 2023 (= CHI23 in device_power.md): https://www.electronics-cooling.com/?p=35382 (trade magazine, free).
18. **[HE19]** L. He, Y. Lee, E. Kim, K. G. Shin. "Environment-Aware Estimation of Battery State-of-Charge for Mobile Devices." *ACM/IEEE ICCPS 2019*, Montreal. Author PDF: https://rtcl.eecs.umich.edu/rtclweb/assets/publications/2019/ICCPS-liang19.pdf (free author copy; the DOI field in the PDF is a placeholder, so no DOI given).
19. **[LEE20]** Y. Lee, L. He, K. G. Shin. "Causes and Fixes of Unexpected Phone Shutoffs." *ACM MobiSys 2020*, Toronto. DOI: 10.1145/3386901.3389024. Copies: https://rtcl.eecs.umich.edu/rtclweb/assets/publications/2020/youngmoon-mobisys20.pdf and NSF PAR https://par.nsf.gov/servlets/purl/10301364 (free public-access copy, © ACM).
20. **[MAL99]** H. Maleki, S. A. Hallaj (as listed by OSTI), J. R. Selman, R. B. Dinwiddie, H. Wang. "Thermal Properties of Lithium-Ion Battery and Components." *J. Electrochem. Soc.* 146(3), 947–954 (1999). DOI: 10.1149/1.1391704. OSTI record: https://www.osti.gov/biblio/345380 (abstract; © ECS).
21. **[SS18]** J. Schmalstieg, D. U. Sauer. "Full Cell Parameterization of a High-Power Lithium-Ion Battery for a Physico-Chemical Model: Part II. Thermal Parameters and Validation." *J. Electrochem. Soc.* 165(16), A3811 (2018). DOI: 10.1149/2.0331816jes. **CC BY 4.0.**
22. **[BER85]** D. Bernardi, E. Pawlikowski, J. Newman. "A General Energy Balance for Battery Systems." *J. Electrochem. Soc.* 132(1), 5–12 (1985). DOI: 10.1149/1.2113792 (© ECS; equation form only).
23. **[SCH14]** J. Schmalstieg, S. Käbitz, M. Ecker, D. U. Sauer. "A holistic aging model for Li(NiMnCo)O2 based 18650 lithium-ion batteries." *J. Power Sources* 257, 325–334 (2014). DOI: 10.1016/j.jpowsour.2014.02.012 (© Elsevier; original not accessed, equations via CA21).
24. **[CA21]** M. Clarke, J. J. Alonso. "Lithium–Ion Battery Modeling for Aerospace Applications." *Journal of Aircraft* (2021, Article in Advance). DOI: 10.2514/1.C036209. PDF: https://suave.stanford.edu/publications/2021_Battery_Degradation_Modeling.pdf (free author copy, © AIAA). Restates SCH14 Eqs. (14a–d), (15a–b).
25. **[SCM18]** M. Schimpe, M. E. von Kuepach, M. Naumann, H. C. Hesse, K. Smith, A. Jossen. "Comprehensive Modeling of Temperature-Dependent Degradation Mechanisms in Lithium Iron Phosphate Batteries." *J. Electrochem. Soc.* 165(2), A181–A193 (2018). DOI: 10.1149/2.1181714jes. **CC BY 4.0.**
26. **[KEIL16]** P. Keil, S. F. Schuster, J. Wilhelm, J. Travi, A. Hauser, R. C. Karl, A. Jossen. "Calendar Aging of Lithium-Ion Batteries: I. Impact of the Graphite Anode on Capacity Fade." *J. Electrochem. Soc.* 163(9), A1872–A1880 (2016). DOI: 10.1149/2.0411609jes. **CC BY 4.0.** PDF: https://mediatum.ub.tum.de/doc/1651485/document.pdf
27. **[WU17]** Y. Wu, P. Keil, S. F. Schuster, A. Jossen. "Impact of Temperature and Discharge Rate on the Aging of a LiCoO2/LiNi0.8Co0.15Al0.05O2 Lithium-Ion Pouch Cell." *J. Electrochem. Soc.* 164(7), A1438 (2017). DOI: 10.1149/2.0401707jes. **CC BY 4.0.**
28. **[NEL18]** K. J. Nelson, J. E. Harlow, J. R. Dahn. "A Comparison of NMC/Graphite Pouch Cells and Commercially Available LiCoO2/Graphite Pouch Cells Tested to High Potential." *J. Electrochem. Soc.* 165(3), A456 (2018). DOI: 10.1149/2.0041803jes. **CC BY 4.0.**
29. **[RAM02]** P. Ramadass, B. Haran, R. E. White, B. N. Popov. "Capacity fade of Sony 18650 cells cycled at elevated temperatures: Part I. Cycling performance." *J. Power Sources* 112(2), 606–613 (2002). DOI: 10.1016/S0378-7753(02)00474-3. Record: https://scholarcommons.sc.edu/eche_facpub/443 (© Elsevier). Fade percentages via a search snippet of the abstract. Cell and protocol confirmed in B. N. Popov, P. Ramadass, B. S. Haran, R. E. White, "Capacity Fade Studies of Li-ion Cells Cycled at Different Temperatures," ECS meeting abstract (2002; meeting 201 per the URL path), https://www.electrochem.org/dl/ma/201/pdfs/0134.pdf
30. **[RAM04]** P. Ramadass, B. Haran, P. M. Gomadam, R. White, B. N. Popov. "Development of First Principles Capacity Fade Model for Li-Ion Cells." *J. Electrochem. Soc.* 151(2), A196–A203 (2004). DOI: 10.1149/1.1634273. USC Scholar Commons copy (© ECS, redistribution per ECS terms).
31. **[WANG11]** J. Wang, P. Liu, J. Hicks-Garner, E. Sherman, S. Soukiazian, M. Verbrugge, H. Tataria, J. Musser, P. Finamore. "Cycle-life model for graphite-LiFePO4 cells." *J. Power Sources* 196(8), 3942–3948 (2011). DOI: 10.1016/j.jpowsour.2010.11.134 (© Elsevier; not accessed). Constants via FAR24. Author list as given in FAR24's reference list.
32. **[FAR24]** A. Farakhor, D. Wu, P. Chen, J. Wang, Y. Wang, H. Fang. "Economic Optimal Power Management of Second-Life Battery Energy Storage Systems." arXiv:2408.16197 (2024), Table I. https://arxiv.org/abs/2408.16197 (arXiv).
33. **[WAL14]** T. Waldmann et al. "Temperature dependent ageing mechanisms in Lithium-ion batteries – A Post-Mortem study." *J. Power Sources* 262, 129–135 (2014). ADS: https://ui.adsabs.harvard.edu/abs/2014JPS...262..129W (abstract only; © Elsevier).
34. **[CHOI02]** S. S. Choi, H. S. Lim. "Factors that affect cycle-life and possible degradation mechanisms of a Li-ion cell based on LiCoO2." *J. Power Sources* 111(1), 130–136 (2002) (abstract via search; © Elsevier).
35. **[YM12]** T. Yuksel, J. Michalek. "Development of a Simulation Model to Analyze the Effect of Thermal Management on Battery Life." SAE Technical Paper 2012-01-0671 (2012). DOI: 10.4271/2012-01-0671. Author copy: https://www.cmu.edu/me/ddl/publications/2012-SAE-Yuksel-Michalek-Thermal-Mgmt.pdf (© SAE).
36. **[ZIL19]** I. Zilberman, S. Ludwig, A. Jossen. "Cell-to-cell variation of calendar aging and reversible self-discharge in 18650 nickel-rich, silicon–graphite lithium-ion cells." *J. Energy Storage* 26, 100900 (2019). DOI: 10.1016/j.est.2019.100900. Abstract via TUM portal: https://portal.fis.tum.de/en/publications/cell-to-cell-variation-of-calendar-aging-and-reversible-self-disc/ (access licence not confirmed).
37. **[BU802b]** Battery University (Cadex). "BU-802b: What does Elevated Self-discharge Do?" Updated 2 Nov 2021. https://batteryuniversity.com/article/elevating-self-discharge (industry, not peer-reviewed).
38. **[HOQ16]** M. A. Hoque, M. Siekkinen, J. Koo, S. Tarkoma. "Accurate Online Full Charge Capacity Modeling of Smartphone Batteries." arXiv:1604.05689 (2016). https://arxiv.org/abs/1604.05689 (arXiv).
39. **[GUO17]** Y. Guo, C. Wang, X. Chen. "Understanding Application-Battery Interactions on Smartphones: A Large-Scale Empirical Study." *IEEE Access* 5, 13387–13400 (2017). DOI: 10.1109/ACCESS.2017.2728620 (open-access journal; exact licence not checked). PDF: https://yaoguopku.github.io/papers/Guo-Access-17.pdf
40. **[TOM19]** A. Tomaszewska et al. (20 authors, incl. G. Offer and B. Wu). "Lithium-ion battery fast charging: A review." *eTransportation* 1, 100011 (2019). DOI: 10.1016/j.etran.2019.100011. **CC BY** (per search summary). Used at abstract level only.
41. **[CALCE]** CALCE Battery Research Group, University of Maryland. Battery Data (CS2 1100 mAh and CX2 1350 mAh LiCoO2 prismatic cells; CX2_4 cycled at 25/35/45/55 °C). https://calce.umd.edu/battery-data ("open access", cite He et al. 2011, Xing et al. 2013, Williard et al. 2013; no explicit licence).
42. **[NASA-PCoE]** B. Saha, K. Goebel. "Battery Data Set." NASA Prognostics Data Repository, NASA Ames (2007). Catalog: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets (access level public; "License not specified"). 18650 cells to 30 % fade; charge, discharge and EIS at different temperatures.
43. **[BH19]** G. Bhat, S. Gumussoy, U. Y. Ogras. "Power and Thermal Analysis of Commercial Mobile Platforms: Experiments and Case Studies." *DATE 2019*; arXiv:1904.09814. Cross-referenced from device_power.md; not re-read here.

### Source-key → reference number map

AP-TEMP = 1, AP-BATT = 2, AP-OBC = 3, AP-MAX = 4, AP-WHY = 5, AP-FAST = 6, AP-2017 = 7, PRESS-17 = 8, PRESS-16 = 9, GOOG = 10, SAM = 11, AOSP-TH = 12, TI-JEITA = 13, UL-TS = 14, KANG19 = 15, PU23 = 16, CHI15 = 17, HE19 = 18, LEE20 = 19, MAL99 = 20, SS18 = 21, BER85 = 22, SCH14 = 23, CA21 = 24, SCM18 = 25, KEIL16 = 26, WU17 = 27, NEL18 = 28, RAM02 = 29, RAM04 = 30, WANG11 = 31, FAR24 = 32, WAL14 = 33, CHOI02 = 34, YM12 = 35, ZIL19 = 36, BU802b = 37, HOQ16 = 38, GUO17 = 39, TOM19 = 40, CALCE = 41, NASA-PCoE = 42, BH19 = 43.

*Downloaded working copies (untrusted, not committed) are in /tmp/claude-0/research_dl/thermal_aging/.*
