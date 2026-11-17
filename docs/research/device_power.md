# Smartphone subsystem power consumption: research note

Scope: parameters and model forms for the component power model used by our continuous-time battery model (2026 MCM Problem A).

```
P_total(t) = P_base + P_screen(b, A, content, f_refresh) + P_SoC(u, f) + P_net(mode, throughput, signal, tail)
             + P_GNSS + P_cam + P_audio + P_sensors + P_bg
I_batt(t)  = P_total(t) / V_term(t)        (constant-power load seen at the battery terminals)
```

Conventions used below:

* **Where power is measured.** Unless a row says otherwise, powers are **at the battery terminals**, so they already include PMIC/DC-DC conversion loss. Component ("rail") powers from instrumented boards are marked "rail".
* **Converting mA to mW.** Values reported as current in mA are converted with V = 3.7 V (Ding et al. use 3.7 V) or 3.8 V. Use 3.7–3.85 V for nominal Li-ion. The nominal voltage itself was **not verified** in this session.
* **Confidence levels.** H means read directly from the text or tables of a peer-reviewed or primary source. M means read off a figure (±10–20 %), taken from a secondary report of a primary source, or measured on an old device. L means derived or extrapolated by us, or an unverified background value.
* **DERIVED rows.** Rows tagged **DERIVED** are our own arithmetic on cited numbers. The arithmetic is shown so it can be re-checked.
* **Source keys.** Keys such as `CH13` resolve to the numbered reference list in §5.

---

## 0. Executive summary: reference parameter set

This set is for a 2021–2024 phone: 6.1–6.7" OLED screen, ~5000 mAh battery. It is a **DERIVED** synthesis of the detailed table in §1, intended as model defaults with ranges for uncertainty quantification (UQ).

| Quantity | Symbol | Default | Range for UQ | Unit | Basis (see §1) | Conf. |
|---|---|---|---|---|---|---|
| Deep-sleep floor (screen off, suspended, radio idle/paging) | P_sleep | 30 | 15–60 | mW | CH13 24–27 mW suspend; HU12 11.4 mW base + ~20 mW LTE paging; CH15 7–15 mA | M |
| Extra for 5G attachment while idle | ΔP_5G,idle | ×1.3–2 on cellular idle part | 0–150 | mW | GU24 (5G ≈ 2× 4G cellular standby on Pixel 7a) | M |
| "Awake but idle" (screen off, wakelock held) | ΔP_awake | 100 | 50–250 | mW | CH10 (Freerunner 269 vs 69 mW); AOSP example cpu.awake 50 mA (illustrative) | L |
| Display fixed overhead when on (driver, touch, panel bias) | P_d0 | 0.15 | 0.1–0.25 | W | CH13 S3: 223 mW constant term; AOSP screen.on example 200 mA | M |
| Display at 50 % brightness, typical UI, 6.1–6.7" OLED | P_disp(0.5) | 0.4 | 0.25–0.9 | W | GU24 Pixel 7a ≈0.30 W, S23+ ≈0.45 W; CH13 S3 0.3–0.6 W at 4.8" (area-scaled upper bound 0.43–1.04 W) | M |
| Display, full-white screen at max (non-HBM) brightness, 6.1–6.7" | P_disp,max | 1.5 | 1.0–2.0 | W | CH13 1124 mW at 4.8"; CH10 N1 white max +1313 mW; area scaling (DERIVED) | L |
| Brightness law (OLED, content-dependent part) | g(b) | affine: g(min)=0.16, g(0.5)=1, g(max)=2.0 | — | × | CH13 | H |
| 120 Hz refresh penalty (fixed-rate panel) | ΔP_RR | 0.17 | 0.10–0.20 (0.03–0.10 with LTPO, unverified) | W | AT20 via secondary: ~160–180 mW on S20 Ultra | M/L |
| SoC, light interactive (scrolling, UI) | P_SoC,light | 0.4 | 0.2–0.8 | W | GU24 compute+memory ≈0.4–0.5 W in Instagram | M |
| SoC, sustained heavy (3D game), thermally capped | P_SoC,sust | 3.0 | 2.2–3.8 | W | CHI23 sustainable dissipation at 45 °C skin | M |
| SoC burst peak | P_SoC,peak | 5 | 3–7+ | W | CH13 core 2.85 W + GPU 1.42 W (S3); device >7 W | M |
| Wi-Fi idle (associated, PSM beacons) | P_wifi,idle | 5 | 2–50 | mW | HU12 beacon avg ≈1.6 mW; CH15 1.1 mA; CH13 says +50 mW (awake) | M |
| Wi-Fi active, streaming-type load | P_wifi,act | 0.3 | 0.13–0.8 | W | HU12 β=133 mW + 137 mW/Mbps DL; SU14 240–795 mW; GU24 ≤0.66 W (upper bound) | M |
| LTE idle paging average | P_LTE,idle | 20 | 5–35 | mW | HU12 DRX (DERIVED 19.7 mW); CH15 2.3–8.3 mA | M |
| LTE active, moderate DL | P_LTE,act | 0.9 | 0.6–1.5 | W | GU24 Pixel 7a cellular 0.73–0.95 W (YouTube/Instagram) | M |
| LTE tail power / timer (modern) | P_tail, T_tail | 0.15 W / 10 s | 0.07–1.06 W / 5–11.6 s | W, s | NA21 66–178 mW, 5–10.2 s; HU12 (2012) 1.06 W, 11.6 s | H |
| 5G NSA low-band tail / timer | P_tail,5G | 0.25 W / 10.3 s | 0.25–0.6 W / 10.2–21.4 s | W, s | NA21; XU20 (20 s chained tails) | H |
| 5G mmWave tail | P_tail,mmW | 1.09 W / 10.5 s | — | W, s | NA21 | H |
| Weak-signal multiplier on cellular active energy (−85 → −105 dBm) | m_sig | 1.5 | 1.05–3 | × | DI13 (+52 % energy; Tx/Rx current ×1.2–3.1) | M |
| GNSS on (1 Hz fixes) | P_GNSS | 0.3 | 0.23–0.45 | W | KW19 232–318 mJ per fix; CH13 386–433 mW | H |
| Camera pipeline (preview/record), whole device | P_cam,dev | 2.0 | 1.4–2.6 | W | CH13; CX13 | M |
| Audio playback (screen off, headphones) | P_audio | 0.05 above sleep (total ≈0.23 W, S3) | 0.03–0.2 | W | CH13; HE20 0.11 W | M |
| Sensors at 50 Hz (accel / gyro / mag) | P_sens | 5 / 30 / 12 | 1–30 | mW | CH13 | H |
| PMIC conversion efficiency (rail ← battery) | η_PMIC | 0.83 | 0.75–0.92 | — | CH13 | M |

---

## 1. Detailed parameter table

### 1.1 Display

| # | Name | Symbol | Value | Range | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|---|---|
| D1 | OLED display power model, Galaxy S III | P_disp = P_d0 + k·P_O | P_d0 = 223, k = 1.57 | — | mW, — | S3 4.8" 720×1280 Super AMOLED. P_O is the one OLED supply they could measure. <3 % error over several images and brightness levels | CH13 | H |
| D2 | OLED display at 50 % brightness | P_disp(0.5) | 300–600 | — | mW (rail ×1.2) | S3, all scenarios; "highly dependent on image content" | CH13 | H |
| D3 | Brightness multiplier on content-dependent OLED part | g(b) | 0.16 (min), 1 (50 %), 2.0 (max) | — | × | S3; apply to the OLED component only. Nearly affine in slider position | CH13 | H |
| D4 | Max display power | P_disp,max | 1124 | — | mW | S3, white screen, max brightness | CH13 | H |
| D5 | Colour-scheme saving | ΔP_color | up to 300 (medium), up to 600 (max) | — | mW | S3 | CH13 | H |
| D6 | Nexus One AMOLED (3.7", 480×800) | ΔP_white | +194 (min brightness), +1313 (max) | — | mW | Black screen power is fixed whatever the brightness. White screen adds the listed power | CH10 | H |
| D7 | OLED pixel power is linear in **linear** RGB | P = Σ_pixels (a·R + b·G + c·B) | R² > 0.95 | — | — | Nexus One, Galaxy S, Nokia N85. Colour coefficients differ per device | DZ11 | H |
| D8 | Gamma 2.2 makes power nonlinear in sRGB grey level | P ∝ (grey/255)^2.2 per channel | — | — | — | Samsung AMOLED. Per-unit-area efficiency "almost the same" across generations (Nexus S → S3). Pixel-level efficiency gain only 5.7–29.5 % over 2 generations | CX13 | H |
| D9 | Dark-mode saving vs light mode | s_dark | 3–9 % (30–50 % brightness); 39–47 % (100 %) | — | % of battery power | Pixel 2, Moto Z3, Pixel 4, Pixel 5; 6 Google apps, 60 s each | DH21 (press release) | M |
| D10 | Equivalence: light mode at 20 % brightness ≈ dark mode at 50 % | — | — | — | — | Pixel 5, Google News | DH21 (press) | M |
| D11 | Typical indoor auto-brightness | b_auto | 30–40 | — | % | Statement in DH21 press release | DH21 (press) | M |
| D12 | Fraction of white pixels in mobile web pages | APL proxy | 65 % white pixels | — | % | LiveLab traces, 25 users, 3 months | DZ11 | H |
| D13 | Colour transforms (Dark / Green / Inversion / Arbitrary) | ΔP_disp | −25 / −34 / −66 / −72 | — | % display power | Nexus One web pages; system power −41 % for browsing | DZ11 | H |
| D14 | OLED display share in video and games | P_disp | 74.6–374.2 (full-screen video); 237–363 (games); 50.3–446.4 (range across games) | — | mW | Samsung Nexus S, Galaxy S/S2/S3, Galaxy Nexus (2010–12) | CX13 | H |
| D15 | Display rail on modern phones | P_disp | Pixel 7a ≈ 300; S23+ ≈ 450 | ±15 % | mW | Pixel 7a (6.1") ODPM rails; S23+ (6.6") battery-counter subtraction. Brightness "same nominal" but value not stated. **Read off Fig. 2/4/5** | GU24 | M |
| D16 | Area-normalised OLED power (S3) | p_A | 17.7 (max white); 4.7–9.4 (50 %, typical content) | — | mW/cm² | **DERIVED**: D2/D4 ÷ 63.5 cm² (4.8", 16:9) | CH13 + calc | L |
| D17 | Area-scaled to modern sizes (upper bound, 2012 panel efficiency) | P_disp | 6.1" (91.3 cm²): 0.43–0.86 (50 %), 1.62 (max white). 6.7" (110.2 cm²): 0.52–1.04 (50 %), 1.95 (max) | — | W | **DERIVED**: D16 × area, 19.5:9. Over-estimates modern panels (compare D15) | calc | L |
| D18 | LCD backlight vs brightness level | P_BL(b) | 7.8 (min) → 75 (level 143/255) → 414 (max) | — | mW | Openmoko Freerunner 2.8" LCD. **Convex** in level | CH10 | H (old) |
| D19 | LCD panel (no backlight), content effect | P_LCD | 33.1 (white) – 74.2 (black) | — | mW | Freerunner. G1: content effect ≤ 17 mW; button and keyboard backlights 189 mW | CH10 | H (old) |
| D20 | LCD brightness model is quadratic, not linear | P = Fdisp·(c0 + κH + λH²) | max brightness contribution ≈ 0.79 W | — | W | Fairphone 2 (IPS TFT LCD), video streaming model with k = 11 features | HE20 | H |
| D21 | Screen at 100 % (LCD, 2011 LTE phone) | P_screen | 847.2 ± 2.7 | — | mW | HTC LTE phone (MSM8655, Android 2.2.1); Monsoon | HU12 | H |
| D22 | 120 Hz vs 60 Hz penalty | ΔP_RR | ~160–180 | — | mW | Galaxy S20 Ultra (6.9" 3200×1440 AMOLED), black static screen, measured by AnandTech (original page offline; figure from secondary reports) | AT20 | M/L |
| D23 | Web-browsing battery life at 60 → 120 Hz | t_web | 14.05 → 11.33 h (Snapdragon); 12.28 → 10.18 h (Exynos) | — | h | S20 Ultra, AnandTech web test (via Ubergizmo). 9to5Google quotes 13.04 → 10.02 h for the same source, so the reports **conflict** | AT20 | M |
| D24 | Implied average power in web test | P̄_web | 1.37 (60 Hz) → 1.70 (120 Hz) | 1.3–1.9 | W | **DERIVED**: 5000 mAh × 3.86 V (V not verified) ÷ D23 hours. ΔP ≈ 0.33 W includes extra GPU/CPU rendering, not only the panel | calc | L |

### 1.2 SoC: CPU, GPU, memory, thermal

| # | Name | Symbol | Value | Range | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|---|---|
| S1 | Linux Energy-Model dynamic power law | P_dyn = C·V²·f | C in µW/MHz/V², f in MHz, V in V | — | µW | Used by the kernel's Energy-Aware Scheduling. C is "calculated from power measurements or derived by analysis" | LX | H |
| S2 | Per-core dynamic-power coefficient | C | SDM845 (Snapdragon 845): little 154, big 442. SM8250 (Snapdragon 865): little 105, big 379, prime 444. MT8183: A53 84, A73 211. RK3399: A53 100, A72 436 | — | µW/MHz/V² | Mainline device-tree values (vendor-supplied estimates) | LX | H (as data) |
| S3 | DVFS curve, big core | (f, V) → P | MT8183 A73: 793 MHz @ 0.700 V → 82 mW; 1417 MHz @ 0.850 V → 216 mW; 1989 MHz @ 1.050 V → 463 mW per core | — | MHz, V, mW | **DERIVED** from DTS OPP table × S2. V(f) ≈ affine, slope 0.29 V/GHz. Effective P ∝ f^1.88 | LX + calc | M |
| S4 | DVFS curve, little core | — | MT8183 A53: 28 mW (793 MHz, 0.65 V) → 184 mW (1989 MHz, 1.05 V) per core | — | mW | **DERIVED** | LX + calc | M |
| S5 | Big-core top OPP example | — | RK3399 A72: 1.8 GHz @ 1.2 V → ≈1.13 W per core | — | W | **DERIVED** | LX + calc | M |
| S6 | Maximum component power (S3) | P_max | CPU cores 2845 (AnTuTu); GPU 1415 (AnTuTu 3D); RAM 208; display 1124; 3G 1137 | — | mW (rail) | "The S3 can consume in excess of 7 W" | CH13 | H |
| S7 | Gaming power (whole device) | P_game | Angry Birds 1516; Need for Speed (3D) 2425 | — | mW | S3, 50 % brightness | CH13 | H |
| S8 | Memory and storage | P_RAM, P_flash | RAM ≈ 6 % of total on average (>100 mW in heavy scenes; up to 20 % in suspend). Flash < 8 mW (peak 30 mW in video recording) | — | % / mW | S3 | CH13 | H |
| S9 | Idle application processor | P_AP,idle | "negligible"; static power "still mostly a non-issue in phones" | — | — | 2013 finding; may not hold for hot, modern nodes | CH13 | M |
| S10 | Hardware vs software video decode | ΔP | HQ 720p: 1270 (HW) vs 2329 (SW). LQ: 1084 vs 1571. Saving > 1 W (45 %) | — | mW | S3, MX Player | CH13 | H |
| S11 | Software-decode power per frame-rate (normalised model) | P_fv | max contribution 0.35 W at 60 fps (0.175 W at 30 fps); constant offset 0.90 W; bitrate term 0.38 W; Wi-Fi bitrate extra 0.21 W | — | W | Fairphone 2, streaming model (k = 11) | HE20 | H |
| S12 | Sustainable heat dissipation at 45 °C skin limit | P_sust | 2.2 → 3.5 W (larger spreader); ≤ 3.8 W with ideal spreading; 0.16 W per °C of skin temperature | — | W | Simulation and test of phone designs; trade-journal article | CHI23 | M |
| S13 | Throttling impact | — | performance loss up to 34 % in popular apps; package reaches 50 °C without throttling | — | % | Nexus 6P (Snapdragon 810) | BH19 | H |
| S14 | Software (fuel-gauge) power readings vs Monsoon | SW/HW | 80.9–94.9 % (always under-reads); idle screen off: 80.9 % at 1 Hz, 94.9 % at 10 Hz | — | % | Samsung 5G phones | NA21 | H |

### 1.3 Network: Wi-Fi, LTE, 5G, signal strength

**Huang et al. 2012 (HU12)**: HTC LTE phone (MSM8655, Android 2.2.1), Monsoon monitor. All power readings in their Table 3 include the screen-off base power of 11.4 mW.

| # | Name | Symbol | Value | Unit | Notes | Conf. |
|---|---|---|---|---|---|---|
| N1 | Screen-off base power | P_base | 11.4 ± 0.4 | mW | | H |
| N2 | LTE promotion (IDLE → CONNECTED) | P_pro, T_pro | 1210.7 ± 85.6 mW; 260.1 ± 15.8 ms | mW, ms | | H |
| N3 | LTE tail (RRC_CONNECTED inactivity) | P_tail, T_tail | 1060.0 ± 3.3 mW; 11 576 ± 26 ms | mW, ms | Timer resets on every packet | H |
| N4 | LTE C-DRX on-duration (connected) | — | 1680 mW; T_on 1 ms; short and long cycles 20 / 40 ms | mW, ms | | H |
| N5 | LTE idle DRX (paging) | — | 594.3 ± 8.7 mW for 43.2 ms every 1280.2 ms | mW, ms | **DERIVED** average increment (594.3 − 11.4)·43.2/1280.2 = 19.7 mW | H / L(derived) |
| N6 | Data-transfer power, linear model P = α_u·t_u + α_d·t_d + β (t in Mbps) | α_u, α_d, β | LTE: 438.39, 51.97, 1288.04 (α_u/α_d = 8.44). 3G: 868.98, 122.12, 817.88. Wi-Fi: 283.17, 137.01, 132.86 | mW/Mbps, mW | Validated on apps with −1.1 to −5.5 % error; LTE uplink up to 3300 mW | H |
| N7 | Promotion + tail energy for one packet | E_1pkt | LTE 12.76; 3G 7.38; Wi-Fi 0.04 | J | | H |
| N8 | 3G (UMTS) state powers | — | promotion 659.4 mW / 582 ms; DCH tail 803.9 mW / 8088 ms; FACH tail 601.3 mW / 824 ms; idle DRX 374.2 mW / 55.4 ms every 5112 ms | | | H |
| N9 | Wi-Fi state powers | — | promotion 124.4 mW / 79.1 ms; tail 119.3 mW / 238.1 ms; beacon 77.2 mW / 7.6 ms every 308.2 ms (**DERIVED** average ≈ 1.6 mW) | | | H |
| N10 | Tail-energy share | — | LTE and 3G energy efficiency is dominated by tail energy, not data transfer. LTE is up to 23× less power-efficient than Wi-Fi on user traces | — | | H |

**Narayanan et al. 2021 (NA21)**: Samsung Galaxy S20 Ultra (S20U), S10 and Pixel 5. Verizon and T-Mobile, US. Monsoon at 5 kHz.

| # | Name | Value | Unit | Notes | Conf. |
|---|---|---|---|---|---|
| N11 | Tail power (average over tail incl. DRX) | Verizon 4G 178; T-Mobile 4G 66; Verizon NSA low-band (DSS) 249; Verizon NSA mmWave 1092; T-Mobile NSA low-band 260; T-Mobile SA low-band 593 | mW | Table 2 | H |
| N12 | 4G → 5G switch power | Verizon NSA-LB 799; Verizon NSA-mmWave 1494; T-Mobile NSA-LB 699; T-Mobile SA 245 | mW | Table 2 | H |
| N13 | RRC timers (UE inactivity / long DRX / idle DRX / 4G promo / 5G promo) | T-Mobile SA: 10 400 / 40 / 1250 / – / 341. T-Mobile NSA: 10 400 (12 120) / 320 / 1200 / 210 / 1440. Verizon NSA-HB: 10 500 / 320 / 1280 / 396 / 1907. Verizon NSA-LB: 10 200 (18 800) / 400 / 1100 / 288 / –. T-Mobile 4G: 5000 / 400 / 1300 / 190 / –. Verizon 4G: 10 200 / 300 / 1280 / 265 / – | ms | Table 7. SA also has RRC_INACTIVE for ≈ 5 s | H |
| N14 | **DERIVED** tail energy per burst (P_tail·T_tail) | Verizon 4G 1.8; T-Mobile 4G 0.33; NSA-LB 2.5–2.7; SA 6.2; mmWave 11.5. Compare HU12 2012 LTE: 12.1 | J | calc | L |
| N15 | Power vs throughput | Linear for 4G and 5G, both directions. mmWave slope smallest. Crossover where mmWave beats 4G / low-band: 187 / 189 Mbps (DL), 40 / 123 Mbps (UL) | Mbps | Text | H |
| N16 | Approximate fitted lines (S20U, Verizon, LoS) | Intercepts ≈ 2.8 W (4G), ≈ 3.3 W (LB), ≈ 5.3 W (mmWave). DL slope ≈ 14 mW/Mbps (4G ≈ LB), ≈ 1.6 mW/Mbps (mmWave). UL slope ≈ 68 (4G), ≈ 27 (LB), ≈ 9.5 (mmWave) mW/Mbps | W, mW/Mbps | **Read off Fig. 11** (±20 %). Device totals during iPerf with screen power subtracted; include CPU/network-stack load | M/L |
| N17 | Energy efficiency, 5G vs 4G | At low throughput 5G is up to 79 % (DL) / 74 % (UL) less efficient. At high throughput it is up to 5× (DL) / 2× (UL) more efficient | — | Text | H |
| N18 | Share of total power due to data transfer | mmWave 48–76 % (DL), 46–66 % (UL). 4G 21–53 % (DL), 20–66 % (UL) | % | | H |
| N19 | Signal strength effect, mmWave | Energy per bit decreases monotonically as NR-SS-RSRP rises. Median ≈ 5× lower at [−80, −75] than at [−110, −105] dBm | × | **Read off Fig. 14** (S20U) | M/L |

**Xu et al. 2020 (XU20)**: ZTE Axon10 Pro 5G (Snapdragon 855 + X50 modem), China, NSA sub-6 GHz. Battery-API logger "pwrStrip". Signal held at RSRP −80 to −60 dBm.

| # | Name | Value | Notes | Conf. |
|---|---|---|---|---|
| N20 | Radio share of energy | 5G module 55.18 % of total on average across 4 apps vs screen 30.73 %. 4G module 24.22–50.20 % | Fig. 21 text | H |
| N21 | 5G vs 4G power | 5G power is 2–3× that of 4G. At saturated throughput, 5G energy per bit is a fraction of 4G's (≈ 1/4 as extracted from the text) | Abstract and §6 | H / M |
| N22 | 5G NSA timers | paging DRX 1280 ms; on-duration 10 ms; LTE promotion 623 ms; LTE → NR 1238 ms; NR promotion 1681 ms; inactivity 80 / 100 ms; long C-DRX 320 ms; **tail 10 720 ms (4G) / 21 440 ms (5G)** | Table 7. 5G tail ≈ 20 s because NR and LTE tails chain | H |
| N23 | App-level energy (J), web / video / file | LTE 85.44 / 227.13 / 357.67. NR NSA 113.94 / 140.19 / 157.29. Dynamic 4G ↔ 5G switching 85.41 / 133.66 / 150.80 (−24.8 % vs NSA) | Table 4 (trace-driven) | H |

**Gupta et al. 2024 (GU24)**: Google Pixel 7a, Google On-Device Power Rails Monitor (ODPM, Perfetto), 1-minute uses. S23+ by battery-counter subtraction. **All values read off bar charts (±10 %).**

| # | Name | Value (mW) | Conf. |
|---|---|---|---|
| N24 | Cellular rail, screen-on standby (Fig. 4) | Pixel 7a: 4G ≈ 185, 5G ≈ 335. S23+: 4G ≈ 140, 5G ≈ 185. The text says 5G standby is "about two times" 4G | M |
| N25 | Radio-attributed power, 720p YouTube (Fig. 5, battery-subtraction method) | Pixel 7a: Wi-Fi ≈ 660, 4G ≈ 895, 5G ≈ 1050. S23+: Wi-Fi ≈ 480, 4G ≈ 600, 5G ≈ 675. Upper bound: includes non-radio processing | M |
| N26 | Cellular rail per app (Fig. 2) | Standby ≈ 310; Instagram ≈ 950; YouTube ≈ 730; Clash of Clans ≈ 730; Google Meet video call ≈ 2580 (uplink-heavy) | M |

**Ding et al. 2013 (DI13)**: Nexus One, Motorola Atrix 4G, Sony Xperia S. Currents in mA at 3.7 V. 3G on T-Mobile / AT&T.

| # | Name | Value | Conf. |
|---|---|---|---|
| N27 | 3G state currents vs RSSI (−85 / −95 / −105 dBm), Nexus One | Rx 180 / 192 / 333. Tx 198 / 281 / 390. DCH tail 150 / 160 / 260 mA (3.6 s). FACH tail 97 mA (3.3 s, constant) | H |
| N28 | Same, Atrix | Rx 217 / 287 / 520. Tx 417 / 497 / 530. Tail1 151 / 188 / 411 (4.0 s). Tail2 114 (10.0 s) | H |
| N29 | Same, Xperia S | Rx 155 / 215 / 488. Tx 303 / 425 / 512. Tail1 120 / 190 / 424 (4.0 s). Tail2 80 (10.0 s) | H |
| N30 | Wi-Fi Rx / Tx current vs RSSI (Nexus One, −50 … −90 dBm) | Rx 117 → 89. Tx 165 → 251 (peaks near −70 dBm). PSM tail 70 (Nexus One) / 68 (Atrix) / 50 (Xperia S) mA, independent of signal | H |
| N31 | Energy penalty, 100 KB download (30 ms RTT) | 3G: +6.6 % at −95 and +52.0 % at −105 vs −85 dBm. Wi-Fi: +113.3 % at −85 and +810.5 % at −90 vs −50 dBm | H |
| N32 | Prevalence of weak signal (3785 users) | 3G below −91.7 dBm 47 % of the time on average. Wi-Fi at −80 dBm or below 23–25 % of the time. 43 % (3G) and 21 % (Wi-Fi) of foreground data transferred under poor signal | H |
| N33 | Model-error reduction with a signal-aware model | Wi-Fi 61.0 % → 5.4 %; 3G 52.1 % → 7.2 % | H |

**Other network rows**

| # | Name | Value | Source | Conf. |
|---|---|---|---|---|
| N34 | 802.11n active Rx power spread | 240–795 mW (median 434 mW) across MCS 0–7, TCP/UDP. Energy per bit 13.97–1902 nJ/bit (median 60.15) | SU14 (Nexus S) | H |
| N35 | LTE uplink Tx-RF power vs UE transmit power S_Tx | Slope 5.50 mW/dB for −30 ≤ S_Tx ≤ 10 dBm. **117 mW/dB** for 10 < S_Tx ≤ 23 dBm. Total 2.1 → 3.8 W over that upper range. Rx-RF slope −1.11 mW/dB; Rx BB 2.89 mW/Mbps; Tx BB 0.87 mW/Mbps | JE12 (2012 LTE USB dongle, **not a phone**; use for shape only) | M |
| N36 | 3G cellular data vs Wi-Fi in apps | Web / email: Wi-Fi about 4× more efficient than 3G (≈ 20 % less total energy). Idle connection: 3G ≈ 10 mW vs Wi-Fi ≈ 50 mW | CH13 | H |
| N37 | 3G microbenchmark (S3) | UL 547 ± 149 kbps at 1137 ± 372 mW. DL 1317 ± 814 kbps at 768 ± 64 mW. Efficiency 481 / 1715 kb/J | CH13 | H |

### 1.4 GNSS, Bluetooth, audio, camera, sensors, haptics

| # | Name | Value | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|
| G1 | GPS above idle (S3) | acquisition 386 ± 19.5; tracking 433 ± 21.5 | mW | Idle = 573 mW (screen on, Data Monitor app). Cold fix ≈ 400 s and 154 J; warm (5–10 min old) ≈ 30 s | CH13 | H |
| G2 | GPS module only (u-blox ANTARIS 4) | 143.1 (internal antenna) / 166.1 (external) | mW | Independent of satellite count and signal. "Worst case" | CH10 | H (old) |
| G3 | Energy per 1-s position update | Mi 8 (dual-frequency L1+L5) 318 ± 32; Redmi Note 7 (L1) 232 ± 20 | mJ | ≈ 0.32 / 0.23 W at 1 Hz. Dual-frequency costs +37 % outdoors, +28 % indoors | KW19 | H |
| B1 | Bluetooth (classic) stereo headset audio | +36.0 (30 cm) / +44.9 (10 m) | mW | HTC Dream (G1) | CH10 | H (old) |
| B2 | BLE scanning / advertising | no verified smartphone measurement found | — | — | — | — |
| A1 | Audio playback, headphones, display off | total 226 (S3); 322 (N1); 460 (G1). Max volume adds +5 mW | mW | Whole device | CH13 | H |
| A2 | Audio branch in streaming model | ≈ 0.11 (constant) | W | Fairphone 2 incl. speaker | HE20 | H |
| A3 | Loudspeaker at high volume | no verified source | — | — | — | — |
| C1 | Camera, whole device (S3) | still 2256; video 2614 | mW | 50 % brightness. "Among the highest energy consumers" | CH13 | H |
| C2 | Camera recording, whole device | 1400 (640×360) / 1650 (1280×720) average. AMOLED share 170 / 180. Spikes ≈ 1620 / 2050 | mW | Galaxy S2, 4-min recording, 30 fps | CX13 | H |
| C3 | Camera rail during video call | ≈ 260 | mW | Pixel 7a, Google Meet (**fig. read**) | GU24 | M |
| SN | Sensors at 50 Hz, above "dummy" | accelerometer 5 ± 2.3; gyroscope 30 ± 1.3; light 3 ± 1.7; magnetometer 12 ± 0.6; barometer 1 ± 0.7; proximity 7 ± 2.2 | mW | S3 | CH13 | H |
| HP | Haptics (vibration motor) | no source found; plausible 50–200 mW while active (unverified). Duty cycle is tiny, so ignore | mW | — | unverified | L |

### 1.5 Standby, background, wakelocks

| # | Name | Value | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|
| BG1 | Suspend power (airplane mode) | G1 27, N1 25, S3 24 (suggested "intrinsic lower bound" ≈ 25 mW) | mW | Brightness-dependent display removed | CH13 | H |
| BG2 | Cost of keeping an unused data connection (S3) | 2G/3G ≈ 5–10 mW ("similar to the suspend case"). Wi-Fi ≈ +50 mW above airplane mode | mW | Stated with the idle (awake) measurements; the paper's analysis also gives 3G ≈ 10 mW vs Wi-Fi 50 mW | CH13 | H |
| BG3 | Suspend breakdown (S3, Galaxy S4) | SoC suspend 3.8 / 5.1 mA; cellular paging 8.3 / 2.3 mA; Wi-Fi beacon 1.1 mA (screen off), 3.3 mA (screen on) | mA | **DERIVED** total 7–15 mA ≈ 27–55 mW | CH15 (via course slides) | M |
| BG4 | Daily energy breakdown (1520 users, Galaxy S3/S4) | suspended-state 24.6 %; cellular paging 14.6 %; Wi-Fi beacon 2.7 %; cellular data 11.8 % vs Wi-Fi 1.5 %; GPU 1.7 %; screen 27.4 % (component); CPU idle 9.0 % (CPU idle 75.4 % of CPU time); apps and services 38.1 %; screen-off app/service 12.6 %, screen-on 23.9 % | % | Wi-Fi: 1.4× more time, 1.8× more bytes, 4.2× less energy than cellular | CH15 (via slides) | M |
| BG5 | Foreground app drain rates | 92.6 % of apps draw 200–400 mA on average (≈ 0.75–1.5 W) | mA | Galaxy S3/S4 | CH15 (via slides) | M |
| BG6 | Awake idle vs suspend | Freerunner: idle (backlight off) 268.8 vs suspend 68.6 mW. S3 idle (airplane, display on 50 %) 805 mW; with the brightness part removed, 666 mW | mW | | CH10, CH13 | H (old) |
| BG7 | Ads and tracking in free apps | 65–75 % of app energy in third-party ad modules. Core function only 10–30 %. "Wakelock bugs" found | % | | PA12 | M (abstract) |
| BG8 | AOSP example power-profile values (illustrative only, **not a device**) | cpu.idle (suspend) 3 mA; cpu.awake 50 mA; screen.on 200 mA; screen.full 100–300 mA; wifi.on 2 mA; wifi.active 31 mA; wifi.scan 100 mA; radio.active 100–300 mA; radio.scanning / radio.on 1.2 mA; gps.on 50 mA; audio ≈ 10 mA; video ≈ 50 mA; camera.avg 600 mA; flashlight 200 mA; ambient.on ≈ 100 mA | mA | Documentation examples | AOSP | L |
| BG9 | **DERIVED** standby drain on a 19.3 Wh (5000 mAh × 3.86 V) pack | 30 mW → 0.16 %/h; 60 mW → 0.31 %/h; 100 mW → 0.52 %/h; 200 mW → 1.04 %/h | %/h | calc | — | L |

### 1.6 Power conversion

| # | Name | Value | Notes | Source | Conf. |
|---|---|---|---|---|---|
| PM1 | SMPS efficiency (battery → rails) | ≈ 83 % (scaling factor 1.2 applied to rail power) | S3, determined experimentally | CH13 | H (one device) |
| PM2 | Light-load buck efficiency | falls below ≈ 80 % at < 20 mW load in one smartphone buck design; peak ≈ 92.7 % at 30 mA | Seen only in search snippets of a 2023 preprint (preprints.org 202312.0199). Not opened, so **unverified** | — | L |
| PM3 | Linear regulators for RF rails | Radios partly fed by LDOs (less efficient, lower noise) | CH13 | H (qualitative) |

### 1.7 Device-level totals by activity (measured)

| Activity | Device (year) | P_total | Conditions | Source | Conf. |
|---|---|---|---|---|---|
| Suspend, screen off | G1 / N1 / S3 (2008–12) | 27 / 25 / 24 mW | airplane mode | CH13 | H |
| Screen-on idle | S3 (2012) | 805 mW | 50 % brightness, airplane | CH13 | H |
| Screen-on standby | Pixel 7a (2023) | ≈ 0.88 W (cellular ≈ 0.31, display ≈ 0.35) | ODPM, Fig. 2 | GU24 | M |
| Screen-on standby, 4G / 5G | Pixel 7a; S23+ | ≈ 0.53 / 0.69 W; ≈ 0.72 / 0.77 W | Fig. 4 | GU24 | M |
| Web browsing | S3 | 1275 (Wi-Fi) / 1479 (3G) mW | BBC mobile, 180 s, 50 % brightness | CH13 | H |
| Web browsing | S20 Ultra (2020) | ≈ 1.37 (60 Hz) / 1.70 (120 Hz) W | **DERIVED** from endurance | AT20 | L |
| Email | S3 | 1264 (Wi-Fi) / 1543 (3G) mW | 140 s | CH13 | H |
| Local video (HW decode) | S3 | 1084–1270 mW (SW decode 1571–2329) | 720p H.264 | CH13 | H |
| Streaming 720p YouTube | Pixel 7a | ≈ 1.31 W (cellular, Fig. 2). 1.01 (Wi-Fi) / 1.25 (4G) / 1.40 (5G) W (Fig. 5) | | GU24 | M |
| Streaming 720p YouTube | S23+ | ≈ 1.07 (Wi-Fi) / 1.19 (4G) / 1.26 (5G) W | Fig. 5 | GU24 | M |
| Social feed scrolling (Instagram) | Pixel 7a | ≈ 1.88 W | cellular | GU24 | M |
| Casual game (Clash of Clans) | Pixel 7a | ≈ 1.54 W | cellular | GU24 | M |
| 2D / 3D game | S3 | 1.52 / 2.43 W | Angry Birds / Need for Speed | CH13 | H |
| Games, 2010–12 Samsung phones | various | 1.14–1.75 W (no input); 1.64–2.22 W (with input) | | CX13 | H |
| Video call (Google Meet) | Pixel 7a | ≈ 4.16 W (cellular ≈ 2.58, compute ≈ 0.50, camera ≈ 0.26) | cellular uplink-heavy | GU24 | M |
| Voice call | S3 | 854–865 mW | 3G | CH13 | H |
| Audio, screen off | S3 | 226 mW | headphones | CH13 | H |
| Camera still / video | S3 | 2.26 / 2.61 W | | CH13 | H |
| GPS tracking, screen on | S3 | ≈ 1.0 W (573 + 433 mW) | Data Monitor app, no map | CH13 | H |
| Navigation (map + GNSS + data, screen on) | modern (est.) | 1.2–2.5 W | **DERIVED** sum: display 0.4–0.8, GNSS 0.25–0.45, SoC render 0.3–0.8, cellular 0.1–0.5 W | calc | L |
| Peak device power | S3 | > 7 W | sum of component maxima | CH13 | H |

---

## 2. Key equations and model forms (ODE-ready)

### 2.1 Coupling to the battery ODE

Phone loads are regulated by the PMIC, so the battery sees an approximately **constant-power** load. This is a standard modelling assumption, not separately verified here.

```
I(t)      = P_total(t) / V_term(SOC, I, T)
dSOC/dt   = − I(t) / (3600 · Q_eff)
P_total   = Σ_i P_rail,i / η_i(load)    (η ≈ 0.83, CH13), or use battery-terminal measurements directly
```

All §1.7 totals and the HU12, NA21, DI13 and GU24 numbers are battery-terminal values, so **do not divide them by η again**.

### 2.2 Display

**OLED** (CH13, CH10, DZ11, CX13, DH21):

```
P_OLED(b, content) = s_on · [ P_d0 + A · p_A · g(b) · APL_lin ] · (1 + ρ_RR)
APL_lin = (1/N) Σ_pixels Σ_c w_c (x_c/255)^2.2      (linear-light average picture level, 0..1)
g(b)    ≈ 0.16 + 1.84·b   (b = brightness-slider fraction; matches 0.16 / ≈1.0 / 2.0 of CH13)
```

* P_d0 ≈ 0.1–0.25 W.
* p_A is the content-dependent power density at 50 % brightness and APL = 1 (S3: ≈ 9–18 mW/cm²; modern panels lower).
* A is the active area in cm²: 6.1" ≈ 91 cm², 6.7" ≈ 110 cm² at 19.5:9.
* Typical APL: light UI and web ≈ 0.6–0.8 (65 % white pixels, DZ11); dark mode ≈ 0.1–0.2 (assumption).
* Dash & Hu (DH21) show that a per-pixel linear model is **not monotone/accurate** on modern wide-gamut panels and needs a piecewise model. The APL form above is a lumped approximation.

**LCD** (CH10, HE20): power is almost content-independent.

```
P_LCD(b) = s_on · [ P_panel + c1·b + c2·b² ],   with c2 > 0 (convex)
```

**Refresh rate** (DERIVED from AT20):

```
ΔP_RR ≈ κ_RR · (f_r − 60 Hz),   κ_RR ≈ 2.7–3.0 mW/Hz on a fixed-rate 2020 AMOLED panel
```

LTPO variable-refresh panels reduce this. Vendor claims of 16–46 % lower display power are unverified.

### 2.3 SoC (CPU/GPU) with DVFS

Linux Energy Model (LX):

```
P_CPU(t) = Σ_k [ n_k(t) · u_k(t) · C_k · V_k(f_k)² · f_k  +  n_k · P_leak,k(V_k, T) ] + P_uncore
V_k(f)   ≈ V_0 + s·f      (MT8183 A73: 0.70 V at 0.793 GHz → 1.05 V at 1.989 GHz, s ≈ 0.29 V/GHz)
```

Hence along the DVFS curve P_dyn ∝ f·V(f)² ≈ f^1.9, while performance is ∝ f. Energy per operation therefore rises roughly ∝ f^0.9 at high OPPs.

GPU: same form, or lump it as P_GPU = P_GPU,max · u_GPU^γ. S3 max ≈ 1.4 W; GPU energy is only 1.7 % of daily energy for typical users (CH15).

**Thermal cap** (CHI23, BH19):

```
P_SoC(t) ≤ P_sust(T_amb) ≈ P_sust,25°C − 0.16 W/°C · (T_amb − 25 °C)
```

This is our linearisation of the 0.16 W/°C skin sensitivity. When the cap binds, throttling cuts f (performance −34 % in BH19). If thermal dynamics are modelled, use a first-order lumped RC: C_th·dT/dt = P_total − (T − T_amb)/R_th.

### 2.4 Network: radio resource control (RRC) state machine, ODE-friendly

Exact state machine (HU12, NA21):

* **IDLE**: average P_idle,DRX ≈ P_on,idle · T_on,i / T_DRX,i. LTE 2012 ≈ 19.7 mW; modern measurements in BG3.
* **PROMOTION**: P_pro for T_pro (LTE 0.26 s; 5G NSA 1.4–1.9 s; SA 0.34 s).
* **CONNECTED, active**: P = β + α_u·t_u + α_d·t_d (HU12 linear; NA21 also linear per band).
* **TAIL**: P_tail for T_tail after the last packet. The timer resets on every packet.

**Continuous surrogate** (our derivation). Let bursts arrive as a Poisson process of rate λ (bursts/s) with short durations. The fraction of time the radio sits in CONNECTED (tail) is 1 − e^(−λ·T_tail). The mean radio power is then:

```
P̄_net = P_idle + (P_tail − P_idle)·(1 − e^(−λ T_tail)) + λ·E_pro + Σ_dir α_dir·t̄_dir + (β − P_tail)·d_active
```

Here d_active is the duty cycle of actual transfer. For an ODE state, use r(t) ∈ [0, 1] with dr/dt = −r/T_tail + k_burst·u(t)·(1 − r). An exponential with time constant τ = T_tail reproduces the per-burst tail energy P_tail·T_tail exactly.

This makes background sync frequency an explicit, physically meaningful parameter. Example (HU12 numbers): λ = 1/60 s⁻¹ gives 1 − e^(−11.6/60) = 18 % of time in tail, about +190 mW.

**Signal-strength penalty** (DI13, NA21, JE12; uplink power control background from 3GPP TS 36.213, **not verified in this session**):

```
S_Tx  = min(P_max, P_0 + α_PL·(RS_power − RSRP) + 10·log10(M_PRB))      (dBm; P_max = 23 dBm class-3 UE, unverified)
P_TxRF ≈ a1·S_Tx (S_Tx ≤ 10 dBm) ;  ≈ P_10 + 0.117 W/dB·(S_Tx − 10) (10 < S_Tx ≤ 23 dBm)       (JE12 dongle shape)
P_cell,active(RSRP) ≈ P_cell,active,good · m_sig(RSRP),
m_sig = 1 + γ·max(0, RSRP_ref − RSRP)/10 dB
```

Calibrate γ so that m_sig(−105) / m_sig(−85) ≈ 1.5 (DI13 3G energy +52 %; current ratios 1.2–3.5×). Retransmissions and lower MCS also stretch transfer time, so energy per bit grows faster than power. Wi-Fi shows a cliff below ≈ −85 dBm (+810 % at −90 dBm).

### 2.5 GNSS

```
P_GNSS(t) = P_acq for t < TTFF  (TTFF: cold ≈ 400 s, warm ≈ 30 s, CH13), then P_track
```

* P_acq ≈ 0.39 W and P_track ≈ 0.43 W (S3, including app/CPU).
* Modern energy per fix ≈ 0.23–0.32 J at 1 Hz (KW19).
* Duty cycling with D = fix interval / on-time scales the average power ∝ 1/D plus the re-acquisition cost.

### 2.6 Background and standby

```
P_bg(t) = P_sleep + P_paging(RAT) + P_beacon + Σ_j ν_j·E_wake,j + ν_sync·(E_pro + P_tail·T_tail) + s_wl(t)·ΔP_awake
```

* ν_j are wake-up rates.
* s_wl is a wakelock indicator, which captures "no-sleep" bugs (PA12).
* With HU12 LTE numbers, a single sync costs ≈ 12.8 J. With modern 4G tails (NA21) it costs ≈ 0.3–1.8 J.

---

## 3. Key findings

1. **Display and cellular radio dominate energy; on 5G phones the radio can exceed the screen.**
   * Screen is ≈ 27 % of daily energy by component in 1520 Galaxy S3/S4 users, and 58.8 % of energy is spent during screen-on periods (CH15).
   * On a 2019 5G phone the 5G module averaged 55 % vs the screen's 31 % (XU20).
   * On a 2023 Pixel 7a, cellular was the largest rail in every app tested (GU24).
2. **OLED power is strongly content- and brightness-dependent; LCD is not.**
   * N1: white adds 0.19 W at minimum and 1.31 W at maximum brightness, while black costs only the fixed overhead (CH10).
   * Brightness scales the content part roughly affinely, 0.16 → 2.0× (CH13).
   * Dark mode saves **only 3–9 % at typical 30–50 % brightness but 39–47 % at 100 %** (DH21). It is a strong lever mainly outdoors or at high brightness. This "surprisingly small" saving is a good discussion point.
3. **A 120 Hz refresh rate costs ≈ 0.16–0.18 W on a fixed-rate panel even on a black screen.** That shortened web-browsing life by ≈ 19 % on the S20 Ultra (AT20). Adaptive (LTPO) refresh is the OS-level fix.
4. **CPU power is superlinear in frequency.** Along a real DVFS table, P ∝ f^≈1.9 (MT8183 A73: 2.5× frequency → 5.6× power). At equal V and f, little cores have a coefficient C ≈ 1/2.5–1/4.4 that of big cores (LX). Idle application-processor power is negligible (CH13), so it is energy-optimal to schedule light work on little cores at low OPPs. "Race-to-idle" pays only if it shortens radio and screen-on time.
5. **Peak power far exceeds what a phone can sustain.** The S3 can exceed 7 W, but a phone can shed only ≈ 2.2–3.8 W at a 45 °C skin limit (CHI23). Sustained gaming is therefore thermally capped, and throttling cut performance by up to 34 % (BH19). A cap of P_sust is needed in any long-gaming time-to-empty (TTE) scenario.
6. **Tail energy, not bytes, dominates small transfers.**
   * Old LTE: one packet cost 12.76 J vs 0.04 J on Wi-Fi (HU12).
   * Modern 4G tails are much cheaper: 66–178 mW for 5–10 s vs 1060 mW for 11.6 s in 2012 (NA21 vs HU12). **Using 2012 LTE constants would overstate modern background drain by ≈ 5–35×.**
   * 5G tails are again expensive: mmWave 1.09 W; SA 0.59 W; NSA tails can chain to ≈ 21 s (XU20).
   * Periodic small syncs over 5G are especially wasteful. "Stay on 4G unless high throughput is needed" saved 24.8 % in XU20.
7. **Weak signal is a large hidden multiplier.**
   * 3G energy for the same download rose 52 % from −85 to −105 dBm, and Wi-Fi energy rose 810 % from −50 to −90 dBm (DI13).
   * Users spend about half their 3G time below −91.7 dBm (DI13).
   * mmWave energy per bit improves ≈ 5× from poor to good RSRP (NA21).
   * Weak signal is a prime candidate for "conditions that drain most".
8. **Wi-Fi is far cheaper than cellular per byte.** Cellular drained 4.2× more energy while carrying fewer bytes (CH15). Wi-Fi was ≈ 4× more efficient for web and email (CH13). Preferring Wi-Fi and deferring bulk background transfers to good-signal or Wi-Fi periods is a robust recommendation.
9. **Things that matter surprisingly little** (CH13, CH10):
   * Environmental sensors: accelerometer 5 mW, gyroscope 30 mW.
   * Flash storage: < 8 mW.
   * RAM: ≈ 6 % of total.
   * Audio volume: +5 mW at max.
   * Classic Bluetooth audio: ≈ 40 mW.
   * GPU on non-game apps: 1.7 % of daily energy (CH15).
   * Display power during video and camera use is also a smaller share than people expect (CX13).
10. **Things that cost surprisingly much:**
    * Camera use: 2.3–2.6 W on the S3.
    * Video calls: ≈ 4.2 W on a Pixel 7a over cellular, mostly uplink radio.
    * Cold GPS fixes: 154 J.
    * Cellular paging in standby: 14.6 % of daily energy (CH15).
    * Third-party ad modules: 65–75 % of free-app energy (PA12).
11. **The standby floor has stayed at ≈ 25 mW since 2008** (CH13). Modern standby drain is set mainly by radio paging and background wakeups: ≈ 0.15–0.5 %/h on a 19 Wh pack (DERIVED).
12. **Power conversion.** About 17 % of rail power is lost in DC-DC conversion on the S3 (CH13). Light-load efficiency is worse, which slightly raises the effective standby cost. Battery-terminal measurements already include this loss.
13. **Measurement method matters.** Android fuel-gauge (software) readings under-report power by 5–20 % (NA21). Subtraction methods attribute cross-component effects to the subtracted component; uplink is poorly captured (GU24).

---

## 4. Caveats

* **Device generation.** Many high-confidence component numbers come from 2008–2013 phones (Freerunner, G1, N1, S3, HTC LTE, Nexus One). Panels, modems and SoC nodes have improved since. Use them for **model structure and ratios**; anchor absolute levels to modern sources (GU24, NA21, XU20) and widen UQ ranges.
* **Figure-read values.** GU24 bar heights and NA21 Fig. 11/14 slopes are read off figures (±10–20 %). The GU24 poster reports no tabulated numbers, and its brightness and network conditions are only partly specified.
* **Secondary sources.**
  * CH15 numbers come from a WPI course presentation summarising the SIGMETRICS paper; the original PDF was not retrievable (ACM 403).
  * DH21 numbers come from Purdue's press release, not the paper.
  * AnandTech 120 Hz numbers (AT20) come from secondary news reports, which disagree with each other (14.05 h vs 13.04 h at 60 Hz). The original site is offline.
  * Treat all of these as M or L.
* **Carrier and network configuration.** RRC timers and tail powers depend on the carrier and network configuration (NA21 Table 7: T-Mobile 4G 5 s vs Verizon 10.2 s), and on modem generation (GU24: Exynos 5300 vs Qualcomm X70 differ markedly). Treat tail parameters as scenario inputs.
* **Brightness scale.** Our b is the **slider fraction**. Android maps the slider to luminance nonlinearly, and high-brightness mode (HBM, outdoor) can exceed "max" substantially. HBM is not covered by any source here.
* **Kernel energy-model coefficients.** The dynamic-power-coefficient values (LX) are vendor-supplied estimates for scheduling. They exclude leakage, uncore, interconnect and DRAM, so true per-core power at a given OPP is higher. The Snapdragon 865 (SM8250) voltages are not in mainline device trees, so only MediaTek and Rockchip OPP voltages were used for the DERIVED P(f) curves.
* **Thermal numbers.** CHI23 is a trade-journal article (Electronics Cooling, 2023) with design-specific results. A widely quoted "≈ 3.5 W sustainable envelope" appears in a Southampton accepted manuscript that we could not open (HTTP 403), so it is not cited as a value.
* **Missing data.** No verified smartphone measurements were found for haptics, loudspeaker at high volume, BLE scanning, NFC, or modern Wi-Fi 6/7 power versus throughput. Values for these are marked unverified or omitted.
* **Battery and nominal voltage.** The S20 Ultra's 5000 mAh capacity comes from secondary spec listings. 3.86 V nominal is an assumption (unverified); it affects only DERIVED rows.
* **Data licensing.** Most measured values come from copyrighted papers (ACM/IEEE). They are used as **published measurements with citation**, which the problem statement allows. They are not redistributed datasets.
  * Fully open-licensed sources: HE20 (CC BY 4.0), Linux kernel device trees (GPL-2.0 / BSD-3-Clause / MIT, per the SPDX headers), and arXiv preprints (arXiv licence).
  * The NA20 "5Gophers" dataset is reported as CC BY 4.0 by its IEEE DataPort listing (seen only in a search snippet, unverified).
* **Our own derivations.** The Poisson tail formula, the m_sig signal multiplier and the refresh-rate slope κ_RR are our constructions on top of cited data. Label them as modelling assumptions in the paper.

---

## 5. References

Access and licence notes describe what we observed; "©" means the publisher holds copyright and only cited numbers are used.

1. **[CH10]** A. Carroll, G. Heiser. "An Analysis of Power Consumption in a Smartphone." *USENIX ATC 2010*, pp. 271–284. PDF: https://www.usenix.org/legacy/event/usenix10/tech/full_papers/Carroll.pdf (USENIX open-access proceedings).
2. **[CH13]** A. Carroll, G. Heiser. "The Systems Hacker's Guide to the Galaxy: Energy Usage in a Modern Smartphone." *APSys 2013*, Singapore. DOI: 10.1145/2500727.2500734. Author PDF: https://www.trustworthy.systems/publications/nicta_full_text/7044.pdf (free author copy, © ACM).
3. **[HU12]** J. Huang, F. Qian, A. Gerber, Z. M. Mao, S. Sen, O. Spatscheck. "A Close Examination of Performance and Power Characteristics of 4G LTE Networks." *MobiSys 2012*, pp. 225–238. DOI: 10.1145/2307636.2307658. Author PDF: https://feng-qian.github.io/paper/lte_mobisys12.pdf (free author copy, © ACM).
4. **[NA21]** A. Narayanan, X. Zhang, R. Zhu, A. Hassan, S. Jin, et al. "A Variegated Look at 5G in the Wild: Performance, Power, and QoE Implications." *ACM SIGCOMM 2021*, pp. 610–625. DOI: 10.1145/3452296.3472923. NSF PAR copy: https://par.nsf.gov/servlets/purl/10291906. Artifact: https://github.com/SIGCOMM21-5G/artifact (free; © ACM; artifact licence not checked).
5. **[NA20]** A. Narayanan, E. Ramadan, J. Carpenter, Q. Liu, Y. Liu, F. Qian, Z.-L. Zhang. "A First Look at Commercial 5G Performance on Smartphones." *WWW 2020*, pp. 894–905. DOI: 10.1145/3366423.3380169. arXiv: https://arxiv.org/abs/1909.07532. Background only; no power values used. Dataset "5Gophers" reported as CC BY 4.0 (unverified).
6. **[XU20]** D. Xu, A. Zhou, X. Zhang, G. Wang, X. Liu, C. An, Y. Shi, L. Liu, H. Ma. "Understanding Operational 5G: A First Measurement Study on Its Coverage, Performance and Energy Consumption." *ACM SIGCOMM 2020*. DOI: 10.1145/3387514.3405882. PDF: https://xyzhang.ucsd.edu/papers/DXu_SIGCOMM20_5Gmeasure.pdf (free author copy, © ACM).
7. **[GU24]** A. Gupta, A. Heidari, A. Kalipattapu, I. K. Jain, D. Bharadia. "3 W's of Smartphone Power Consumption: Who, Where and How Much Is Draining My Battery?" *ACM MobiCom 2024* (poster/abstract, 3 pp.). DOI: 10.1145/3636534.3695905. PDF: https://par.nsf.gov/servlets/purl/10598540. Project: https://wcsng.ucsd.edu/ue-power/ (free; © owner/authors).
8. **[DI13]** N. Ding, D. Wagner, X. Chen, A. Pathak, Y. C. Hu, A. Rice. "Characterizing and Modeling the Impact of Wireless Signal Strength on Smartphone Battery Drain." *ACM SIGMETRICS 2013*, pp. 29–40. DOI: 10.1145/2494232.2466586 (DOI as given by a BibSonomy record). PDF: https://www.cl.cam.ac.uk/~acr31/pubs/ding-signalstrength.pdf (free author copy, © ACM).
9. **[CH15]** X. Chen, N. Ding, A. Jindal, Y. C. Hu, M. Gupta, R. Vannithamby. "Smartphone Energy Drain in the Wild: Analysis and Implications." *ACM SIGMETRICS 2015*, pp. 151–164. DOI: 10.1145/2745844.2745875. Numbers taken from the course presentation by M. Zhan, WPI CS528 (2016): https://web.cs.wpi.edu/~emmanuel/courses/cs528/S16/student_talks/smart_energy_drain_Miao_Zhan_wk12.pdf (secondary; the original paper was not accessible).
10. **[PA12]** A. Pathak, Y. C. Hu, M. Zhang. "Where Is the Energy Spent Inside My App? Fine Grained Energy Accounting on Smartphones with Eprof." *EuroSys 2012*, pp. 29–42. DOI: 10.1145/2168836.2168841. PDF: https://www.cs.columbia.edu/~nieh/teaching/e6998_s13/papers/eprof_eurosys2012.pdf (abstract-level values only; © ACM).
11. **[DZ11]** M. Dong, L. Zhong. "Chameleon: A Color-Adaptive Web Browser for Mobile OLED Displays." *MobiSys 2011*; Rice tech report 2010-1212, arXiv:1101.1240. https://arxiv.org/abs/1101.1240 (arXiv, free to read).
12. **[DH21]** P. Dash, Y. C. Hu. "How Much Battery Does Dark Mode Save? An Accurate OLED Display Power Profiler for Modern Smartphones." *MobiSys 2021*. Values from the Purdue University news release, 28 July 2021: https://www.purdue.edu/newsroom/releases/2021/Q3/dark-mode-may-not-save-your-phones-battery-life-as-much-as-you-think,-but-there-are-a-few-silver-linings.html (secondary/press). The related patent US 11,532,271 describes the piecewise per-pixel model.
13. **[CX13]** X. Chen, Y. Chen, Z. Ma, F. C. A. Fernandes. "How Is Energy Consumed in Smartphone Display Applications?" *HotMobile 2013*. PDF: https://hotmobile.org/2013/papers/full/17.pdf (free; © ACM).
14. **[HE20]** C. Herglotz, S. Coulombe, C. Vazquez, A. Vakili, A. Kaup, J.-C. Grenier. "Power Modeling for Video Streaming Applications on Mobile Devices." *IEEE Access* 8:70234–70244, 2020. DOI: 10.1109/ACCESS.2020.2986580. **CC BY 4.0 (open access).**
15. **[KW19]** B. Karki, M. Won. "Characterizing Power Consumption of Dual-Frequency GNSS of Smartphone." arXiv:1910.13041 (v3, 2021; reported as IEEE GLOBECOM 2020, venue unverified). https://arxiv.org/abs/1910.13041 (arXiv).
16. **[SU14]** L. Sun, R. K. Sheshadri, W. Zheng, D. Koutsonikolas. "Modeling WiFi Active Power/Energy Consumption in Smartphones." *IEEE ICDCS 2014*, pp. 41–51. DOI: 10.1109/ICDCS.2014.13 (© IEEE; read from a copy hosted at typeset.io).
17. **[JE12]** A. R. Jensen, M. Lauridsen, P. Mogensen, T. B. Sørensen, P. Jensen. "LTE UE Power Consumption Model for System Level Energy and Performance Optimization." *IEEE VTC-Fall 2012*. Aalborg University repository: https://vbn.aau.dk/ws/files/73426863/vtcFall2012.pdf (accepted manuscript, © IEEE). Measured on an LTE USB dongle.
18. **[BH19]** G. Bhat, S. Gumussoy, U. Y. Ogras. "Power and Thermal Analysis of Commercial Mobile Platforms: Experiments and Case Studies." *DATE 2019*; arXiv:1904.09814. https://arxiv.org/abs/1904.09814 (arXiv).
19. **[CHI23]** V. Chiriac. "A Figure of Merit for Smart Phone Thermal Management." *Electronics Cooling*, 19 Sept 2023. https://www.electronics-cooling.com/?p=35382 (trade magazine, free to read).
20. **[AT20]** A. Frumusanu. "Samsung Galaxy S20+ & Ultra (Snapdragon & Exynos) Battery Life Preview." AnandTech, 23 Mar 2020 (https://www.anandtech.com/show/15647, now offline/redirected), and the S20 review, display page /show/15603/…/10. Numbers via secondary reports: 9to5Google, 24 Mar 2020, https://9to5google.com/2020/03/24/galaxy-s20-exynos-battery-120hz-display/; Ubergizmo, 25 Mar 2020, https://www.ubergizmo.com/2020/03/samsung-galaxy-s20-120hz-display-battery/. The 160–180 mW figure comes from search snippets of those secondary reports; **verify before quoting in the paper**.
21. **[LX]** Linux kernel v6.6:
    * Documentation/power/energy-model.rst
    * Documentation/devicetree/bindings/arm/cpus.yaml (`dynamic-power-coefficient`, µW/MHz/V²)
    * arch/arm64/boot/dts/qcom/sdm845.dtsi (GPL-2.0)
    * arch/arm64/boot/dts/qcom/sm8250.dtsi (BSD-3-Clause)
    * arch/arm64/boot/dts/mediatek/mt8183.dtsi (GPL-2.0 OR MIT)
    * arch/arm64/boot/dts/rockchip/rk3399.dtsi and rk3399-opp.dtsi (GPL-2.0+ OR MIT)

    Source: https://github.com/torvalds/linux/tree/v6.6 (**open licence**).
22. **[AOSP]** Android Open Source Project. "Power profiles for Android" (power_profile.xml values): https://source.android.com/docs/core/power/values. Documentation example values only. Licence assumed to be the standard AOSP docs licence, not checked.
23. **[SPEC-S20U]** Samsung Galaxy S20 Ultra battery: 5000 mAh typical / 4855 mAh rated, via retailer and spec aggregators (e.g., https://ao.com/brands/features-we-love/samsung-galaxy-s20). Secondary; used only for a DERIVED estimate.

### Source-key → reference number map

CH10 = 1, CH13 = 2, HU12 = 3, NA21 = 4, NA20 = 5, XU20 = 6, GU24 = 7, DI13 = 8, CH15 = 9, PA12 = 10, DZ11 = 11, DH21 = 12, CX13 = 13, HE20 = 14, KW19 = 15, SU14 = 16, JE12 = 17, BH19 = 18, CHI23 = 19, AT20 = 20, LX = 21, AOSP = 22, SPEC-S20U = 23.

*Downloaded working copies (untrusted, not committed) are in /tmp/claude-0/research_dl/device_power/.*
