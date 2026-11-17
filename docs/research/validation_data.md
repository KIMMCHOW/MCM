# Open validation datasets and published device-level benchmarks: research note

Research note for 2026 MCM Problem A (smartphone battery drain).

**Scope.** This note covers:

* open-licensed cell-level and phone-level datasets that can **validate** (not fit) our continuous-time SOC/TTE model;
* published device battery-life measurements and specifications that can be cited;
* a check of whether each source can be downloaded from this container;
* a concrete validation plan.

Compiled 2026-11-17.

**Where the numbers come from.** Every number either comes from a source in the reference list (section e),
cited by key, or was **computed in this session** (tag `[RUN]`) from data we downloaded. The scripts are
`docs/research/scripts/validation_*.py` and their outputs are in `docs/research/data/validation/` (see the
README there). Values that rest only on background knowledge are marked **unverified**.

**Conventions.**

* I > 0 means discharge.
* "TTE" is time from full charge to the device's shutdown, or to the test's cutoff voltage.
* P_avg = E/TTE is the mean battery-side power implied by a runtime. Here E is the rated pack energy.

**Confidence levels.**

* **H**: read in a primary or official document, or computed by us from the raw open data.
* **M**: from a secondary source (press or review site, review article), or a derived value that depends on a
  secondary input such as a regulatory Wh figure reported by press.
* **L**: our inference, a heuristic, or unverified.

---

## 0. Executive summary

**Verdict: the best open validation sources.**

1. **Bonn iPhone sensor-battery dataset** [BONN26, SCH26]. This is the best phone-level source. It is
   CC BY 4.0 and 112.5 MB, and we downloaded it here. It holds iOS battery level (1 % steps) against time for
   9 iPhones (3 × 13 Pro, 6 × 14 Pro) under controlled scenarios: idle, GPS, rear camera, camera + GPS, and
   motion sensors. In every run the OLED screen is on and shows black, and airplane mode is on. These are the
   only open, documented, full discharge curves of modern phones we found.
2. **LG HG2 cell dataset** (Kollmeyer et al.) [KOL20]. This is the best cell-level source for temperature
   effects. It is CC BY 4.0, and we downloaded a subset here. The cell is a 3 Ah 18650 tested at
   −20/−10/0/10/25/40 °C with C/20, 0.5C and 2C discharges, HPPC pulses and automotive drive cycles. Use it to
   validate the OCV/ECM/Arrhenius temperature part of the model.
3. **NASA Ames PCoE Battery Data Set** [NASA07]. This is the best cell-level source for aging and temperature.
   It is a NASA public release with no licence stated (see caveats). It is 209.7 MB and we downloaded it in
   full. It has 34 cells of 2 Ah 18650 at 4/24/43 °C, with capacity fade over cycling, EIS-fitted Re/Rct, and
   square-wave pulsed loads.

Supplementary sources:

* the Oxford Degradation Dataset 1 (ODbL; a 740 mAh pouch cell under drive-cycle load) [OX17];
* the Panasonic 18650PF set (CC BY 4.0) [KOL18];
* CALCE (open access, no licence stated) [CALCE].

Published runtimes from Apple, Google, Notebookcheck and Laptop Mag (section a.4) give a cross-device sanity
band. They are cited facts, not open data.

**Headline numbers computed in this session `[RUN]`:**

* **Black-screen idle, iPhone.** With the screen on and showing black and the radios off, the median TTE is
  **55.0 h on the iPhone 14 Pro** and **50.8–51.5 h on the 13 Pro**. That is P_avg ≈ **0.225 W** and
  **0.236 W**.
* **Rear camera.** The camera (frames not rendered) cuts TTE to **8.4 h / 9.5 h**. That is P ≈ **1.47 / 1.26 W**.
* **Component powers add up, and capacity does not depend on the load rate here.** The camera + GPS run
  follows 1/TTE_combo = 1/TTE_cam + 1/TTE_GPS − 1/TTE_idle to within **0.15 %** on both models. This holds
  across a 6.5× range of load power. It strongly supports an additive component power model feeding an
  energy-bucket or ECM battery at phone-level loads.
* **iOS holds "100 %" for a long time.** The display stays at "100 %" for **≈ 6–7 % of the run**. After that,
  displayed SOC falls linearly: the slope changes by ≤ 6 % between the 99–50, 50–20 and 20–1 % windows.
* **Cold costs capacity, and more so at high rate.** On the LG HG2 cell at 0.5C, the delivered capacity at
  −20 °C is **0.73 ×** the 25 °C value and the energy is **0.67 ×**. At 2C the figures are **0.49 ×** and
  **0.45 ×**. At 0 °C they are 0.91 × / 0.89 × (0.5C).
* **Apparent resistance follows Arrhenius.** The apparent DC resistance (C/20 vs 0.5C) rises from
  **0.038 Ω at 25 °C to 0.25 Ω at −20 °C**. That is an Arrhenius fit with **E_a ≈ 26 kJ/mol**
  (E_a/R_g ≈ 3170 K).
* **NASA 18650 cells at 2 A.** They deliver **1.84–2.03 Ah in 55–61 min at 24 °C**, against **1.29–1.38 Ah in
  ≈ 39–42 min at 4 °C** (median of the first 5 full cycles). B0005 loses **28.7 %** of its capacity over
  168 cycles.

---

## (a) Parameter tables

### a.1 Cell-level validation quantities (from open datasets)

| Name | Symbol | Value | Range | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Delivered capacity, fresh NASA 18650, 24 °C, ~2 A CC (1C) | Q_24 | 1.86 (B0005, first cycle) | 1.84–2.03 (median of first 5 full cycles, B0005/6/7/18) | A h | CC-CV charge 1.5 A to 4.2 V, cut at 20 mA. Discharge 2 A to 2.7/2.5/2.2/2.5 V | NASA07 [RUN] | H |
| Time to cutoff, same runs | TTE_24 | 55.8 (B0005, first cycle) | 55–61 | min | as above | NASA07 [RUN] | H |
| Delivered energy, same runs | E_24 | 6.61 (B0005) | 6.56–7.20 | W h | as above | NASA07 [RUN] | H |
| Energy-weighted mean discharge voltage, same runs | V̄ | 3.55 | 3.54–3.56 | V | as above | NASA07 [RUN] | H |
| Cell surface temperature rise during a 1C discharge in the chamber | ΔT_cell | 14.1 (B0005) | 13.8–16.1 | K | ambient 24 °C, from 24.3 °C up to 39.0 °C (B0005) | NASA07 [RUN] | H |
| Delivered capacity at 4 °C, ~2 A | Q_4 | 1.30 | 1.29–1.38 (B0053–56, median of first 5) | A h | ambient 4 °C. Cutoffs 2.0–2.7 V. The README warns of unexplained low-capacity runs | NASA07 [RUN] | M |
| Delivered capacity at 4 °C, ~1 A | Q_4,1A | 1.53 | 1.50–1.56 (B0046–48) | A h | ambient 4 °C | NASA07 [RUN] | M |
| Capacity ratio 4 °C / 24 °C at ~2 A | Q_4/Q_24 | 0.71 | 0.64–0.75 | — | DERIVED from the two rows above. Cell and cutoff differ between groups | NASA07 [RUN] | M |
| Delivered capacity at 43 °C, ~4 A (2C) | Q_43 | 1.85 | 1.79–1.87 (B0029–32, median of first 5) | A h | TTE 27–28 min. Cell temperature rises 15–18 K | NASA07 [RUN] | H |
| EIS-fitted electrolyte resistance, fresh, 24 °C | R_e | 0.045 (B0005) | 0.038–0.065 | Ω | NASA fitted parameter, EIS 0.1 Hz–5 kHz. May include fixture | NASA07 [RUN] | M |
| EIS-fitted charge-transfer resistance, fresh, 24 °C | R_ct | 0.070 (B0005) | 0.062–0.096 | Ω | as above | NASA07 [RUN] | M |
| EIS-fitted R_ct at 4 °C (first EIS) | R_ct,4 | 0.22 | 0.15–0.25 (B0041–56) | Ω | ambient 4 °C | NASA07 [RUN] | M |
| Capacity fade over 168 cycles, B0005 | 1−Q_N/Q_0 | 28.7 % | 24.1–41.1 % (B0005/6/7); 27.2 % for B0018 (132 cycles) | — | 1C cycling at 24 °C to 2.7–2.2 V | NASA07 [RUN] | H |
| Re growth over the same life | ΔR_e/R_e | +12 % (B0005) | +1 % to +76 % | — | Very cell dependent, while capacity fell 24–41 % | NASA07 [RUN] | M |
| LG HG2 C/20 capacity at 25 °C (≈ thermodynamic capacity) | Q_C/20 | 2.78 | — | A h | 0.15 A to 2.8 V. Mean V 3.705 V, E = 10.30 Wh | KOL20 [RUN] | H |
| HG2 0.5C capacity at 25 °C | Q_0.5C | 2.75 | — | A h | 1.5 A to 2.8 V, 110 min, E = 10.05 Wh | KOL20 [RUN] | H |
| HG2 capacity ratio Q(T)/Q(25 °C), 0.5C | f_Q(T) | 40 °C 0.97; 10 °C 0.93; 0 °C 0.91; −10 °C 0.86; −20 °C 0.73 | — | — | Single cell, tested in the order 25→40→10→0→−10→−20 °C, so some aging is mixed in (≈ 3 % at C/20, 40 °C) | KOL20 [RUN] | H |
| HG2 energy ratio E(T)/E(25 °C), 0.5C | f_E(T) | 0.98 / 0.93 / 0.89 / 0.83 / 0.67 | — | — | same temperatures as above | KOL20 [RUN] | H |
| HG2 capacity ratio, 2C (6 A) | f_Q,2C(T) | 0.98 / 0.93 / 0.88 / 0.73 / 0.49 | — | — | Cell self-heats by up to 6 K at −20 °C | KOL20 [RUN] | H |
| HG2 capacity ratio, C/20 | f_Q,C/20(T) | 0.97 / 0.95 / 0.93 / 0.92 / 0.90 | — | — | Shows that the cold-temperature loss is mostly kinetic (rate-dependent) | KOL20 [RUN] | H |
| HG2 apparent DC resistance (V̄_C/20 − V̄_0.5C)/ΔI | R_app(T) | 25 °C 0.038; 10 °C 0.066; 0 °C 0.098; −10 °C 0.158; −20 °C 0.250; 40 °C 0.028 | — | Ω | DERIVED. Lumps ohmic, polarisation and SOC-window effects | KOL20 [RUN] | M |
| Arrhenius activation energy of R_app (T ≤ 25 °C) | E_a | 26.4 (E_a/R_g = 3169 K) | single-temperature values 3096–3203 K | kJ mol⁻¹ | fit through the 25 °C point. The 40 °C point gives 1808 K (weaker) | KOL20 [RUN] | M |
| Oxford Kokam 740 mAh pouch, Artemis drive-cycle discharge at 40 °C | — | net 0.483 A h in 3144 s; I_mean 0.55 A, I_rms 1.21 A, peak 5.0 A (6.8C) discharge and 1.6 A regen | — | — | 1 s logging. V fell from 4.18 to 3.63 V, cell 39.9–41.1 °C | OX17 [RUN] | H |
| Bonn iPhone 14 Pro initial battery capacity (as reported in the paper) | Q_14P | 3321–3326 | — | mAh | measured or reported by the authors' setup | SCH26 | M |
| Bonn iPhone 13 Pro initial battery capacity (pre-used units) | Q_13P | 3142–3181 | up to 3.05 % decrease during the study | mAh | as above | SCH26 | M |

### a.2 Device-level validation quantities (phone as a constant-power load)

E is the rated pack energy: iPhone 14 Pro 12.38 Wh, iPhone 13 Pro 11.97 Wh [N9M22], iPhone 14 Pro Max
16.68 Wh [NBC22]. For Android phones we use the published mAh × 3.87 V. That voltage is the nominal value implied
by Apple's Wh/mAh pairs: 12.38/3.200 = 3.869 V, 11.97/3.095 = 3.867 V and 16.68/4.323 = 3.858 V (DERIVED, M).
The P values below are DERIVED (`device_benchmarks.json`). The rated-energy uncertainty is about ±4 %.

| Name | Symbol | Value | Range | Unit | Conditions | Source | Conf. |
|---|---|---|---|---|---|---|---|
| Black-screen idle TTE, iPhone 14 Pro | TTE_idle | 55.0 | 54.4–57.3 (6 devices) | h | iOS 16.6, OLED on and black, airplane mode, Wi-Fi/BT off, auto-brightness off, 9 s logging, in a dark box | BONN26 [RUN]; SCH26 Table 6 (55:01:32) | H |
| Black-screen idle TTE, iPhone 13 Pro | TTE_idle | 50.8 | 50.6–51.4 | h | same; the units were pre-used | BONN26 [RUN]; SCH26 (50:48:53) | H |
| Implied idle power, 14 Pro / 13 Pro | P_idle | 0.225 / 0.236 | 0.21–0.24 | W | E/TTE | DERIVED | M |
| GPS "best" accuracy TTE, 14 Pro / 13 Pro | TTE_GPS | 48.5 / 29.8 | 45.1–51.3 / 28.9–31.2 | h | same box, 9 s fixes requested | BONN26 [RUN] | H |
| GPS increment over idle, 14 Pro / 13 Pro | ΔP_GPS | 0.030 / 0.167 | — | W | DERIVED: E(1/TTE_GPS − 1/TTE_idle) | DERIVED | M |
| Rear camera 30 fps TTE (frames not rendered), 14 Pro / 13 Pro | TTE_cam | 8.43 / 9.53 | 8.34–8.61 / 9.23–9.68 | h | 30, 60, 120 and 240 fps all within ±0.2 h (SCH26 Table 5) | BONN26 [RUN] | H |
| Camera increment over idle, 60 fps | ΔP_cam | 1.22 / 1.04 | — | W | DERIVED | DERIVED | M |
| Camera + GPS TTE, 14 Pro / 13 Pro | TTE_cxl | 8.38 / 8.33 | 8.11–8.53 / 8.04–8.42 | h | camera 60 fps + GPS best | BONN26 [RUN] | H |
| Additivity residual, 1/TTE_cxl vs 1/TTE_cam60 + 1/TTE_GPS − 1/TTE_idle | ε_add | +0.15 % (14 Pro), −0.10 % (13 Pro) | — | — | energy-free test | BONN26 [RUN] | H |
| Run time at a displayed "100 %" before it first shows 99 % | t_hold | 3.2 (idle, 14 Pro); 0.55 (camera) | — | h | ≈ 6–7 "percent-equivalents" of the steady drain rate | BONN26 [RUN] | H |
| Steady displayed-SOC drain rate, idle 14 Pro | dSOC/dt | 31.1 (99–50 %), 30.8 (50–20 %), 32.4 (20–1 %) | — | min per % | linear to within 6 % (GPS and camera runs similar) | BONN26 [RUN] | H |
| Audio playback, display not visible (Apple lab) | P_audio | 0.165 (14 Pro), 0.160 (13 Pro) | — | W | 75 h claim; BT headphones, Wi-Fi associated | APL14P, APL13P, APLBAT; DERIVED | M |
| Local video playback (Apple lab) | P_video | 0.54 | 0.54 (14 Pro 23 h; 13 Pro 22 h) | W | looped 2 h 23 min movie, auto-brightness off, brightness not disclosed | APL14P, APL13P; DERIVED | M |
| Streamed video playback (Apple lab) | P_stream | 0.62 / 0.60 | — | W | 20 h claim; looped HDR movie over Wi-Fi | APL14P, APL13P; DERIVED | M |
| Streaming network overhead (streamed − local) | ΔP_net,stream | 0.06–0.08 | — | W | DERIVED from the two rows above | DERIVED | L |
| Wi-Fi web browsing at 150 cd/m² | P_web,WiFi | 1.06 (14 Pro Max, S24); 1.31 (Pixel 8) | 1.05–1.31 | W | Notebookcheck scripted browsing | NBC22, NBC23, NBC24; DERIVED | M |
| Cellular web browsing at 150 nits | P_web,cell | 1.19 (14 Pro), 1.00 (13 Pro) | — | W | Laptop Mag, "web surfing at 150 nits on a mobile network", T-Mobile eSIM, mean of 3 runs | LMAG22; DERIVED | M |
| Local H.264 1080p video (review lab) | P_video,NBC | 0.59 (14 Pro Max), 0.78 (Pixel 8), 0.76 (S24) | — | W | Big Buck Bunny loop | NBC22–24; DERIVED | M |
| Screen-on idle at minimum brightness, WLAN off | P_idle,min | 0.29 (14 Pro Max), 0.39 (S24), 0.71 (Pixel 8) | — | W | Notebookcheck "Reader/Idle" | NBC22–24; DERIVED | M |
| Sustained max load (3D/"Burnout" benchmark) at max brightness | P_load | 2.4 (14 Pro Max), 5.2 (Pixel 8), 5.8 (S24 Exynos) | — | W | gaming-like upper bound; thermally limited | NBC22–24; DERIVED | M |
| Typical mixed use, median-user profile | P_mixed | 0.56 | — | W | Pixel 8a: "average battery life during testing was at least 31 hours", 4492 mAh typical | GPX; DERIVED | M |
| Standby (Apple claim, iPhone 7) | TTE_sb | 240 (up to 10 days) | — | h | i.e. P_sb ≈ 0.42 % of E per hour. iPhone 7 Wh not verified | APL7 | M |

### a.3 Dataset inventory, licences and reachability tests

Reachability was tested from this container on 2026-11-17 with `curl -I`, the API, or a download.

| # | Dataset (key) | What it contains | Dynamic loads? / temperatures | Licence (as verified) | Size | Reachability test result |
|---|---|---|---|---|---|---|
| 1 | NASA PCoE "Battery Data Set" [NASA07] | 34 × 2 Ah 18650 cells (B0005–B0056). Charge, discharge and EIS cycles: t, V, I, T_cell, capacity, Re, Rct | CC 1/2/4 A; **0.05 Hz 4 A square wave** (B0025–28); **4, 24, 43/44 °C**; cycled to 20–30 % fade | data.nasa.gov: "License not specified", public access. The PCoE page asks for acknowledgement (no licence text) | 209,708,670 B zip (6 nested zips, .mat) | **HEAD 200 and full download OK**: `https://phm-datasets.s3.amazonaws.com/NASA/5.+Battery+Data+Set.zip` (md5 7fbc095e…3fc7) |
| 2 | NASA "Randomized Battery Usage" [NASARW14] | 28 × LCO 18650 (~2.2 Ah), random-walk current profiles with periodic reference cycles | **yes, random-walk loads**; room temperature and 40 °C (per DOSREIS21) | as row 1 | **1,065,821,095 B** | HEAD 200. **Not downloaded** (over 300 MB) |
| 3 | NASA "Accelerated Battery Life Testing" [NASAALT23] | 26 two-cell 18650 packs under constant and variable loads, plus second-life packs | yes (variable load levels) | catalog.data.gov lists `usa.gov/government-works` | not stated | `data.nasa.gov/docs/legacy/battery_alt_dataset.zip`: 302 to a presigned S3 URL (reachable). Not downloaded |
| 4 | LG HG2 (Kollmeyer et al.) [KOL20] | 1 × 3 Ah NMC 18650: C/20 OCV, 0.5C, 2C, HPPC, pause cycles, UDDS/HWFET/LA92/US06 and "Mixed" drive cycles. Digatron logs of t, V, I, T, Ah, Wh | **yes, drive cycles**; **−20, −10, 0, 10, 25, 40 °C** | **CC BY 4.0** (Mendeley API) | 421 files, 817.8 MB total; single files 1 kB–6 MB | API 200. **18 CSVs downloaded (2.1 MB)**. Files 302 to S3 |
| 5 | Panasonic 18650PF (Kollmeyer) [KOL18] | 1 × 2.9 Ah NCA 18650: HPPC, EIS, US06/HWFET/UDDS/LA92/NN drive cycles | **yes**; 25, 10, 0, −10, −20 °C | **CC BY 4.0** (Mendeley API) | 259 files, 197.5 MB | API 200. File URL 302 to S3 (reachable). Not downloaded |
| 6 | Oxford Battery Degradation Dataset 1 [OX17] | 8 × Kokam 740 mAh pouch cells: 1C and pseudo-OCV (40 mA) characterisation every 100 drive cycles; example drive cycle at 1 Hz | **Artemis urban drive cycle**; 40 °C only | **ODbL v1.0 / DbCL v1.0** (readme and ORA page) | main .mat 253.8 MB; example 70.9 kB | Readme and example **downloaded**. Main file reachable but **~30 kB/s**, so the download was aborted at 5.7 MB |
| 7 | CALCE (U. Maryland) [CALCE] | CS2/CX2 **LCO prismatic** (1.1/1.35 Ah) cycling; PL **LCO pouch** 1.5 Ah partial cycling; INR 18650-20R and A123 with DST/FUDS/US06/BJDST | **yes**; 18650-20R at 0/25/45 °C; A123 −10…50 °C; CX2 one cell cycled at 25–55 °C | "open access"; cite the CALCE articles; **no licence stated** | e.g. SP2_25C_DST.zip 3,189,603 B; SP2_45C_US06.zip 3,945,329 B; A123_DST-US06-FUDS-N10.zip 9,213,655 B | **HEAD 200** for all 4 files tested |
| 8 | Severson et al. 2019 / TRI [SEV19] | 124 A123 LFP 1.1 Ah cells, fast-charge cycling to EOL at 30 °C | CC discharge only; one temperature | "CC BY 4.0" per DOSREIS21 (not seen on the JS-only site) | not determined | `data.matr.io/1/` HTML 200 (JS app). Files not tested |
| 9 | Sandia (SNL) via batteryarchive.org [PREGER20, BA] | 86 × 18650 LFP/NCA/NMC cells; 15/25/35 °C; DOD and C-rate matrix | CC cycling; 3 temperatures | batteryarchive: no site licence; data "uploaded with permission"; cite the article. "Double attribution" per DOSREIS21 | not determined | Pages HTML 200. The home page says full CSVs are sent **on request by email** |
| 10 | HNEI via batteryarchive [BA, DOSREIS21] | 15 × LG ICR18650 C2 (NMC-LCO blend) cells, 1.5C cycling to ~1000 cycles | no; 25 °C | "CC BY 4.0 + source attribution" per DOSREIS21 | — | not tested |
| 11 | **Bonn iPhone sensor battery dataset** [BONN26] | 29 experiments × 9 iPhones (3 × 13 Pro, 6 × 14 Pro), iOS 16.6. JSON of `batteryLevel` (0–1, 1 % steps) vs Unix time, plus sensor values and a config file | scenario-based (idle, GPS ×4 accuracies, camera ×4 fps, accelerometer, gyroscope, pedometer, TrueDepth, camera + GPS); room conditions | **CC BY 4.0** (Dataverse API and readme) | 30 files, 112.5 MB | API 200. **14 zips downloaded (27 MB)**: Idle, GPS, CAM, CXL |
| 12 | AndroWatts [AW24] | 1000 tests on one Android phone (model not verified). Per-rail energy (display, CPU clusters, GPU, Wi-Fi/BT, cellular, GPS…), screen RGB and brightness, battery discharge, SoC temperature | stimulus app, not usage traces | **CC BY 4.0** (Zenodo API) | open_data.zip 192,671,263 B | HEAD 200. **aggregated.csv (0.8 MB) extracted by HTTP range read** (0.84 MB transferred). Units not verified |
| 13 | GreenHub / BatteryHub [GH19, GH21] | 23.6 M crowdsourced Android samples (July 2019). One sample per 1 % change: battery level, state, voltage, temperature, screen on, brightness, network, CPU, apps | real-world, uncontrolled | **Data licence not verified.** The app code is Apache-2.0. The dataset is now on Kaggle (`hmatalonga/greenhub-farmer`) | not determined | greenhubproject.org 200. **farmer.greenhubproject.org: connection reset.** Kaggle page 404 to curl (needs JS/login) |
| 14 | Carat top-1000 users [CARAT13] | 2014–2018 app usage plus battery level and status, 18.1 M records | real-world | "research purposes only. Commercial use is strictly prohibited". Password-protected zip | **6,375,824,106 B** | HEAD 200. Not downloaded (size; not an open licence) |
| 15 | Device Analyzer (Cambridge) [DA14] | 16 k+ Android users, usage incl. charging | real-world | access by contacting the team (not openly licensed) | — | not tested |

### a.4 Published device-level benchmarks (cite as published measurements or specifications)

Implied power uses E from a.2. "Brightness" is the stated test brightness. Apple does not disclose the brightness
used in its video tests; only "Auto-Brightness and True Tone off" is stated [APLBAT].

| # | Device (display) | Battery | Test type | Conditions | Reported runtime | Implied P_avg | Source |
|---|---|---|---|---|---|---|---|
| 1 | iPhone 14 Pro (6.1" OLED, Always-On, 1000 nits max typ.) | 12.38 Wh / 3200 mAh (rated) [N9M22] | video (local) | Apple lab, Jul–Aug 2022, preproduction units, Wi-Fi associated, BT headphones | up to 23 h | 0.54 W | APL14P, APLBAT |
| 2 | iPhone 14 Pro | same | video (streamed) | repeated 3 h 01 min HDR movie | up to 20 h | 0.62 W | APL14P |
| 3 | iPhone 14 Pro | same | audio | 358 AAC tracks, display not visible (AOD on) | up to 75 h | 0.165 W | APL14P, APLBAT |
| 4 | iPhone 13 Pro (6.1" OLED) | 11.97 Wh / 3095 mAh [N9M22] | video / streamed / audio | Apple lab, Aug 2021 | 22 / 20 / 75 h | 0.54 / 0.60 / 0.16 W | APL13P |
| 5 | iPhone 15 Pro (6.1" OLED) | Wh not verified | video / streamed / audio | Apple lab, 2023 | 23 / 20 / 75 h | — | APL15P |
| 6 | iPhone 16 Pro (6.3" OLED) | Wh not verified | video / streamed / audio | Apple lab, Jul–Aug 2024 | 27 / 22 / 85 h | — | APL16P |
| 7 | iPhone 14 Pro | 12.38 Wh | web browsing, cellular | Laptop Mag: continuous browsing at 150 nits, T-Mobile eSIM, mean of 3 runs | 10 h 22 min | 1.19 W | LMAG22 |
| 8 | iPhone 13 Pro | 11.97 Wh | web browsing, cellular | same | 12 h 00 min | 1.00 W | LMAG22 |
| 9 | iPhone 14 Pro Max (6.7" OLED) | 16.68 Wh / 4323 mAh | Wi-Fi browsing | Notebookcheck, fixed 150 cd/m² | 15 h 48 min | 1.06 W | NBC22 |
| 10 | iPhone 14 Pro Max | same | video H.264 1080p / idle (min brightness, no WLAN) / max load | Notebookcheck | 28 h 08 min / 57 h 53 min / 6 h 56 min | 0.59 / 0.29 / 2.41 W | NBC22 |
| 11 | Google Pixel 8 (6.2") | 4575 mAh typ. / 4485 min [GPX] | Wi-Fi browsing / video / idle / load | Notebookcheck WiFi v1.3 (150 nits script) / Big Buck Bunny / Reader / Burnout at max brightness | 13 h 34 min / 22 h 43 min / 24 h 53 min / 3 h 25 min | 1.31 / 0.78 / 0.71 / 5.2 W | NBC23 |
| 12 | Samsung Galaxy S24, Exynos 2400 (6.2", 2340×1080) | 4000 mAh (typ.) | Wi-Fi browsing / video / idle / load | Notebookcheck: Chrome 122, 150 cd/m²; Big Buck Bunny; Reader; Burnout | 14 h 40 min / 20 h 27 min / 39 h 58 min / 2 h 39 min | 1.06 / 0.76 / 0.39 / 5.8 W | NBC24 |
| 13 | Google Pixel 8a (6.1") | 4492 mAh typ. / 4404 min | typical mixed use | Google: median-user profile (talk/data/standby/other), US carrier, California, late 2023 | "average … at least 31 hours" (claim "24+ h") | 0.56 W | GPX |
| 14 | Google Pixel 8 / 8a / 7a | as above | extreme battery saver | Google: same profile with Extreme Battery Saver (5G off, etc.) | up to 72 h | ≈ 0.24 W (Pixel 8a) | GPX |
| 15 | iPhone 14 Pro / 13 Pro | as row 1 / row 4 | black-screen idle; GPS best; rear camera; camera + GPS | Bonn lab, iOS 16.6, airplane mode, OLED black (section a.2) | 55.0 / 48.5 / 8.4 / 8.4 h (14 Pro); 50.8 / 29.8 / 9.5 / 8.3 h (13 Pro) | 0.225 / 0.255 / 1.47 / 1.48 W (14 Pro) | BONN26 [RUN] |
| 16 | iPhone 7 (4.7" LCD) | Wh not verified | standby / Wi-Fi internet / video / audio / talk (3G) | Apple lab, Aug 2016 | up to 10 days / 14 h / 13 h / 40 h / 14 h | — | APL7, APLBAT |

Note on the Notebookcheck "WiFi v1.3" test: a browser script opens a rotating set of websites (some with video)
with the display dimmed to 150 nits. This description comes from a Notebookcheck 2026 laptop comparison and was
seen only as a search snippet (L). Notebookcheck's own power readings (e.g. iPhone 14 Pro Max idle avg 1.7 W,
load avg 3.85 W) are taken at the charger. Do **not** compare them with battery-side powers.

---

## (b) Key equations and model forms used for validation

1. **Runtime to mean power (energy bucket, constant-power load).** With usable energy E_u:
   `TTE = E_u / P̄` ⇒ `P̄ = E_u / TTE`.
   For a published runtime we use E_u ≈ E_rated = Q_rated·V_nom, with V_nom ≈ 3.87 V (from a.2). The error in
   E is about ±4 %.

2. **Energy-free component increments.** Within one device, for two scenarios:
   `ΔP_x = P_x − P_idle = E_u (1/TTE_x − 1/TTE_idle)`, and `TTE_x / TTE_idle = P_idle / P_x` (E cancels).

3. **Additivity test.** If component powers add and E_u does not depend on load:
   `1/TTE_{a+b} = 1/TTE_a + 1/TTE_b − 1/TTE_idle`.
   Bonn data: the residual is +0.15 % (14 Pro) and −0.10 % (13 Pro) for camera 60 fps + GPS [RUN].

4. **Displayed SOC vs model SOC (iOS).** The phone does not report true SOC; it reports a fuel-gauge value with a
   top hold. A minimal map consistent with the Bonn traces:
   `SOC_disp = min(1, SOC_model / (1 − h))`, with h ≈ 0.06–0.07, rounded down to 1 % steps (L).
   Validate on TTE and on the slope over 99→20 %, not on the first or last percent.

5. **Temperature scaling of usable capacity and resistance** (to check our Arrhenius-ECM against HG2):
   `R(T) = R(T_ref) · exp[(E_a/R_g)(1/T − 1/T_ref)]`, with E_a/R_g ≈ 3170 K for the *apparent* resistance
   (E_a ≈ 26 kJ/mol) [RUN].
   Capacity ratios f_Q(T) at 0.5C and 2C (a.1) are the targets the ECM must reproduce from
   `V_term = OCV(z) − I·R(T) − V_RC` hitting the cutoff earlier.

6. **Coulomb counting on cell data** (for computing Q and E from the logs):
   `Q = ∫ I dt / 3600`, `E = ∫ I V dt / 3600` (trapezoid rule) [RUN].

7. **Validation metrics** (proposed):
   - Relative TTE error: `ε_TTE = (TTE_model − TTE_obs)/TTE_obs`.
   - Voltage RMSE over a discharge: `RMSE_V = sqrt(mean((V_model − V_obs)²))`.
   - Capacity-ratio error: `|f_Q,model(T) − f_Q,obs(T)|`.
   - Coverage of the predicted uncertainty band: the fraction of the 6 (14 Pro) or 3 (13 Pro) Bonn device
     runtimes inside the model's 90 % interval.

---

## (c) Key findings

1. **The Bonn dataset is the one open phone-level data source that suits a physics model.** It is CC BY 4.0,
   well documented, and has controlled settings (Table 1 of SCH26). It provides full 100→0 % discharges of
   modern phones with repeated devices, so cross-device spread can be quantified: CV is 0.4–4.9 % depending
   on the scenario. Each run's mean power is close to constant, which suits a constant-power battery-side
   load in our ODE.
2. **Two linearity checks pass.** Component powers add, and usable energy does not depend on load from 0.23 W
   to 1.48 W. The camera + GPS runtime is predicted from the single-sensor runtimes to within 0.15 % [RUN].
   This supports building P(t) as a sum of components. It also says that at phone-level C-rates
   (≈ 0.02C–0.12C for a ~12 Wh pack) and room temperature, rate-dependent capacity loss is negligible.
   This matches the HG2 cell: at 25 °C, 2C still delivers 0.97 × the 0.5C capacity.
3. **GPS matters surprisingly little.** On the iPhone 14 Pro, GPS at "best" accuracy adds only ≈ 0.03 W, or
   13 % over the black-screen idle. On the older 13 Pro it adds ≈ 0.17 W, a 1.7× power ratio. The camera
   pipeline is the dominant sensor cost at ≈ 1.0–1.2 W [RUN]. Frame rate (30–240 fps) makes almost no
   difference to camera runtime (SCH26 Table 5).
4. **Rough activity ranking from published runtimes** (battery-side P_avg, DERIVED):
   | Activity | P_avg |
   |---|---|
   | standby | ≈ 0.4 % of E per hour |
   | audio, screen off | 0.16 W |
   | black-screen idle | 0.23 W |
   | screen-on idle at minimum brightness | 0.3–0.7 W |
   | local video | 0.54–0.78 W |
   | streamed video | 0.60–0.62 W |
   | typical mixed use | ≈ 0.56 W |
   | web browsing at 150 nits | 1.0–1.3 W |
   | camera | 1.3–1.5 W |
   | sustained max load or gaming at max brightness | 2.4–5.8 W, thermally limited |

   Our component model (device_power.md) should reproduce each row within its spread. The order is robust
   across vendors, but absolute values differ by up to 2.4× between phones for the same nominal test (idle
   at minimum brightness: 0.29 W on the iPhone 14 Pro Max vs 0.71 W on the Pixel 8).
5. **Cold hurts mainly through kinetics.** On the HG2 cell, the C/20 capacity at −20 °C is still 0.90 × the
   25 °C value. At 0.5C it is 0.73 × and at 2C 0.49 ×, and energy falls further because of the lower mean
   voltage. So a phone in the cold loses runtime through I·R(T), and more so at high power. An
   Arrhenius-ECM with E_a ≈ 26 kJ/mol for the lumped resistance reproduces the trend [RUN].
6. **The NASA 18650 data cross-checks the temperature trend.** At ~2 A they give a capacity ratio of ≈ 0.71
   between 4 °C and 24 °C, which is a stronger cold penalty than HG2. Those cells are older and of a
   different chemistry. The NASA data also give capacity fade versus cycles (24–41 % over 168 cycles at 1C)
   and only modest growth of the EIS-fitted resistance for most cells. That is useful for the SOH part, but
   it is an accelerated 1C/100 % DOD protocol, much harsher than phone use.
7. **No open dataset has phone-level runs at controlled cold or hot ambient.** For temperature, the plan
   must combine cell-level data (HG2, NASA, Panasonic) with the device-level power model.
8. **Crowd datasets (GreenHub, Carat, Device Analyzer) are real-world, but none has a verified open licence**
   (Carat is explicitly research-only). GreenHub's data host was unreachable from this container. They can
   be cited as related work, not used as data under the problem's open-licence rule.

### Recommended validation plan

The parameters come from device_power.md and cell_ecm.md. The data below are used **only** for validation, never
for fitting, unless marked "calibration".

1. **Phone-level TTE, primary.** Source: Bonn [BONN26], CC BY 4.0.
   1. Build iPhone 14 Pro and 13 Pro configurations: E = 12.38 / 11.97 Wh ± 4 %. The screen is on and black,
      so only P_d0 of the display counts; radios are off; GPS and camera use the device_power.md defaults.
   2. Predict TTE for idle, GPS-best and camera-60 fps. Compare with the medians and min–max ranges in a.2.
      Target |ε_TTE| ≤ 15 %, and the observed device spread should fall inside the predicted 90 % interval.
   3. **Holdout test:** predict camera + GPS from the model alone and compare with 8.38 h / 8.33 h.
   4. Plot model SOC (with the top-hold map from b.4) against the median curves in
      `bonn_iphone_summary_soc_curves.csv`. Report the slope error over 99→20 %.
   5. Optional calibration variant: fit P_idle from idle runs only, then predict the other scenarios. This
      shows how far one baseline calibration propagates.
2. **Cell-level temperature, primary.** Source: LG HG2 [KOL20], CC BY 4.0.
   1. Using our ECM (OCV from cell_ecm.md, or recalibrated on the 25 °C C/20 curve, "calibration"), predict
      the 0.5C and 2C discharges at −20 … 40 °C.
   2. Compare Q, E, V̄ and the V(t) RMSE. Targets: |Δf_Q| ≤ 0.05 for T ≥ 0 °C; reproduce the −20 °C
      collapse (0.73 at 0.5C, 0.49 at 2C) within ±0.1; RMSE_V ≤ 50 mV for T ≥ 0 °C.
   3. Dynamic check: run one UDDS or US06 file at 25 °C and one at 0 °C (1.6–1.8 MB each in the same dataset).
3. **Aging and SOH, secondary.** Source: NASA PCoE [NASA07].
   1. Use B0005/6/7/18 at 24 °C: capacity vs cycle, Re/Rct growth, and the decimated fresh and aged B0005
      curves in `data/validation/`.
   2. Check that our SOH model (capacity fade plus resistance growth) gives the observed fresh-to-aged change
      in TTE at 2 A: 55.8 → 39.7 min, from `nasa_curves_index.json`.
   3. Use B0029 (43 °C, 4 A) and B0047/B0055 (4 °C) as temperature sanity checks.
4. **Dynamic load on a small pouch cell, optional.** Source: Oxford ExampleDC [OX17], ODbL.
   1. Drive the ECM with the 1 Hz Artemis current trace at 40 °C and compare V(t) (RMSE).
   2. This tests the RC dynamics under fast pulsed loads, which phones also produce (radio bursts, CPU
      boost).
5. **Cross-device sanity band, cited specifications.**
   1. Compute model TTE for the a.4 scenarios (audio, video, browsing at 150 nits, max load) on a
      generic 12–18 Wh phone.
   2. Report a table of model vs published values with ±25 % tolerance, given the undisclosed brightness and
      the vendor differences.
   3. Present this as plausibility evidence, not as validation on open data.

---

## (d) Caveats

* **Licences.**
  * The NASA PCoE set carries no explicit licence. data.nasa.gov says "License not specified", and NASA only
    asks for acknowledgement. Treating it as a US-Government work is background knowledge and was
    **not verified** here.
  * The CALCE data are "open access" with a citation requirement but **no licence**.
  * The TRI/Severson "CC BY 4.0" and Sandia "double attribution" labels are from the review DOSREIS21, not
    read on the hosting sites (M).
  * Manufacturer and review-site runtimes are copyrighted publications. We use the numbers as cited facts,
    which the problem statement allows for published specifications with citation. We do not redistribute
    them as a dataset.
* **Bonn dataset.**
  * Device model is not stored in the files. The 13 Pro / 14 Pro split is **inferred** from idle runtimes
    (3 shortest = 13 Pro). It matches the paper's 3 + 6 split and Table 6 medians (our 55.03 h vs their
    55:01:32), and it stays consistent across all experiments.
  * Battery level is the iOS fuel-gauge value at 1 % resolution, not true SOC. Logging stops at 1 %.
  * The 13 Pro units were pre-used and their capacity drifted by up to 3 % during the study [SCH26].
  * The black screen means **display brightness was not exercised** (SCH26 lists brightness as future work).
    So the Bonn data validate the baseline and the sensors, not the display model.
  * The paper's §4.2.1 text says GPS-best runtime is "∼30 h and 20 h" for the new and old devices. The data
    give 48.5 h / 29.8 h, which agrees with the paper's own "≈ 30 min per percent (14 Pro) and 20 min (13 Pro)"
    in §4.2.3/3.6. We trust the data.
* **Rated energy is an approximation of usable energy.**
  * Rated Wh for iPhones comes from press reports of Apple's Chemtrec filing and a Chinese regulatory
    database [N9M22] (M).
  * Android Wh uses an assumed V_nom = 3.87 V.
  * Real usable energy differs from the rated value through the cutoff, aging, temperature and OS reserves.
    For a fresh device we use ±4 % in UQ (L).
* **Published runtimes are not standardised.**
  * Apple does not disclose its test brightness.
  * Notebookcheck and Laptop Mag use 150 nits, but over Wi-Fi vs cellular, on different networks.
  * Review units differ (e.g. Exynos vs Snapdragon S24).
  * Use these runtimes as a band, not as point targets.
* **Cell-data quirks.**
  * NASA 4 °C cells contain unexplained low-capacity runs (README). We filtered to runs that start at ≥ 4.1 V
    and last more than 5 min, and use the median of the first 5. Values at 4 °C are therefore M.
  * B0050/B0052 have corrupt late EIS fits (Re ≈ −5000 Ω). We ignored them.
  * The HG2 data come from a single cell tested sequentially, so there is a small aging confound (≈ 3 %).
  * The HG2 logger writes the first row about 60 s into each step. We fill that minute with the first logged
    values, an error of < 1 % of Q.
  * HG2 apparent resistance lumps ohmic, polarisation and OCV-window effects. Do not equate it with R0.
* **Chemistry mismatch.** None of the open cells is a modern phone cell (high-voltage LCO-based pouch;
  background, unverified). Closest are the Kokam pouch (Oxford) and the CALCE LCO cells. Validation of *trends* (temperature,
  rate, aging) transfers better than validation of absolute values.
* **Process notes.**
  * Listing the nested NASA zip by HTTP range reads moved ≈ 751 MB over the network before we switched to a
    single 210 MB download. Nothing larger than 300 MB was stored.
  * The Oxford main file (254 MB) was abandoned at ≈ 30 kB/s.
  * GSMArena battery tests were not accessible: a tollbit redirect blocks automated access.
  * GreenHub's Farmer host refused connections.

---

## (e) References

Licence and access notes are as verified on 2026-11-17 unless marked.

1. **[NASA07]** B. Saha and K. Goebel (2007). "Battery Data Set". NASA Prognostics Data Repository, NASA Ames
   Research Center, Moffett Field, CA.
   - Download: https://phm-datasets.s3.amazonaws.com/NASA/5.+Battery+Data+Set.zip
   - Repository: https://www.nasa.gov/intelligent-systems-division/discovery-and-systems-health/pcoe/pcoe-data-set-repository/
   - Catalogue: https://data.nasa.gov/dataset/li-ion-battery-aging-datasets
   - *Licence: not specified; public access; acknowledgement requested. README files read in this session.*
2. **[NASARW14]** B. Bole, C. Kulkarni, M. Daigle. "Randomized Battery Usage Data Set". NASA Prognostics Data
   Repository. https://phm-datasets.s3.amazonaws.com/NASA/11.+Randomized+Battery+Usage+Data+Set.zip.
   Associated paper: Bole, Kulkarni, Daigle, "Adaptation of an Electrochemistry-based Li-Ion Battery Model to
   Account for Deterioration Observed Under Randomized Use", PHM Society 2014.
   - *Licence as in [NASA07]; 1.07 GB; not downloaded.*
3. **[NASAALT23]** K. Fricke, R. Nascimento, M. Corbetta, C. Kulkarni, F. Viana. "Accelerated Battery Life
   Testing Dataset". NASA Prognostics Data Repository / UCF.
   - Catalogue (listed there as "Randomized and Recommissioned Battery Dataset", same description: 26 two-cell
     packs, UCF testbed): https://catalog.data.gov/dataset/randomized-and-recommissioned-battery-dataset
   - Download: https://data.nasa.gov/docs/legacy/battery_alt_dataset.zip
   - Related paper: Fricke et al. 2023, PHM Conf., doi:10.36001/phmconf.2023.v15i1.3463.
   - *Licence: U.S. Government Works (usa.gov/government-works).*
4. **[DOSREIS21]** G. dos Reis, C. Strange, M. Yadav, S. Li (2021). "Lithium-ion battery data and where to find
   it". *Energy and AI* 5, 100081. doi:10.1016/j.egyai.2021.100081.
   - Accepted manuscript: https://www.pure.ed.ac.uk/ws/files/206117798/Batteries_Li_ion_battery_data_and_where_to_find_it_final_draft_.pdf
   - *Read in this session: dataset descriptions and licence labels.*
5. **[BA]** BatteryArchive.org (Sandia National Laboratories). Study summaries.
   https://www.batteryarchive.org/study_summaries.html
   - *No site licence. Data uploaded with permission; cite the original article and batteryarchive.org.*
6. **[PREGER20]** Y. Preger et al. (2020). "Degradation of Commercial Lithium-Ion Cells as a Function of
   Chemistry and Cycling Conditions". *J. Electrochem. Soc.* 167, 120532. doi:10.1149/1945-7111/abae37.
   - *Data on batteryarchive.org.*
7. **[SEV19]** K. A. Severson et al. (2019). "Data-driven prediction of battery cycle life before capacity
   degradation". *Nature Energy* 4, 383–391. doi:10.1038/s41560-019-0356-8.
   - Data: https://data.matr.io/1/
   - *Licence "CC BY 4.0" per DOSREIS21, not seen on the site.*
8. **[CALCE]** CALCE Battery Research Group, University of Maryland. "Battery Data".
   https://calce.umd.edu/battery-data (files under https://web.calce.umd.edu/batteries/data/)
   - *Open access; cite the CALCE articles; no licence stated.*
9. **[OX17]** C. Birkl, D. Howey (2017). "Oxford Battery Degradation Dataset 1". University of Oxford.
   doi:10.5287/bodleian:KO2kdmYGg. https://ora.ox.ac.uk/objects/uuid:03ba4b01-cfed-46d3-9b1a-7d4a7bdf6fac
   - The authors also ask users to cite C. R. Birkl, "Diagnosis and Prognosis of Degradation in Lithium-Ion
     Batteries", PhD thesis, University of Oxford, 2017 (thesis not read).
   - *Licence: ODC ODbL v1.0; contents DbCL v1.0.*
10. **[KOL18]** P. Kollmeyer (2018). "Panasonic 18650PF Li-ion Battery Data". Mendeley Data v1.
    doi:10.17632/wykht8y7tg.1.
    - *Licence: CC BY 4.0.*
11. **[KOL20]** P. Kollmeyer, C. Vidal, M. Naguib, M. Skells (2020). "LG 18650HG2 Li-ion Battery Data and
    Example Deep Neural Network xEV SOC Estimator Script". Mendeley Data v3. doi:10.17632/cp3473x7xv.3.
    - *Licence: CC BY 4.0.*
12. **[BONN26]** F. Schweizer, J. Yu, E. Mille, L. M. Reimer, M. Kapsecker, J. Klinker, S. M. Jonas (2026).
    "Battery Impact of Individual Hardware Sensors in Modern iPhone Devices" [dataset]. bonndata.
    doi:10.60507/FK2/LSJGWW.
    - *Licence: CC BY 4.0.*
    - Test-app code: doi:10.5281/zenodo.18800055.
13. **[SCH26]** F. Schweizer, J. Yu, E. Mille, L. M. Reimer, M. Kapsecker, J. Klinker, S. Jonas (2026).
    "Smartphone Sensor Battery Consumption: A Standardized and Reproducible Test Protocol". *Sensors* 26(10),
    2923. doi:10.3390/s26102923.
    - *Open access; PDF read in this session.*
14. **[AW24]** E. Guegain, R. Raes, N. Chachignot, C. Quinton, R. Rouvoy (2024). "AndroWatts: Unpacking the
    Power Consumption of Mobile Device's Components" [dataset]. Zenodo. doi:10.5281/zenodo.14314943.
    - Paper at MOBILESoft 2025.
    - *Licence: CC BY 4.0.*
15. **[GH19]** H. Matalonga, B. Cabral, F. Castor, M. Couto, R. Pereira, S. Melo de Sousa, J. P. Fernandes
    (2019). "GreenHub Farmer: Real-world data for Android energy mining". MSR 2019, IEEE
    (ieeexplore.ieee.org/document/8816755).
    - *Data licence not verified.*
16. **[GH21]** R. Pereira, H. Matalonga, M. Couto, F. Castor, B. Cabral, P. Carvalho, S. Melo de Sousa,
    J. P. Fernandes (2021). "GreenHub: a large-scale collaborative dataset to battery consumption analysis of
    android devices". *Empirical Software Engineering* 26, 38. doi:10.1007/s10664-020-09925-5.
    - Dataset page: https://greenhubproject.org/dataset (points to Kaggle `hmatalonga/greenhub-farmer`).
17. **[CARAT13]** A. J. Oliner, A. P. Iyer, I. Stoica, E. Lagerspetz, S. Tarkoma (2013). "Carat: Collaborative
    Energy Diagnosis for Mobile Devices". ACM SenSys 2013.
    - Data: https://www.cs.helsinki.fi/group/carat/data-sharing
    - *Research use only, no commercial use.*
18. **[DA14]** D. T. Wagner, A. Rice, A. R. Beresford. "Device Analyzer: Understanding smartphone usage".
    MobiQuitous 2013, LNICST 131 (2014), doi:10.1007/978-3-319-11569-6_16.
    - *Data access by request.*
19. **[APLBAT]** Apple Inc. "iPhone — Battery Test Information / General Battery Information".
    https://www.apple.com/iphone/battery.html (accessed 2026-11-17).
20. **[APL14P]** Apple Inc. "iPhone 14 Pro — Technical Specifications". https://support.apple.com/en-us/111849
21. **[APL13P]** Apple Inc. "iPhone 13 Pro — Technical Specifications". https://support.apple.com/en-us/111871
22. **[APL15P]** Apple Inc. "iPhone 15 Pro — Tech Specs". https://support.apple.com/en-us/111829
23. **[APL16P]** Apple Inc. "iPhone 16 Pro — Tech Specs". https://support.apple.com/en-us/121031
    - iPhone 16 (22/18/80 h): https://support.apple.com/en-us/121029
    - iPhone 14 (20/16/80 h): https://support.apple.com/en-us/111850
24. **[APL7]** Apple Inc. "iPhone 7 — Technical Specifications". https://support.apple.com/en-us/111943
25. **[N9M22]** 9to5Mac (11 Sep 2022). Article on iPhone 14 battery sizes.
    https://9to5mac.com/2022/09/11/iphone-14-battery-sizes/
    - Secondary report of Apple's Chemtrec filing (Wh) and a Chinese regulatory database (mAh): 14 Pro
      12.38 Wh / 3200 mAh; 13 Pro 11.97 Wh / 3095 mAh; 14 12.68 Wh / 3279 mAh; 14 Pro Max 16.68 Wh / 4323 mAh.
26. **[GPX]** Google. "Pixel phone hardware tech specs". https://support.google.com/pixelphone/answer/7158570
    - Pixel 8: 4575 mAh typ.; Pixel 8a: 4492 mAh typ., "at least 31 hours"; Pixel 7a: 4385 mAh typ.;
      Pixel 9: 4700 mAh typ.
27. **[NBC22]** Notebookcheck (2022). "Apple iPhone 14 Pro Max review: A gigantic brawny smartphone".
    https://www.notebookcheck.net/Apple-iPhone-14-Pro-Max-review-A-gigantic-brawny-smartphone.659750.0.html
28. **[NBC23]** Notebookcheck (2023). "Google Pixel 8 smartphone review: Compact and with 7 years of updates".
    https://www.notebookcheck.net/Google-Pixel-8-smartphone-review-Compact-and-with-7-years-of-updates.768436.0.html
29. **[NBC24]** Notebookcheck (2024). "Samsung Galaxy S24 review – The best small Android smartphone gets
    cheaper but not better".
    https://www.notebookcheck.net/Samsung-Galaxy-S24-review-The-best-small-Android-smartphone-gets-cheaper-but-not-better.809474.0.html
30. **[LMAG22]** Laptop Mag (2022). "iPhone 14 Pro battery life results are in — and we're not impressed".
    https://www.laptopmag.com/news/iphone-14-pro-battery-life-results-are-in-and-were-not-impressed
31. **[ODBL]** Open Knowledge Foundation. Open Data Commons Open Database License v1.0.
    https://opendatacommons.org/licenses/odbl/1.0/

---

## Appendix: reproduction

Download each source into its own empty directory, then run the scripts from outside it with `python3 -I`.

```bash
# NASA (210 MB) -> unzip the nested zips -> summary + curves
python3 -I docs/research/scripts/validation_nasa_summary.py <extracted_dir> docs/research/data/validation/nasa_pcoe_summary.json
python3 -I docs/research/scripts/validation_nasa_extract_curves.py <extracted_dir> docs/research/data/validation
# LG HG2 CSVs (Mendeley API file list -> *_C20DisCh.csv, *_Dis_0p5C.csv, *_Dis_2C.csv)
python3 -I docs/research/scripts/validation_lg_hg2_summary.py docs/research/data/validation/lg_hg2_raw docs/research/data/validation/lg_hg2_summary.json
# Oxford example file
python3 -I docs/research/scripts/validation_oxford_example.py <ExampleDC_C1.mat> docs/research/data/validation
# Bonn iPhone zips (https://bonndata.uni-bonn.de/api/access/datafile/<id>) -> extract
python3 -I docs/research/scripts/validation_bonn_iphone_summary.py <unzipped_dir> docs/research/data/validation/bonn_iphone_summary.json
python3 -I docs/research/scripts/validation_device_benchmarks.py docs/research/data/validation/bonn_iphone_summary.json docs/research/data/validation/device_benchmarks.json
# Range-read one member of a large remote zip (used for AndroWatts aggregated.csv)
python3 -I docs/research/scripts/validation_remote_zip.py "https://zenodo.org/api/records/14314943/files/open_data.zip/content" list
```
