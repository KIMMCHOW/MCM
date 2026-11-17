# Core issues found by the validation group (AN6)

No bug in `src/phonebatt` was found. The local L2 cell driver (`analysis/validation_cells.py`) reproduces the core
LSODA path to 1e-6 (HG2 constant power, `phonebatt.run` with `UsageState(P_extra=P)`) and 2e-8 (Ai2020
configuration, isothermal). In the V1 UQ, `run_batch` at h = 10 s matches LSODA to <= 1.4e-6 on 64 samples x 4 cases,
with identical events. Two caveats matter for other groups.

## C1. RK4 stability with sampled tau1: section 6.3's "h <= 49 s" holds only for the default tau1 (spec caveat, not a code bug)

* The section 6.3 bound `|h·lambda| < 2.785` uses tau1 >= 17.6 s. `UQ_SPEC` samples `tau1_ref ~ U(10, 60)` s and
  `Ea_tau ~ U(15, 27)` kJ/mol. The fastest RC mode is then 1/10 s^-1 at 25 °C, which needs h < 27.9 s, and about
  1/5.3 s^-1 at 45 °C, which needs h < 15 s.
* `run_batch` has no guard against this.
* **Reproduction:** `sample(UQ_SPEC, 64, seed=2026)` mapped to the iPhone 14 Pro (`P_awake = 0.090`), Bonn idle case,
  `run_batch(..., h=30.0)`. Sample 1 ends after 0.14 h (`cutoff`, z = 0.997, T_end = 47.5 °C). LSODA gives 38.28 h.
  With h = 10 s all 64 samples agree with LSODA.
* **Recommendation:** in UQ, Sobol or sensitivity runs, keep h <= 10 s (h <= 15 s in hot scenarios). Optionally,
  add an assertion `h < 2.785·min(tau1(T_max))` in `run_batch`.

## C2. `params.nasa_b0005`: the B0005 runs end on the `charge` event, so V3 cannot see resistance growth (parameterisation limit)

* `Q_max` is set to the integrated 2 A capacity, and the Ai2020 4.2 V OCV ends at U(0) = 3.018 V. With
  R_e + R_ct = 0.114 Ω the 2 A drop is about 0.23 V, so V(z = 0) is about 2.79 V, above V_cut = 2.7 V. Every
  B0005 run therefore ends at z = 0, not at the cut-off.
* Effects on V3:
  * The aged/fresh TTE ratio reduces to the capacity ratio: model −28.7 % against −28.8 % observed. R growth alone
    gives −0.03 %.
  * The 4 °C and 43 °C sanity checks give a capacity ratio of exactly 1.00. Observed: 0.74 (B0055), 0.92 (B0047) and
    0.94 (B0029).
* A thermodynamic capacity (for example from a low-rate discharge) together with an OCV that extends below 3.0 V
  would be needed. Neither is in the extracted data.
* The numbers are reported as they are in `results/validation.json` (key `V3`).
