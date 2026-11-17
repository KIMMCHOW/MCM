# Core issues found by group "uq_sens" (AN7 run_uq.py, AN8 run_sensitivity.py)

**No bug found in `src/phonebatt`.** Every default-parameter value used by AN7/AN8 reproduces its [PROTO] target
(29 TTE cells within 0.08 %, the 36 §10.1 elasticities within 0.0005). One API hazard and two model properties
are recorded so that other groups do not mistake them for bugs.

## H1 (API hazard, not a bug): index-paired scrambled Sobol' samples are not independent

- Where: `src/phonebatt/params.py`, `unit_samples` / `sample` (method `"sobol"`).
- Each call is correct on its own. But two calls with different seeds, paired row by row, are deterministically
  coupled dimension by dimension: SciPy's scrambling is a bijection of the same radical inverse of the row
  index. So factor j of one sample is a function of factor j of the other.
- Minimal reproduction:

  ```python
  import numpy as np
  from phonebatt.params import unit_samples
  A = unit_samples(70, 2048, seed=2026); B = unit_samples(70, 2048, seed=2028)
  h, _, _ = np.histogram2d(A[:, 0], B[:, 0], bins=8, range=[[0, 1], [0, 1]])
  print(((h - 32) ** 2 / 32).sum())   # 14336 (maximal); independent pairs give about 63
  ```

- Seen in: AN7's first draft paired the hardware θ (`UQ_SPEC`, seed 2026) with the §4.6 usage-model factors
  (`USAGE_UQ_SPEC`, another Sobol' seed). This coupled Q_n with λ_on, R0_ref with d_game, and so on. A
  pick-freeze decomposition built on it gave a −23 % "interaction".
- Workaround in run_uq.py (exact): the second sample is pseudo-random (`method="random"`), or all factors come
  from one joint Sobol' sample.
- Suggestion for the core (not applied, since the core is frozen): one sentence in the `sample` docstring.
  The same applies to SALib/Saltelli designs only if A and B come from separate Sobol' calls.

## P1 (model property): s_off slightly below 0 for near-zero loads in the cold

- D1 is `s = (z − z_e)/(1 − z_e)` with `z_e` taken at T_ref (FINAL_SPEC §2.11, no clamp).
- At 0 / −10 °C the entropic term raises the OCV (∂U/∂T < 0). A 0.04 W standby load therefore reaches V_cut at
  z slightly below z_e. In uq.json the 5 % quantile of s_off for S1/S2 at 0 / −10 °C is −0.2 / −0.3 %.
- This matches the spec definition. Clip it for display only if the paper prefers.

## P2 (model property): cut-off vs brown-out labels near empty depend on the step

- When the cut-off and the brown-out envelope cross within one step of each other at the end of the
  discharge, the label in `classify` follows whichever is crossed first on the step grid. The TTE is the same.
- Measured: S4 and S7 at 25 °C, h = 2 s vs h = 10 s on 256 θ. TTE agrees to 2e-6 (relative). The label agrees
  on 98.8 % of θ: 3 of 256 flip between brown-out and cut-off.
- Read event shares at 25 °C together with s_off. For S4 at 25 °C, 65 % of θ end in a brown-out, but most of
  these end close to empty: s_off median 1.1 %, 95 % quantile 7.2 %. Only the h-dependent flips above are true
  ties.
