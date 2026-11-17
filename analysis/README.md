# analysis/: scripts → results and paper figures

Every script implements one entry of the analysis registry in `docs/design/FINAL_SPEC.md` §14. It imports the
frozen model package `src/phonebatt`, writes one JSON file to `results/`, and writes its figures to
`paper/figures/` (PDF, plus a PNG preview in `paper/figures/png/`) and its tables to `paper/tables/`.

## How to run

From the repository root:

```bash
PYTHONPATH=src:analysis python3 analysis/run_all.py            # all scripts in §14 order (missing ones skipped)
PYTHONPATH=src:analysis python3 analysis/run_all.py --list     # registry, with [missing] marks
PYTHONPATH=src:analysis python3 analysis/run_all.py --only run_power,AN2
PYTHONPATH=src:analysis python3 analysis/run_power.py          # one script
```

`run_all.py` runs each script in its own process, continues after a failure (`--stop-on-error` to stop), prints
a timing summary, writes `results/run_all_log.json`, and exits with 1 if any existing script failed.

## Registry (FINAL_SPEC §14)

| ID | Script | Results file (main keys) | Figures / tables |
|---|---|---|---|
| AN0 | `run_params.py` | `params.json`: `{symbol: {value, unit, dist, lo, mode, hi, source, conf, section}}`, `_derived`, `_usage_factors`, `_checks` | `tab_params.tex` (T1) |
| AN1 | `run_power.py` | `power.json`: `primitives.{ID}.{P_W, P_capped_W, components{}, groups{}, dPb_W}`, `v6.{case}`, `bands`, `proto_check` | `fig_power_breakdown.pdf`; `tab_scenarios.tex` (T2) |
| AN2 | `run_cell_maps.py` | `cell_maps.json`: `ocv_table`, `R0_grid_mOhm`, `R1_mOhm`, `tau1_s`, `Pcut_W`, `Pcut_sys_W`, `operating_points` | `fig_cell_maps.pdf` |
| AN3 | `run_scenarios.py` | `scenarios.json` | `fig_trajectories.pdf`, `fig_tte_s0.pdf`, `fig_throttle.pdf` |
| AN4 | `run_analytic.py` | `analytic.json`: `F1`, `F2_vs_L2_pct{P}`, `L1_vs_L0_pct`, `a_coeff{ohm, cut, total}`, `a_coeff_T_soh`, `F4_check`, `F4_map_check`, `safe_burst_W`, `F5_cap{Tamb}`, `Pcut_table`, `F2_steady_vs_L2_pct` | `tab_pcut.tex` (optional) |
| AN5 | `run_ctmc.py` | `ctmc.json` | `fig_ctmc_tte.pdf`, `fig_cv_scaling.pdf` |
| AN6 | `run_validation.py` | `validation.json` | `fig_validation_bonn.pdf`, `fig_validation_hg2.pdf`, `fig_validation_band.pdf` |
| AN7 | `run_uq.py` | `uq.json`: `"{S}\|s0=..\|Tamb=..[\|soh=..]".{median_h, q05_h, q95_h, cv, p_lt_target, s_off_q05_q50_q95, ...}`, `variance_split`, `delta_method_check`, `abc_bonn`, `step_check`, `deviations` | `fig_uq_violins.pdf`; `tab_tte_uq.tex` (T3) |
| AN8 | `run_sensitivity.py` | `sensitivity.json`: `elasticities."{S}\|{T}\|{SOH}".{param}`, `elasticity_units`, `proto_check`, `oat.{S}.{factor}.{lo_h, hi_h, lo_pct, hi_pct}`, `oat_ranking` | `fig_tornado.pdf`; `tab_elasticities.tex` |
| AN9 | `run_sobol.py` | `sobol.json` | `fig_sobol.pdf` |
| AN10 | `run_structural.py` | `structural.json` | `tab_structural.tex` (T4); `fig_structural.pdf` |
| AN11 | `run_usage_fluct.py` | `usage_fluct.json` | `fig_usage_fluct.pdf` |
| AN12 | `run_stranded.py` | `stranded.json` | `fig_stranded_map.pdf` |
| AN13 | `run_aging.py` | `aging.json`: `calibration`, `{A1..A7}.{t_yr[], soh[], r_ratio[], years_to_80, tte_S13_h[], tte_S12_h[], s_off_0C_pct[], s_off_m10C_pct[], ...}`, `soh_grid`, `winter_threshold`, `F11_check`, `beta_L_band`, `uq` | `fig_aging.pdf`; `tab_aging.tex` |
| AN14 | `run_recommendations.py` | `recommendations.json`: `base`, `levers.{id}.{s12_pct, ctmc_pct, *_ci90_pct, p_gain_*, rank_*, s12_L2_pct}`, `ranking`, `energy_share_S12`, `scenario_levers`, `os.{O1..O9}` | `fig_levers.pdf` |
| AN15 | `run_generalization.py` | `generalization.json`: `RP63.{Pi1..Pi5, Pi4_cold_aged}`, `classes.{name}.{Pi1..Pi5 {mid, lo, hi}, Pi4_cold_aged, inputs, notes}`, `scaling_rules` | `fig_generalization.pdf` |
| AN16 | `make_architecture.py` | – | `fig_architecture.pdf` |

The JSON keys of AN3 and AN5–AN15 are specified in FINAL_SPEC §14; see each script's docstring for extras.
Each analysis group also writes a short digest `results/<group>_summary.md` for the paper writer
(`basics_summary.md` covers AN0, AN1, AN2, AN4 and AN16; `uq_sens_summary.md` covers AN7 and AN8;
`aging_recs_summary.md` covers AN13, AN14 and AN15). `run_aging.py` (~5 min, `--workers 2`, `--n-uq 128`) and
`run_recommendations.py` (~1 min) accept `--figure-only` to redraw from their JSON; AN14 reads `aging.json`
(AN13) and `stranded.json` (AN12) for the O9/O1 cross-references when present.
`run_uq.py` caches its per-cell parameter samples in `results/cache/uq/*.npz` (git-ignored; `--fresh` recomputes;
`--workers N`, default 3, sets the process pool). A full run takes about 20–30 min of wall time on 3 cores.

## Shared modules

* `style.py`: figure style (read its docstring rules). Call `apply()` first and `save(fig, "fig_xxx")`.
* `basics_common.py`: helpers used by AN0, AN1, AN2, AN4: `meta()` (the top-level `"meta"` block: spec version,
  git SHA, whether `src/` differs from that commit, seed, solver, step), `write_json()` (NumPy-safe, NaN → null),
  `write_tex()`, `load_proto()` (imports `docs/design/scripts/final_proto.py` for regression comparisons).
* `structural_variants.py`: AN10's vectorised "union" model (core equations plus the extra terms of the
  structural variants X2–X14, X18; per-column switches, X0 = core) with an rk4_batch driver that follows
  `run_batch`'s step and event conventions; `summary_global.py`: writes `results/global_summary.md` (AN9 + AN10).
  `run_sobol.py --from-raw/--save-raw` re-analyses stored model outputs; `--figure-only` (both scripts) redraws
  from the JSON.
* `stochastic_common.py`: helpers of AN5/AN11: κ-scaled CTMCs, an event-driven usage simulator with lognormal
  sessions and per-visit jitter (exact L0 TTE, or `Paths` for `run_batch`), `ExpandedCTMC`, `l0_from_paths`,
  `summ()` (TTE moments/quantiles with standard errors). AN11 reuses `rate_scaling` from `results/ctmc.json` when
  AN5 produced it with the same configuration (otherwise it recomputes it).

## Conventions

* JSON numbers are SI unless the key says otherwise (`_h`, `_pct`, `_pp`, `_C`, `_W`, `_Wh`, `_mOhm`).
  Non-finite values are written as `null`.
* `[PROTO]` numbers from the spec are regression targets, not validation. Scripts that reproduce one store the
  target next to the result (`proto*` keys) and the deviation.
* Component colours in `fig_power_breakdown` (fixed categorical slots, following the entity): SoC = blue,
  display = orange, radio = aqua, floor = yellow, camera = magenta, GNSS + sensors = green, audio/BT = violet.
* Temperature lines in `fig_cell_maps` use the single-hue sequential blue (darker = colder).
