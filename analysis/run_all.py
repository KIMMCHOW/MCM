"""Run every analysis script of FINAL_SPEC §14 in registry order (AN0 ... AN16).

Usage (from the repository root)::

    PYTHONPATH=src:analysis python3 analysis/run_all.py               # everything, continue after failures
    PYTHONPATH=src:analysis python3 analysis/run_all.py --only run_power,run_cell_maps
    PYTHONPATH=src:analysis python3 analysis/run_all.py --from run_uq --stop-on-error
    PYTHONPATH=src:analysis python3 analysis/run_all.py --list

Each script runs in its own Python process (cwd = repository root, PYTHONPATH = src:analysis), writes one JSON
file to results/ and its figures/tables to paper/.  Missing scripts are skipped with a warning.  The exit code is
1 if any script that exists failed, 0 otherwise.  A timing summary is printed at the end and written to
results/run_all_log.json.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / "analysis"

# FINAL_SPEC §14 analysis registry, in order: (id, script, main outputs)
REGISTRY = [
    ("AN0", "run_params.py", "results/params.json; paper/tables/tab_params.tex"),
    ("AN1", "run_power.py", "results/power.json; fig_power_breakdown; tab_scenarios.tex"),
    ("AN2", "run_cell_maps.py", "results/cell_maps.json; fig_cell_maps"),
    ("AN3", "run_scenarios.py", "results/scenarios.json; fig_trajectories, fig_tte_s0, fig_throttle"),
    ("AN4", "run_analytic.py", "results/analytic.json"),
    ("AN5", "run_ctmc.py", "results/ctmc.json; fig_ctmc_tte, fig_cv_scaling"),
    ("AN6", "run_validation.py", "results/validation.json; fig_validation_bonn, fig_validation_hg2, "
                                 "fig_validation_band"),
    ("AN7", "run_uq.py", "results/uq.json; fig_uq_violins; tab_tte_uq.tex"),
    ("AN8", "run_sensitivity.py", "results/sensitivity.json; fig_tornado"),
    ("AN9", "run_sobol.py", "results/sobol.json; fig_sobol"),
    ("AN10", "run_structural.py", "results/structural.json; tab_structural.tex; fig_structural"),
    ("AN11", "run_usage_fluct.py", "results/usage_fluct.json; fig_usage_fluct"),
    ("AN12", "run_stranded.py", "results/stranded.json; fig_stranded_map"),
    ("AN13", "run_aging.py", "results/aging.json; fig_aging"),
    ("AN14", "run_recommendations.py", "results/recommendations.json; fig_levers"),
    ("AN15", "run_generalization.py", "results/generalization.json; fig_generalization"),
    ("AN16", "make_architecture.py", "fig_architecture"),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="comma-separated script names or AN ids to run (others skipped)")
    ap.add_argument("--from", dest="start", default="", help="start at this script name or AN id")
    ap.add_argument("--stop-on-error", action="store_true", help="stop at the first failing script")
    ap.add_argument("--list", action="store_true", help="list the registry and exit")
    ap.add_argument("--timeout", type=float, default=3600.0, help="per-script timeout in seconds (default 3600)")
    a = ap.parse_args(argv)

    def key(s: str) -> str:
        return s.strip().removesuffix(".py")

    if a.list:
        for an, script, outs in REGISTRY:
            mark = "" if (HERE / script).exists() else "   [missing]"
            print(f"{an:5s} {script:24s} {outs}{mark}")
        return 0
    only = {key(s) for s in a.only.split(",") if s.strip()}
    todo = REGISTRY
    if a.start:
        ks = [i for i, (an, sc, _) in enumerate(REGISTRY) if key(a.start) in (an, key(sc))]
        if not ks:
            print(f"unknown --from {a.start!r}", file=sys.stderr)
            return 2
        todo = REGISTRY[ks[0]:]
    if only:
        todo = [r for r in todo if r[0] in only or key(r[1]) in only]

    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(ROOT / "src"), str(HERE)] + ([env["PYTHONPATH"]] if env.get(
        "PYTHONPATH") else []))
    env.setdefault("MPLBACKEND", "Agg")
    log = []
    t_all = time.time()
    failed = False
    for an, script, outs in todo:
        path = HERE / script
        if not path.exists():
            print(f"WARNING: {an} {script} not found - skipped", flush=True)
            log.append(dict(id=an, script=script, status="missing", seconds=0.0))
            continue
        print(f"=== {an} {script} ===", flush=True)
        t0 = time.time()
        try:
            rc = subprocess.run([sys.executable, str(path)], cwd=ROOT, env=env, timeout=a.timeout).returncode
            status = "ok" if rc == 0 else f"failed (exit {rc})"
        except subprocess.TimeoutExpired:
            status = f"failed (timeout {a.timeout:.0f} s)"
        dt = time.time() - t0
        log.append(dict(id=an, script=script, status=status, seconds=round(dt, 1)))
        print(f"--- {an} {script}: {status} in {dt:.1f} s", flush=True)
        if status != "ok":
            failed = True
            if a.stop_on_error:
                break
    total = time.time() - t_all
    print("\nsummary")
    for e in log:
        print(f"  {e['id']:5s} {e['script']:24s} {e['status']:24s} {e['seconds']:8.1f} s")
    print(f"  total {total:.1f} s")
    try:
        (ROOT / "results").mkdir(exist_ok=True)
        with open(ROOT / "results" / "run_all_log.json", "w") as fh:
            json.dump(dict(meta=dict(spec="FINAL_SPEC 2026-11-17", script="analysis/run_all.py", seed=None,
                                     solver=None, h_s=None),
                           finished_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), total_s=round(total, 1),
                           runs=log), fh, indent=1)
    except OSError as e:
        print(f"could not write run_all_log.json: {e}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
