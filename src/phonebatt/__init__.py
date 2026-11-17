"""phonebatt: coupled electro-thermal-aging model of smartphone battery drain (2026 MCM Problem A).

Implements docs/design/FINAL_SPEC.md.  Main entry points::

    from phonebatt import RP63, scenario, run, run_batch, CTMC, sample
    r = run(scenario("S4").schedule, 25.0)            # L2 TTE, LSODA reference solver
    theta = sample(["R0_ref", "Ea_R0"], 1024)         # Sobol' parameter samples
    b = run_batch(scenario("S4").schedule, -10.0, theta, h=2.0)

Modules: ``params`` (§4), ``cell`` (§2.2–2.6, 2.9, 2.11), ``power`` (§2.1, 2.8, 3), ``usage`` (§1.4, 5),
``model`` (§2.13, 6), ``integrate`` (RK4 batch, segmented solve_ivp), ``analytic`` (§7), ``aging`` (§2.12, 11),
``validation`` (§8 loaders and cases).
"""
from . import aging, analytic, cell, integrate, model, params, power, usage, validation
from .aging import AGING_SCENARIOS, calibrate, day_profile, g_stress, run_years
from .analytic import (a_coefficient, brownout_root, ctmc_moments, hub_closed_form, ig_cdf, ig_quantiles,
                       lever_effects, mc_L0, safe_burst, semi_markov_sigma2, thermal_cap, tte_L0, tte_L1)
from .cell import (R0, R1, a_deg, a_R, current, display_soc, nominal_energy, ocv, peak_voltage, soc_to_z, tau1,
                   usable_energy, z_empty)
from .integrate import rk4_batch, solve_ivp_segmented
from .model import EVENTS, Result, Segment, event_batch, rhs_batch, run, run_batch
from .params import (RP63, UQ_SPEC, USAGE_UQ_SPEC, Params, aged, default, device, galaxy_s24, hg2_cell, iphone13pro,
                     iphone14pro, iphone14promax, nasa_b0005, pixel8, pixel8a, sample, with_theta)
from .power import power_parts, steady_power, throttle, x_eff
from .usage import CTMC, LEVERS, PRIMITIVES, Schedule, UsageState, scenario, s12_schedule, s9_schedule

SPEC_VERSION = "FINAL_SPEC 2026-11-17"
__version__ = "1.0.0"

__all__ = [
    "aging", "analytic", "cell", "integrate", "model", "params", "power", "usage", "validation",
    "AGING_SCENARIOS", "calibrate", "day_profile", "g_stress", "run_years",
    "a_coefficient", "brownout_root", "ctmc_moments", "hub_closed_form", "ig_cdf", "ig_quantiles", "lever_effects",
    "mc_L0", "safe_burst", "semi_markov_sigma2", "thermal_cap", "tte_L0", "tte_L1",
    "R0", "R1", "a_deg", "a_R", "current", "display_soc", "nominal_energy", "ocv", "peak_voltage", "soc_to_z",
    "tau1", "usable_energy", "z_empty", "rk4_batch", "solve_ivp_segmented",
    "EVENTS", "Result", "Segment", "event_batch", "rhs_batch", "run", "run_batch",
    "RP63", "UQ_SPEC", "USAGE_UQ_SPEC", "Params", "aged", "default", "device", "galaxy_s24", "hg2_cell",
    "iphone13pro", "iphone14pro", "iphone14promax", "nasa_b0005", "pixel8", "pixel8a", "sample", "with_theta",
    "power_parts", "steady_power", "throttle", "x_eff",
    "CTMC", "LEVERS", "PRIMITIVES", "Schedule", "UsageState", "scenario", "s12_schedule", "s9_schedule",
    "SPEC_VERSION",
]
