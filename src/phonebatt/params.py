"""Parameters of the reference phone RP-63 and their uncertainty distributions (FINAL_SPEC §4).

Design
------
``Params`` is a frozen dataclass holding every symbol of FINAL_SPEC §4.  Every field is a plain float by
default, but any numeric field may be replaced by a NumPy array of shape ``(n,)`` (one value per parameter
sample)::

    theta = sample(UQ_SPEC, 1024)            # dict name -> ndarray (n,)
    pb = with_theta(RP63, theta)             # Params whose fields are arrays
    ocv(z, T, pb)                            # broadcasts over the trailing sample axis

Units are SI except: capacities in A·h, temperatures of *limits* in °C (fields ending ``_C``), and
the display content density ``p_c`` in W/cm² (the spec quotes mW/cm²).  Activation energies are in J/mol.

Distributions (§4): ``U(a,b)`` uniform, ``T(a,m,b)`` triangular, ``LU(a,b)`` log-uniform, ``DIR`` Dirichlet.
``UQ_SPEC`` maps a field name to ``(dist, lo, mode, hi)``; ``mode`` is the default value.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, fields, replace
from typing import Any, Iterable, Mapping

import numpy as np

RG = 8.314  # J mol^-1 K^-1
T_REF = 298.15  # K
K0 = 273.15  # K offset

RATS = ("wifi", "lte", "nr")
RADIO_FIELDS = ("P_idle", "P_tail", "T_tail", "P_act", "T_act", "E_pro", "e_b", "P_rt")


@dataclass(frozen=True, eq=False)
class Params:
    """All model parameters (FINAL_SPEC §4).  Defaults are the RP-63 reference phone."""

    # ---- §4.1 cell and pack ----------------------------------------------------------------------
    Q_n: Any = 4.70  # A h, nominal (new) capacity
    b0: Any = 0.60117  # V    OCV coefficients (E5, CRM06 form fitted to Ai2020)
    b1: Any = 36.426  # -
    a0: Any = 3.6193  # V
    a1: Any = 0.44311  # V
    a2: Any = -0.48605  # V
    a3: Any = 0.63127  # V
    dV_hv: Any = 0.242  # V    high-voltage stretch amplitude
    p_hv: Any = 3.86  # -    high-voltage stretch exponent
    dU: Any = 0.0  # V    hysteresis / cell-to-cell OCV offset
    k_S: Any = 1.0  # -    multiplier on dU/dT
    R0_ref: Any = 0.060  # Ohm  pack R0 at z = 0.5, 25 °C
    A_R: Any = 0.80  # -    low-SOC rise of R0
    z0: Any = 0.139  # -    SOC scale of the R0 rise
    Ea_R0: Any = 18e3  # J/mol
    R1_ref: Any = 0.010  # Ohm
    tau1_ref: Any = 30.0  # s
    Ea_R1: Any = 30e3  # J/mol
    Ea_tau: Any = 21e3  # J/mol
    V_cut: Any = 3.40  # V    fuel-gauge cut-off
    d_sys: Any = 0.10  # V    V_sys = V_cut - d_sys (brown-out level)
    T_shut_C: Any = 68.0  # °C   battery thermal shutdown
    eta_c: Any = 0.99  # -    charge coulombic efficiency (aging layer only)

    # ---- §4.2 thermal -----------------------------------------------------------------------------
    C_th: Any = 150.0  # J/K
    R_th: Any = 3.0  # K/W
    R_skin: Any = 5.5  # K/W
    T_lim_C: Any = 45.0  # °C   skin limit
    dT_thr: Any = 2.0  # K    throttling starts dT_thr below T_lim
    w_thr: Any = 1.0  # K    logistic width
    P_fl: Any = 0.5  # W    throttled SoC floor
    k_env: Any = 1.0  # -    enclosure factor on R_th (and R_skin in F5)
    dT_ch: Any = 8.0  # K    cell temperature rise while charging (aging layer)
    eta_charger: Any = 0.90  # -

    # ---- §4.3 device power (battery side) -------------------------------------------------------
    P_susp: Any = 0.020  # W
    P_awake: Any = 0.090  # W   [CAL] (V1.1)
    E_wake: Any = 0.10  # J
    P_d0: Any = 0.12  # W
    p_c: Any = 0.25 / 97.4  # W/cm²  content power density at b = 0.5, APL_ref
    A_disp: Any = 97.4  # cm²  (6.3" 19.5:9 panel: 2.455 d²)
    APL_ref: Any = 0.70  # -
    gamma_b: Any = 2.0  # -
    k_RR: Any = 0.0028  # W/Hz
    P_SoC_max: Any = 5.0  # W
    v0: Any = 0.446  # -
    P_G: Any = 0.27  # W
    P_cam: Any = 1.15  # W
    P_camf: Any = 0.30  # W
    P_aud: Any = 0.05  # W
    P_BT: Any = 0.04  # W
    m_b: Any = 1.0  # -    common multiplier of burst-envelope amplitudes

    # ---- §4.4 radio (per RAT) ---------------------------------------------------------------------
    wifi_P_idle: Any = 0.004
    wifi_P_tail: Any = 0.12
    wifi_T_tail: Any = 0.24
    wifi_P_act: Any = 0.13
    wifi_T_act: Any = 0.24
    wifi_E_pro: Any = 0.0098
    wifi_e_b: Any = 0.015
    wifi_P_rt: Any = 0.30
    lte_P_idle: Any = 0.015
    lte_P_tail: Any = 0.12
    lte_T_tail: Any = 10.0
    lte_P_act: Any = 0.60
    lte_T_act: Any = 0.10
    lte_E_pro: Any = 0.31
    lte_e_b: Any = 0.030
    lte_P_rt: Any = 1.20
    nr_P_idle: Any = 0.030
    nr_P_tail: Any = 0.255
    nr_T_tail: Any = 10.3
    nr_P_act: Any = 0.90
    nr_T_act: Any = 0.10
    nr_E_pro: Any = 1.0
    nr_e_b: Any = 0.030
    nr_P_rt: Any = 1.40
    k_ul: Any = 5.0  # UL/DL energy per bit
    gamma_sig: Any = 0.26  # per 10 dB
    RSRP_ref: Any = -85.0  # dBm
    m_max: Any = 3.5  # -

    # ---- §4.5 aging ----------------------------------------------------------------------------
    k_c: Any = 0.027 / math.sqrt(365.0)  # d^-1/2
    E_cal: Any = 58e3  # J/mol
    g_hi: Any = 1.75  # high plateau of the KEIL16 SOC step
    k_HV: Any = 1.0  # exponent of the SCH14 high-voltage factor
    k_H: Any = 2.697e-3  # FCE^-1/2 (calibrated, §11.1)
    k_L: Any = 1.288e-3  # FCE^-1/2 (calibrated, §11.1)
    w_L: Any = 0.3  # plating share of cycle loss in the calibration
    E_H: Any = 33e3  # J/mol
    E_L: Any = 55.5e3  # J/mol
    beta_L: Any = 2.64  # h
    c_ref: Any = 1.0  # 1/h
    gamma_R: Any = 2.5  # -
    N_80: Any = 1000.0  # FCE to SOH 0.8 (calibration target)
    self_discharge_per_month: Any = 0.001  # storage scenarios only (A15)

    # ---- state of health within one discharge (§2.10) -----------------------------------------
    q_loss: Any = 0.0  # capacity loss fraction; Q_max = Q_n (1 - q_loss)
    r_growth: Any = None  # None -> gamma_R * q_loss (E20)

    # ---- switches (structural variants) ----------------------------------------------------------
    burst: bool = True  # burst envelope / brown-out event enabled (X7 sets False)
    throttling: bool = True  # skin-temperature throttling enabled
    load_heat: bool = True  # A4: the delivered power P is dissipated inside the enclosure (E7 term P).  False for
    # a bare cell on a tester (hg2_cell, nasa_b0005), whose load is external  [REVIEW 2026-11-17]
    ocv_table: Any = None  # optional (z_grid, U_grid) replacing the CRM06 OCV (e.g. HG2 cell, NMC)
    dudt_table: Any = None  # optional (z_grid, dUdT_grid [V/K]) replacing the Ai2020 table

    # ---- convenience --------------------------------------------------------------------------
    def replace(self, **kw) -> "Params":
        return replace(self, **kw)

    @property
    def V_sys(self):
        return self.V_cut - self.d_sys

    @property
    def T_lim(self):
        return self.T_lim_C + K0

    @property
    def T_shut(self):
        return self.T_shut_C + K0

    @property
    def rg(self):
        """Resistance growth fraction r_growth (E20)."""
        return self.gamma_R * self.q_loss if self.r_growth is None else self.r_growth

    @property
    def Q_max(self):
        return self.Q_n * (1.0 - self.q_loss)

    def radio(self, rat: str) -> dict:
        """Radio parameters of one RAT ('wifi', 'lte', 'nr') as a dict keyed by RADIO_FIELDS."""
        return {f: getattr(self, f"{rat}_{f}") for f in RADIO_FIELDS}

    def n_samples(self) -> int:
        """Length of the sample axis (1 if every field is scalar)."""
        n = 1
        for f in fields(self):
            v = getattr(self, f.name)
            if isinstance(v, np.ndarray) and v.ndim >= 1 and f.name not in ("ocv_table", "dudt_table"):
                n = max(n, v.shape[-1])
        return n


def default() -> Params:
    """The RP-63 reference phone (FINAL_SPEC §4)."""
    return Params()


RP63 = default()


def aged(p: Params, q_loss, r_growth=None) -> Params:
    """Return ``p`` at a fixed state of health: SOH = 1 - q_loss, r_growth = gamma_R*q_loss unless given."""
    return replace(p, q_loss=q_loss, r_growth=r_growth)


def with_theta(p: Params, theta: Mapping[str, Any] | None) -> Params:
    """Return ``p`` with the fields named in ``theta`` replaced (arrays allowed).  Unknown keys raise."""
    if not theta:
        return p
    names = {f.name for f in fields(Params)}
    kw = {}
    for k, v in theta.items():
        if k not in names:
            raise KeyError(f"unknown parameter {k!r}")
        kw[k] = np.asarray(v, dtype=float) if not isinstance(v, (bool, str)) and v is not None else v
    return replace(p, **kw)


def as_dict(p: Params) -> dict:
    return {f.name: getattr(p, f.name) for f in fields(p)}


# ------------------------------------------------------------------------------------------- devices (§4.7)
L_MAX_NITS = 1000.0  # typical maximum luminance used for "150 nit" tests (APL14P; assumption for others)


def nits_to_b(L_nits, p: Params = RP63, L_max=L_MAX_NITS):
    """Slider fraction b giving luminance ``L_nits`` on a panel of maximum ``L_max`` (b = (L/L_max)^(1/γ_b))."""
    return (np.asarray(L_nits, dtype=float) / L_max) ** (1.0 / p.gamma_b)


def device(Q_Ah: float, diag_in: float | None, base: Params = RP63, **overrides) -> Params:
    """Generic validation device (§4.7): R·Q = const scaling (A12) and display area 2.455·d² cm²."""
    kw = dict(Q_n=Q_Ah, R0_ref=base.R0_ref * 4.70 / Q_Ah, R1_ref=base.R1_ref * 4.70 / Q_Ah)
    if diag_in is not None:
        kw["A_disp"] = 2.455 * diag_in ** 2
    kw.update(overrides)
    return replace(base, **kw)


def iphone14pro(**ov) -> Params:
    """iPhone 14 Pro (Bonn, Apple, LMAG22): 3.32 Ah (SCH26 measured), 6.1"."""
    return device(3.32, 6.1, **ov)


def iphone13pro(**ov) -> Params:
    """iPhone 13 Pro (Bonn): 3.16 Ah (SCH26 measured, pre-used), 6.1"."""
    return device(3.16, 6.1, **ov)


def iphone14promax(**ov) -> Params:
    """iPhone 14 Pro Max (NBC22): 4.323 Ah (N9M22), 6.7"."""
    return device(4.323, 6.7, **ov)


def pixel8(**ov) -> Params:
    """Pixel 8 (NBC23): 4.575 Ah typical (GPX), 6.2"."""
    return device(4.575, 6.2, **ov)


def galaxy_s24(**ov) -> Params:
    """Galaxy S24 Exynos (NBC24): 4.000 Ah typical, 6.2"."""
    return device(4.000, 6.2, **ov)


def pixel8a(**ov) -> Params:
    """Pixel 8a (median-user claim, GPX): 4.492 Ah typical, 6.1"."""
    return device(4.492, 6.1, **ov)


_NO_LOADS = dict(P_susp=0.0, P_awake=0.0, E_wake=0.0, P_d0=0.0, p_c=0.0, k_RR=0.0, P_SoC_max=0.0, m_b=0.0,
                 burst=False, throttling=False, load_heat=False)  # cell on a tester: the load is external


def hg2_cell(ocv_table=None, R_sigma=0.038, split=(6, 1), C_th=None, R_th=None, Q_Ah=3.0, **ov) -> Params:
    """LG HG2 18650 cell alone (§4.7, V2): no PCM, no device loads, NMC OCV.

    The [CAL] quantities of V2 are supplied by the caller (see ``validation.hg2_calibration``): the 25 °C
    C/20 OCV table ``ocv_table = (z, U)``, and the cell ``C_th``/``R_th`` fitted to the 25 °C 2C temperature
    rise.  ``R0 + R1 = R_sigma`` (0.038 Ω) split 6:1.  R0 is SOC independent here (A_R = 0) because the
    apparent resistance of the data is an SOC-window average.  All activation energies stay at RP-63 values.
    Defaults for C_th (≈ m·c_p = 46 g × 1.0 J/(g K)) and R_th are placeholders until calibrated.
    """
    s0, s1 = split
    kw = dict(_NO_LOADS, Q_n=Q_Ah, R0_ref=R_sigma * s0 / (s0 + s1), R1_ref=R_sigma * s1 / (s0 + s1), A_R=0.0,
              V_cut=2.8, d_sys=0.0, dV_hv=0.0, ocv_table=ocv_table,
              C_th=46.0 if C_th is None else C_th, R_th=10.0 if R_th is None else R_th, R_skin=10.0)
    kw.update(ov)
    return replace(RP63, **kw)


NASA_B0005 = dict(Q_fresh=1.8622, Q_aged=1.3280, Re_fresh=0.04467, Rct_fresh=0.06946, Re_aged=0.05004,
                  Rct_aged=0.07479)  # nasa_pcoe_summary.json (q_int_Ah, EIS impedance first/last)


def nasa_b0005(aged_cell: bool = False, **ov) -> Params:
    """NASA PCoE B0005 18650 cell (V3).  Measured fade and EIS R_e/R_ct inserted directly.

    R0 <- R_e, R1 <- R_ct (fresh or last EIS), Q_max from the integrated discharge capacity.  A 4.2 V cell, so
    the 4.45 V stretch is removed (dV_hv = 0, Ai2020 4.2 V curve); cut-off 2.7 V (NASA protocol).
    No device loads.  The aged cell is expressed through q_loss with r_growth = 0 and its own R values.
    """
    d = NASA_B0005
    kw = dict(_NO_LOADS, Q_n=d["Q_fresh"], A_R=0.0, dV_hv=0.0, V_cut=2.7, d_sys=0.0, C_th=40.0, R_th=10.0)
    if aged_cell:
        kw.update(q_loss=1.0 - d["Q_aged"] / d["Q_fresh"], r_growth=0.0, R0_ref=d["Re_aged"], R1_ref=d["Rct_aged"])
    else:
        kw.update(R0_ref=d["Re_fresh"], R1_ref=d["Rct_fresh"])
    kw.update(ov)
    return replace(RP63, **kw)


DEVICES = dict(RP63=default, iphone14pro=iphone14pro, iphone13pro=iphone13pro, iphone14promax=iphone14promax,
               pixel8=pixel8, galaxy_s24=galaxy_s24, pixel8a=pixel8a, hg2_cell=hg2_cell, nasa_b0005=nasa_b0005)

# ------------------------------------------------------------------------------------------- UQ spec (§4)
# name: (dist, lo, mode(=default), hi).  'U' uniform, 'T' triangular, 'LU' log-uniform.
UQ_SPEC: dict[str, tuple] = {
    # 4.1 cell and pack
    "Q_n": ("T", 4.56, 4.70, 4.80),
    "dU": ("U", -0.020, 0.0, 0.020),
    "k_S": ("U", 0.5, 1.0, 1.5),
    "R0_ref": ("LU", 0.035, 0.060, 0.130),
    "A_R": ("U", 0.5, 0.80, 2.0),
    "Ea_R0": ("T", 12e3, 18e3, 27e3),
    "R1_ref": ("LU", 0.003, 0.010, 0.025),
    "tau1_ref": ("U", 10.0, 30.0, 60.0),
    "Ea_R1": ("T", 25e3, 30e3, 55e3),
    "Ea_tau": ("U", 15e3, 21e3, 27e3),
    "V_cut": ("U", 3.30, 3.40, 3.50),
    "d_sys": ("U", 0.0, 0.10, 0.20),
    "eta_c": ("U", 0.98, 0.99, 1.00),
    # 4.2 thermal
    "C_th": ("T", 100.0, 150.0, 230.0),
    "R_th": ("T", 2.0, 3.0, 4.0),
    "R_skin": ("T", 4.0, 5.5, 9.0),
    "T_lim_C": ("U", 41.0, 45.0, 48.0),
    "w_thr": ("U", 0.5, 1.0, 3.0),
    "P_fl": ("U", 0.3, 0.5, 1.0),
    "k_env": ("U", 0.8, 1.0, 1.5),
    "dT_ch": ("U", 5.0, 8.0, 12.0),
    "eta_charger": ("U", 0.85, 0.90, 0.95),
    # 4.3 device power
    "P_susp": ("U", 0.010, 0.020, 0.030),
    "P_awake": ("LU", 0.05, 0.090, 0.25),
    "E_wake": ("LU", 0.05, 0.10, 0.30),
    "P_d0": ("U", 0.08, 0.12, 0.25),
    "p_c": ("U", 1.23e-3, 0.25 / 97.4, 5.13e-3),
    "gamma_b": ("U", 1.0, 2.0, 2.9),
    "k_RR": ("U", 1.0e-3, 2.8e-3, 3.3e-3),
    "P_SoC_max": ("U", 3.0, 5.0, 7.0),
    "v0": ("U", 0.35, 0.446, 0.55),
    "P_G": ("T", 0.14, 0.27, 0.45),
    "P_cam": ("U", 0.9, 1.15, 2.0),
    "P_camf": ("U", 0.2, 0.30, 0.5),
    "P_aud": ("U", 0.02, 0.05, 0.15),
    "P_BT": ("U", 0.03, 0.04, 0.05),
    "m_b": ("U", 0.5, 1.0, 2.0),
    # 4.4 radio
    "wifi_P_idle": ("U", 0.002, 0.004, 0.010),
    "wifi_P_act": ("U", 0.10, 0.13, 0.30),
    "wifi_e_b": ("LU", 0.005, 0.015, 0.06),
    "wifi_P_rt": ("U", 0.13, 0.30, 0.8),
    "lte_P_idle": ("U", 0.005, 0.015, 0.035),
    "lte_P_tail": ("U", 0.066, 0.12, 0.178),
    "lte_T_tail": ("U", 5.0, 10.0, 11.6),
    "lte_P_act": ("LU", 0.30, 0.60, 1.30),
    "lte_T_act": ("U", 0.08, 0.10, 0.20),
    "lte_E_pro": ("U", 0.1, 0.31, 0.5),
    "lte_e_b": ("LU", 0.014, 0.030, 0.08),
    "lte_P_rt": ("U", 0.6, 1.20, 2.6),
    "nr_P_idle": ("U", 0.015, 0.030, 0.060),
    "nr_P_tail": ("U", 0.249, 0.255, 0.593),
    "nr_T_tail": ("U", 10.2, 10.3, 21.4),
    "nr_P_act": ("LU", 0.50, 0.90, 1.80),
    "nr_T_act": ("U", 0.08, 0.10, 0.20),
    "nr_E_pro": ("U", 0.5, 1.0, 2.0),
    "nr_e_b": ("LU", 0.014, 0.030, 0.08),
    "nr_P_rt": ("U", 1.2, 1.40, 3.0),
    "k_ul": ("U", 2.0, 5.0, 8.4),
    "gamma_sig": ("LU", 0.10, 0.26, 1.0),
    # 4.5 aging
    "k_c": ("U", 1.05e-3, 0.027 / math.sqrt(365.0), 2.09e-3),
    "E_cal": ("U", 20e3, 58e3, 60e3),
    "g_hi": ("U", 1.3, 1.75, 1.9),
    "k_HV": ("U", 0.0, 1.0, 1.5),
    "w_L": ("U", 0.1, 0.3, 0.5),
    "E_H": ("U", 25e3, 33e3, 45e3),
    "E_L": ("U", 40e3, 55.5e3, 70e3),
    "beta_L": ("U", 0.5, 2.64, 2.64),
    "gamma_R": ("LU", 1.5, 2.5, 7.0),
}

# §4.6 usage-model factors (CTMC).  Multipliers apply to the profile values; 'DIR' entries are Dirichlet
# draws: ('DIR', concentration, base probabilities, None) -> array of shape (k, n).
ACTIVITIES = ("social", "video", "game", "nav", "vcall", "camera", "audio")
USAGE_UQ_SPEC: dict[str, tuple] = {
    "lam_on_mul": ("LU", 0.7, 1.0, 1.3),
    **{f"d_mul_{a}": ("LU", 0.7, 1.0, 1.3) for a in ACTIVITIES},
    "conn_mul_wifi": ("LU", 0.5, 1.0, 2.0),
    "conn_mul_cell": ("LU", 0.5, 1.0, 2.0),
    "p_typical": ("DIR", 50.0, (0.55, 0.12, 0.05, 0.04, 0.04, 0.08, 0.12), None),
    "p_heavy": ("DIR", 50.0, (0.45, 0.20, 0.15, 0.05, 0.05, 0.05, 0.05), None),
}

# §10.3 Sobol factor groups (each group moves jointly by its own §4 distributions).
SOBOL_GROUPS: dict[str, tuple[str, ...]] = {
    "Q_n": ("Q_n",), "R0_ref": ("R0_ref",), "A_R": ("A_R",), "Ea_R0": ("Ea_R0",),
    "RC": ("R1_ref", "tau1_ref", "Ea_R1", "Ea_tau"), "OCV": ("dU", "k_S"), "V_cut": ("V_cut",),
    "d_sys": ("d_sys",), "C_th": ("C_th",), "R_th": ("R_th",), "R_skin": ("R_skin",),
    "throttle": ("T_lim_C", "w_thr", "P_fl"), "k_env": ("k_env",), "floor": ("P_susp", "P_awake", "E_wake"),
    "display": ("P_d0", "p_c", "gamma_b"), "k_RR": ("k_RR",), "SoC": ("P_SoC_max", "v0"),
    "wifi": tuple(k for k in UQ_SPEC if k.startswith("wifi_")),
    "lte": tuple(k for k in UQ_SPEC if k.startswith("lte_")) + ("k_ul",),
    "nr": tuple(k for k in UQ_SPEC if k.startswith("nr_")), "gamma_sig": ("gamma_sig",),
    "peripherals": ("P_G", "P_cam", "P_camf", "P_aud", "P_BT"), "m_b": ("m_b",), "gamma_R": ("gamma_R",),
}

# Units / symbols for the AN0 parameter table.
UNITS = dict(Q_n="Ah", dU="V", R0_ref="Ohm", R1_ref="Ohm", tau1_ref="s", Ea_R0="J/mol", Ea_R1="J/mol",
             Ea_tau="J/mol", V_cut="V", d_sys="V", T_shut_C="degC", C_th="J/K", R_th="K/W", R_skin="K/W",
             T_lim_C="degC", dT_thr="K", w_thr="K", P_fl="W", dT_ch="K", P_susp="W", P_awake="W", E_wake="J",
             P_d0="W", p_c="W/cm2", A_disp="cm2", k_RR="W/Hz", P_SoC_max="W", P_G="W", P_cam="W", P_camf="W",
             P_aud="W", P_BT="W", k_c="d^-1/2", E_cal="J/mol", k_H="FCE^-1/2", k_L="FCE^-1/2", E_H="J/mol",
             E_L="J/mol", beta_L="h", c_ref="1/h", N_80="FCE")
for _r in RATS:
    UNITS.update({f"{_r}_P_idle": "W", f"{_r}_P_tail": "W", f"{_r}_T_tail": "s", f"{_r}_P_act": "W",
                  f"{_r}_T_act": "s", f"{_r}_E_pro": "J", f"{_r}_e_b": "W/(Mbit/s)", f"{_r}_P_rt": "W"})


# ------------------------------------------------------------------------------------------- sampling
def _ppf(dist: str, u: np.ndarray, lo, mode, hi) -> np.ndarray:
    """Inverse CDF of one scalar distribution."""
    if dist == "U":
        return lo + u * (hi - lo)
    if dist == "LU":
        return np.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))
    if dist == "T":
        a, m, b = float(lo), float(mode), float(hi)
        c = (m - a) / (b - a)
        return np.where(u < c, a + np.sqrt(u * (b - a) * (m - a)), b - np.sqrt((1 - u) * (b - a) * (b - m)))
    if dist in ("fixed", "F"):
        return np.full_like(u, float(mode))
    raise ValueError(f"unknown distribution {dist!r}")


def _dims(entry: tuple) -> int:
    return len(entry[2]) if entry[0] == "DIR" else 1


def _resolve_spec(spec) -> dict:
    if spec is None:
        return dict(UQ_SPEC)
    if isinstance(spec, Mapping):
        return dict(spec)
    out = {}
    for k in spec:  # iterable of names
        out[k] = UQ_SPEC[k] if k in UQ_SPEC else USAGE_UQ_SPEC[k]
    return out


def unit_samples(d: int, n: int, method: str = "sobol", seed: int = 2026) -> np.ndarray:
    """(n, d) points in [0, 1)^d: scrambled Sobol' (default), Latin hypercube or plain random."""
    from scipy.stats import qmc

    if method == "sobol":
        with warnings.catch_warnings():  # n need not be a power of two for our use
            warnings.simplefilter("ignore", UserWarning)
            return qmc.Sobol(d, scramble=True, seed=seed).random(n)
    if method == "lhs":
        return qmc.LatinHypercube(d, seed=seed).random(n)
    if method == "random":
        return np.random.default_rng(seed).random((n, d))
    raise ValueError(method)


def map_unit(spec, U: np.ndarray) -> dict[str, np.ndarray]:
    """Map unit-cube points U (n, d) through the inverse CDFs of ``spec`` (in insertion order)."""
    from scipy.stats import gamma

    spec = _resolve_spec(spec)
    out: dict[str, np.ndarray] = {}
    j = 0
    eps = 1e-12
    for name, e in spec.items():
        if e[0] == "DIR":
            alpha = float(e[1]) * np.asarray(e[2], dtype=float)
            k = alpha.size
            u = np.clip(U[:, j:j + k], eps, 1 - eps)
            g = gamma.ppf(u, alpha[None, :])
            out[name] = (g / g.sum(axis=1, keepdims=True)).T  # (k, n)
            j += k
        else:
            dist, lo, mode, hi = e
            out[name] = _ppf(dist, U[:, j], lo, mode, hi)
            j += 1
    return out


def sample(spec=None, n: int = 1024, method: str = "sobol", seed: int = 2026) -> dict[str, np.ndarray]:
    """Joint sample of the factors in ``spec`` (§4 distributions; §9: scrambled Sobol', seed 2026).

    ``spec`` is a dict ``name -> (dist, lo, mode, hi)``, an iterable of names from ``UQ_SPEC``/``USAGE_UQ_SPEC``,
    or None for the whole ``UQ_SPEC``.  Returns ``{name: array (n,)}`` (Dirichlet entries: ``(k, n)``).
    Factors are independent (V_sys = V_cut - d_sys is a derived quantity, so its dependence is automatic).
    """
    spec = _resolve_spec(spec)
    d = sum(_dims(e) for e in spec.values())
    U = unit_samples(d, n, method=method, seed=seed)
    return map_unit(spec, U)


def param_table(p: Params = RP63) -> dict:
    """Resolved parameter table for AN0: {name: {value, unit, dist, lo, mode, hi}}."""
    out = {}
    for f in fields(p):
        v = getattr(p, f.name)
        if f.name in ("ocv_table", "dudt_table", "burst", "throttling", "load_heat", "r_growth"):
            continue
        e = UQ_SPEC.get(f.name, ("fixed", v, v, v))
        out[f.name] = dict(value=v, unit=UNITS.get(f.name, "-"), dist=e[0], lo=e[1], mode=e[2], hi=e[3])
    return out


__all__ = [
    "Params", "RP63", "default", "aged", "with_theta", "as_dict", "device", "iphone14pro", "iphone13pro",
    "iphone14promax", "pixel8", "galaxy_s24", "pixel8a", "hg2_cell", "nasa_b0005", "nits_to_b", "UQ_SPEC",
    "USAGE_UQ_SPEC", "SOBOL_GROUPS", "sample", "unit_samples", "map_unit", "param_table", "RG", "T_REF", "K0",
    "RATS", "RADIO_FIELDS", "ACTIVITIES", "DEVICES",
]
