"""Reference prototype of docs/design/FINAL_SPEC.md (NOT production code).

Computes every number tagged [PROTO] in FINAL_SPEC.md.  Production code lives in src/phonebatt and must
reproduce these values within the tolerances of FINAL_SPEC §15.3.

Run from the repository root:   python3 -I docs/design/scripts/final_proto.py [section ...]
Sections: base scen bonn bench cold aged ctmc mc recs recs2 aging sens   (default: all; full run ≈ 75 s)
"""
import math
import sys
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

RG, TREF, K0 = 8.314, 298.15, 273.15

# --------------------------------------------------------------------------------------------- tables
DUDT_Z = np.round(np.arange(0, 1.0001, 0.05), 2)
DUDT_V = 1e-3 * np.array([-0.4485, -0.7057, -0.5343, -0.4982, -0.4787, -0.4213, -0.3722, -0.3439, -0.3068,
                          -0.2655, -0.2482, -0.2412, -0.2482, -0.2281, -0.2047, -0.1798, -0.1528, -0.1246,
                          -0.0966, -0.0705, -0.0459])  # V/K, docs/research/data/cell_ecm_pybamm/pybamm_ocv_fit.json

RADIO = {  # battery-side, per RAT
    'wifi': dict(idle=0.004, tail=0.12, Ttail=0.24, act=0.13, Tact=0.24, Epro=0.0098, eb=0.015, rt=0.30),
    'lte': dict(idle=0.015, tail=0.12, Ttail=10.0, act=0.60, Tact=0.10, Epro=0.31, eb=0.030, rt=1.20),
    'nr': dict(idle=0.030, tail=0.255, Ttail=10.3, act=0.90, Tact=0.10, Epro=1.0, eb=0.030, rt=1.40),
}


def default_params():
    return dict(
        Qn=4.70, b0=0.60117, b1=36.426, a0=3.6193, a1=0.44311, a2=-0.48605, a3=0.63127,
        dVhv=0.242, phv=3.86, dU=0.0, kS=1.0,
        R0ref=0.060, AR=0.80, z0R=0.139, Ea0=18e3, R1ref=0.010, tau1ref=30.0, Ea1=30e3, Eatau=21e3,
        Vcut=3.40, dsys=0.10,
        Cth=150.0, Rth=3.0, Rskin=5.5, kenv=1.0, Tlim=45.0 + K0, dthr=2.0, wthr=1.0, Pfl=0.5, Tshut=68.0 + K0,
        Psusp=0.020, Pawake=0.090, Pd0=0.12, pc=0.25 / 97.4, Adisp=97.4, APLref=0.70, gb=2.0, kRR=0.0028,
        PSoCmax=5.0, v0=0.446,
        kul=5.0, gsig=0.26, rsrp_ref=-85.0, mmax=3.5, Ewake=0.10,
        PG=0.27, Pcam=1.15, Pcamf=0.30, Paud=0.05, Pbt=0.04, dPbmul=1.0,
        radio={k: dict(v) for k, v in RADIO.items()},
        qloss=0.0, rgrowth=0.0, burst=True,
    )


P0 = default_params()


def a_R(T, E):
    """Arrhenius factor for resistance-like quantities: >1 when colder than T_ref."""
    return math.exp(E / RG * (1.0 / T - 1.0 / TREF))


def dUdT(z):
    return float(np.interp(min(max(z, 0.0), 1.0), DUDT_Z, DUDT_V))


def ocv(z, T, p):
    z = min(max(z, 0.0), 1.0)
    u = (-p['b0'] * math.exp(-p['b1'] * z) + p['a0'] + p['a1'] * z + p['a2'] * z * z + p['a3'] * z ** 3
         + p['dVhv'] * z ** p['phv'] + p['dU'])
    return u + (T - TREF) * p['kS'] * dUdT(z)


def R0f(z, T, p):
    z = max(z, 0.0)
    rinf = p['R0ref'] / (1 + p['AR'] * math.exp(-0.5 / p['z0R']))
    return rinf * (1 + p['AR'] * math.exp(-z / p['z0R'])) * a_R(T, p['Ea0']) * (1 + p['rgrowth'])


def R1f(T, p):
    return p['R1ref'] * a_R(T, p['Ea1']) * (1 + p['rgrowth'])


def tau1f(T, p):
    return p['tau1ref'] * a_R(T, p['Eatau'])


def z_empty(p, T=TREF):
    return brentq(lambda z: ocv(z, T, p) - p['Vcut'], 1e-7, 0.3)


def Eu_Wh(p, z0=1.0):
    ze = z_empty(p)
    return p['Qn'] * (1 - p['qloss']) * quad(lambda z: ocv(z, TREF, p), ze, z0, limit=200)[0]


# --------------------------------------------------------------------------------------------- load model
def msig(rsrp, p):
    return min(p['mmax'], 1.0 + p['gsig'] * max(0.0, (p['rsrp_ref'] - rsrp) / 10.0))


def power_parts(u, p):
    """Return (P_nonSoC, P_SoC_demand, dPb, breakdown) at the battery terminals [W] (E1)."""
    br = {}
    scr = u.get('scr', 0)
    awake = max(scr, u.get('awake', 0))
    br['susp'] = p['Psusp']
    br['awake'] = p['Pawake'] * awake
    b, apl = u.get('b', 0.4), u.get('apl', 0.6) * u.get('aplmul', 1.0)
    br['disp'] = scr * (p['Pd0'] + p['pc'] * p['Adisp'] * (b / 0.5) ** p['gb'] * apl / p['APLref']
                        + p['kRR'] * max(u.get('fr', 60) - 60, 0))
    x = u.get('x', 0.0)
    dem = p['PSoCmax'] * x * (p['v0'] + (1 - p['v0']) * x) ** 2
    R = p['radio']
    net = 0.0
    if u.get('wifi', 0):
        net += R['wifi']['idle']
    cell = u.get('cell', None)
    if cell:
        net += R[cell]['idle']
    via = u.get('via', None)
    lam = u.get('lam', 0.0)
    if via is not None:
        r = R['wifi'] if via == 'wifi' else R[cell]
        m = 1.0 if via == 'wifi' else msig(u.get('rsrp', -85.0), p)
        if u.get('cls', 'burst') == 'rt':
            net += m * r['rt']
        else:
            net += m * ((r['tail'] - r['idle']) * (1 - math.exp(-lam * r['Ttail']))
                        + (r['act'] - r['tail']) * (1 - math.exp(-lam * r['Tact']))
                        + lam * math.exp(-lam * r['Ttail']) * r['Epro']
                        + r['eb'] * (u.get('ddl', 0.0) + p['kul'] * u.get('dul', 0.0)))
    br['net'] = net
    br['wake'] = 0.0 if awake else lam * p['Ewake']
    br['gnss'] = p['PG'] * u.get('gnss', 0.0)
    br['cam'] = p['Pcam'] * u.get('cam', 0.0) + p['Pcamf'] * u.get('camf', 0.0)
    br['aud'] = p['Paud'] * u.get('aud', 0.0) + p['Pbt'] * u.get('bt', 0.0)
    br['sens'] = u.get('sens', 0.0)
    br['extra'] = u.get('Pextra', 0.0)
    br['soc_dem'] = dem
    Pns = sum(v for k, v in br.items() if k != 'soc_dem')
    return Pns, dem, u.get('dPb', 0.0) * p['dPbmul'], br


def throttle(dem, Ts, p):
    if dem <= p['Pfl']:
        return dem
    y = (p['Tlim'] - p['dthr'] - Ts) / p['wthr']
    return p['Pfl'] + (dem - p['Pfl']) / (1.0 + math.exp(-y))


# --------------------------------------------------------------------------------------------- L2 model
class Seg:
    def __init__(self, u, Tamb, p):
        self.u, self.Tamb, self.p = u, Tamb, p
        self.Pns, self.dem, self.dPb, self.br = power_parts(u, p)
        self.Q = p['Qn'] * (1 - p['qloss'])
        self.Vsys = p['Vcut'] - p['dsys']

    def state(self, x):
        p = self.p
        z, v1, T = x
        U = ocv(z, T, p)
        R0 = R0f(z, T, p)
        Ts = self.Tamb + p['Rskin'] / p['Rth'] * (T - self.Tamb)
        Psoc = throttle(self.dem, Ts, p)
        P = self.Pns + Psoc
        Ue = U - v1
        D = Ue * Ue - 4 * R0 * P
        I = 2 * P / (Ue + math.sqrt(D)) if D > 0 else Ue / (2 * R0)
        V = Ue - I * R0
        Pk = P + self.dPb
        Dk = Ue * Ue - 4 * R0 * Pk
        Vk = 0.5 * (Ue + math.sqrt(Dk)) if Dk > 0 else 0.0
        return dict(U=U, R0=R0, I=I, V=V, Vk=Vk, P=P, D=D, Ts=Ts, Psoc=Psoc)

    def rhs(self, t, x):
        p = self.p
        z, v1, T = x
        s = self.state(x)
        I = s['I']
        R1 = R1f(T, p)
        q = I * s['V'] + I * I * s['R0'] + v1 * v1 / R1 - I * T * p['kS'] * dUdT(z)
        return [-I / (3600.0 * self.Q), (R1 * I - v1) / tau1f(T, p),
                (q - (T - self.Tamb) / (p['kenv'] * p['Rth'])) / p['Cth']]

    def gvals(self, x):
        s = self.state(x)
        g = dict(charge=x[0], cutoff=s['V'] - self.p['Vcut'], infeasible=s['D'], thermal=self.p['Tshut'] - x[2])
        if self.p['burst'] and self.dPb > 0:
            g['brownout'] = s['Vk'] - self.Vsys
        return g


def run(schedule, Tamb_C, s0=1.0, T0_C=None, p=P0, tmax_h=1500.0, record=False):
    """schedule: list of (duration_s, u) repeated cyclically. u may contain 'Tamb' (deg C) to override."""
    ze = z_empty(p)
    x = np.array([ze + s0 * (1 - ze), 0.0, (Tamb_C if T0_C is None else T0_C) + K0])
    t, k = 0.0, 0
    tr = []
    while t < tmax_h * 3600:
        dur, u = schedule[k % len(schedule)]
        seg = Seg(u, u.get('Tamb', Tamb_C) + K0, p)
        g0 = seg.gvals(x)
        bad = [n for n, v in g0.items() if not v > 0]
        if bad:
            return finish(t, bad[0], x, ze, seg, tr)
        evs = []
        for name in g0:
            f = (lambda nm: (lambda tt, xx: seg.gvals(xx)[nm]))(name)
            f.terminal, f.direction = True, -1
            f.__name__ = name
            evs.append(f)
        sol = solve_ivp(seg.rhs, (t, t + dur), x, method='LSODA', events=evs, rtol=1e-7,
                        atol=[1e-9, 1e-7, 1e-5], max_step=60.0)
        if record:
            for j in range(sol.t.size):
                st = seg.state(sol.y[:, j])
                tr.append((sol.t[j], sol.y[0, j], sol.y[2, j], st['V'], st['P'], st['Ts']))
        for i, te in enumerate(sol.t_events):
            if len(te):
                return finish(te[0], evs[i].__name__, sol.y_events[i][0], ze, seg, tr)
        x, t = sol.y[:, -1], sol.t[-1]
        k += 1
    return dict(tte_h=np.inf, event='none', s_off=None)


def finish(t, ev, x, ze, seg, tr):
    st = seg.state(x)
    return dict(tte_h=t / 3600, event=ev, s_off=(x[0] - ze) / (1 - ze), T_end=x[2] - K0, Ts_end=st['Ts'] - K0,
                P_end=st['P'], Psoc_end=st['Psoc'], traj=tr)


# --------------------------------------------------------------------------------------------- primitives
def prim():
    U = {}
    sb = dict(scr=0, awake=0, cls='burst', lam=1 / 120, ddl=0.0, dPb=1.0)
    U['SB-W'] = dict(sb, wifi=1, cell='lte', via='wifi')
    U['SB-L'] = dict(sb, wifi=0, cell='lte', via='cell')
    U['SB-5'] = dict(sb, wifi=0, cell='nr', via='cell', lam=60 / 3600)
    U['RD'] = dict(scr=1, b=0.4, apl=0.65, fr=60, x=0.05, wifi=1, cell='lte', via='wifi', lam=1 / 30, ddl=0.05, dPb=3.0)
    U['SO'] = dict(scr=1, b=0.4, apl=0.55, fr=120, x=0.25, wifi=1, cell='lte', via='wifi', lam=2.0, ddl=2.0, dPb=3.0)
    U['SO-L'] = dict(U['SO'], wifi=0, via='cell')
    U['VW'] = dict(scr=1, b=0.4, apl=0.30, fr=60, x=0.10, wifi=1, cell='lte', via='wifi', lam=0.25, ddl=3.0, aud=1, dPb=2.0)
    U['VL'] = dict(U['VW'], wifi=0, via='cell')
    U['V5'] = dict(U['VW'], wifi=0, cell='nr', via='cell')
    U['NV'] = dict(scr=1, b=0.7, apl=0.55, fr=60, x=0.25, wifi=0, cell='lte', via='cell', lam=0.2, ddl=0.1, gnss=1,
                   sens=0.047, dPb=3.0)
    U['GM'] = dict(scr=1, b=0.5, apl=0.40, fr=120, x=0.80, wifi=1, cell='lte', via='wifi', lam=20.0, ddl=0.2, dul=0.05,
                   aud=1, dPb=4.0)
    U['GM-L'] = dict(U['GM'], wifi=0, via='cell')
    U['VC'] = dict(scr=1, b=0.5, apl=0.45, fr=60, x=0.30, wifi=0, cell='lte', via='cell', cls='rt', camf=1, aud=1, dPb=3.0)
    U['VC-W'] = dict(U['VC'], wifi=1, via='wifi')
    U['CM'] = dict(scr=1, b=0.6, apl=0.45, fr=60, x=0.05, wifi=1, cell='lte', via='wifi', lam=0.0, cam=1, dPb=4.0)
    U['AU'] = dict(scr=0, awake=1, wifi=1, cell='lte', via='wifi', lam=1 / 30, ddl=0.256, aud=1, bt=1, dPb=1.0)
    U['AU-L'] = dict(U['AU'], wifi=0, via='cell')
    U['TRK'] = dict(scr=0, awake=1, wifi=0, cell='lte', via='cell', lam=1 / 60, ddl=0.01, gnss=1, sens=0.047, dPb=1.0)
    U['BU'] = dict(scr=1, b=1.0, apl=0.45, fr=60, x=0.50, wifi=0, cell='lte', via='cell', lam=1.0, dul=4.0, cam=1, dPb=0.0)
    return U


PR = prim()


def Pmean(name, p=P0, Ts_C=25.0, **ov):
    u = dict(PR[name], **ov)
    Pns, dem, dPb, br = power_parts(u, p)
    return Pns + throttle(dem, Ts_C + K0, p), br


def S12():
    U = PR
    s = [(30, 'RD'), (45, 'NV')] + [(26, 'SB-W'), (4, 'SO')] * 7 + [(30, 'VW'), (30, 'SO')] + \
        [(26, 'SB-W'), (4, 'SO')] * 8 + [(45, 'NV'), (30, 'VC'), (60, 'GM'), (120, 'VW'), (60, 'SO'), (60, 'RD'),
                                         (480, 'SB-W')]
    assert sum(m for m, _ in s) == 1440
    return [(m * 60.0, U[n]) for m, n in s]


def S9(Tamb):
    return [(570.0, dict(PR['TRK'], Tamb=Tamb)), (30.0, dict(PR['BU'], Tamb=Tamb))]


# --------------------------------------------------------------------------------------------- L1 quadrature
def tte_L1(P, p=P0, T_C=25.0, s0=1.0):
    T = T_C + K0
    Q = p['Qn'] * (1 - p['qloss'])
    ze0 = z_empty(p)
    z0 = ze0 + s0 * (1 - ze0)

    def Rs(z):
        return R0f(z, T, p) + R1f(T, p)

    def Vqs(z):
        U = ocv(z, T, p)
        return 0.5 * (U + math.sqrt(max(U * U - 4 * Rs(z) * P, 0.0)))

    zc = brentq(lambda z: Vqs(z) - p['Vcut'], 1e-6, z0)
    return 3600 * Q / P * quad(Vqs, zc, z0, limit=200)[0] / 3600, zc


# --------------------------------------------------------------------------------------------- CTMC
ACTS = ['idle', 'social', 'video', 'game', 'nav', 'vcall', 'camera', 'audio']
MAPW = dict(idle='SB-W', social='SO', video='VW', game='GM', nav='NV', vcall='VC-W', camera='CM', audio='AU')
MAPL = dict(idle='SB-L', social='SO-L', video='VL', game='GM-L', nav='NV', vcall='VC', camera='CM', audio='AU-L')
PROFILES = {
    'typical': (4.0, dict(social=.55, video=.12, game=.05, nav=.04, vcall=.04, camera=.08, audio=.12),
                dict(social=5, video=15, game=20, nav=20, vcall=10, camera=2, audio=30)),
    'heavy': (4.0, dict(social=.45, video=.20, game=.15, nav=.05, vcall=.05, camera=.05, audio=.05),
              dict(social=8, video=25, game=30, nav=25, vcall=15, camera=3, audio=30)),
}


def gen_activity(prof):
    lam, pj, dmin = PROFILES[prof]
    n = len(ACTS)
    L = np.zeros((n, n))
    for j, a in enumerate(ACTS[1:], 1):
        L[0, j] = lam * pj[a]
        L[j, 0] = 60.0 / dmin[a]
    np.fill_diagonal(L, -L.sum(1))
    return L


def gen_product(prof, mean_w=2.5, mean_c=2.5):
    La = gen_activity(prof)
    Lc = np.array([[-1 / mean_w, 1 / mean_w], [1 / mean_c, -1 / mean_c]])
    return np.kron(La, np.eye(2)) + np.kron(np.eye(len(ACTS)), Lc)


def stationary(L):
    n = L.shape[0]
    A = np.vstack([L.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1
    return np.linalg.lstsq(A, b, rcond=None)[0]


def moments(L, Pv, Eu):
    pi = stationary(L)
    n = len(pi)
    Pi = np.outer(np.ones(n), pi)
    Dm = np.linalg.inv(Pi - L) - Pi
    Pbar = pi @ Pv
    hvec = Dm @ (Pv - Pbar)
    s2 = 2 * np.sum(pi * (Pv - Pbar) * hvec)
    hP = np.sum(pi * Pv * hvec) / Pbar
    return dict(pi=pi, Pbar=Pbar, s2=s2, h=hvec, mu_fluid=Eu / Pbar, mu_dynkin=(Eu - hvec[0] + hP) / Pbar,
                sd=math.sqrt(s2 * Eu / Pbar ** 3))


def state_powers(p=P0, capped=True):
    """Quasi-steady state powers for the CTMC (game at its thermal fixed point at 25 C if capped)."""
    Pv = []
    for a in ACTS:
        for M in (MAPW, MAPL):
            Pv.append(steady_power(M[a], p) if capped else Pmean(M[a], p)[0])
    return np.array(Pv)


def steady_power(name, p=P0, Tamb_C=25.0):
    """Fixed point T_s* = T_amb + k_env R_skin (P(T_s*) + q_batt) ignoring q_batt (F5)."""
    u = PR[name]
    Pns, dem, _, _ = power_parts(u, p)
    f = lambda Ts: Ts - (Tamb_C + K0) - p['kenv'] * p['Rskin'] * (Pns + throttle(dem, Ts, p))
    Ts = brentq(f, Tamb_C + K0 - 1, Tamb_C + K0 + 200)
    return Pns + throttle(dem, Ts, p)


def mc_L0(L, Pv, Eu, N=20000, seed=1, start=0):
    rng = np.random.default_rng(seed)
    rates = -np.diag(L)
    J = (L - np.diag(np.diag(L))) / rates[:, None]
    cum = np.cumsum(J, 1)
    out = np.empty(N)
    for k in range(N):
        s, E, t = start, Eu, 0.0
        while True:
            tau = rng.exponential(1 / rates[s])
            e = Pv[s] * tau
            if e >= E:
                t += E / Pv[s]
                break
            E -= e
            t += tau
            s = int(np.searchsorted(cum[s], rng.random() * cum[s, -1], side='right'))
        out[k] = t
    return out


def mc_L2(L, prof_states, N=400, h=10.0, Tamb_C=25.0, seed=7, tmax_h=120.0, p=P0):
    """Vectorised RK4 of the L2 model under CTMC usage (state held constant over each step)."""
    rng = np.random.default_rng(seed)
    nmax = int(tmax_h * 3600 / h)
    rates = -np.diag(L)
    J = (L - np.diag(np.diag(L))) / rates[:, None]
    cum = np.cumsum(J, 1)
    S = np.zeros((nmax, N), dtype=np.int16)
    for k in range(N):
        s, t, i0 = 0, 0.0, 0
        while i0 < nmax:
            t += rng.exponential(1 / rates[s]) * 3600
            i1 = min(int(t / h), nmax)
            S[i0:i1, k] = s
            i0 = i1
            s = int(np.searchsorted(cum[s], rng.random() * cum[s, -1], side='right'))
    parts = [power_parts(PR[n], p) for n in prof_states]
    Pns = np.array([q[0] for q in parts])
    Dem = np.array([q[1] for q in parts])
    dPb = np.array([q[2] for q in parts])
    TaK = Tamb_C + K0
    ze = z_empty(p)
    z = np.ones(N)
    v1 = np.zeros(N)
    T = np.full(N, TaK)
    dz_tab = np.interp
    Q = p['Qn'] * (1 - p['qloss'])
    rinf = p['R0ref'] / (1 + p['AR'] * math.exp(-0.5 / p['z0R']))

    def ocv_v(z, T):
        zc = np.clip(z, 0, 1)
        u = (-p['b0'] * np.exp(-p['b1'] * zc) + p['a0'] + p['a1'] * zc + p['a2'] * zc ** 2 + p['a3'] * zc ** 3
             + p['dVhv'] * zc ** p['phv'])
        return u + (T - TREF) * np.interp(zc, DUDT_Z, DUDT_V)

    def arr(E, T):
        return np.exp(E / RG * (1 / T - 1 / TREF))

    def st(z, v1, T, s):
        Ts = TaK + p['Rskin'] / p['Rth'] * (T - TaK)
        d = Dem[s]
        sig = 1 / (1 + np.exp(-(p['Tlim'] - p['dthr'] - Ts) / p['wthr']))
        Psoc = np.where(d > p['Pfl'], p['Pfl'] + (d - p['Pfl']) * sig, d)
        P = Pns[s] + Psoc
        R0 = rinf * (1 + p['AR'] * np.exp(-np.maximum(z, 0) / p['z0R'])) * arr(p['Ea0'], T) * (1 + p['rgrowth'])
        Ue = ocv_v(z, T) - v1
        D = Ue * Ue - 4 * R0 * P
        sq = np.sqrt(np.maximum(D, 0))
        I = 2 * P / (Ue + sq)
        V = Ue - I * R0
        Dk = Ue * Ue - 4 * R0 * (P + dPb[s])
        Vk = np.where(Dk > 0, 0.5 * (Ue + np.sqrt(np.maximum(Dk, 0))), 0.0)
        return I, V, Vk, R0, D

    def f(z, v1, T, s):
        I, V, Vk, R0, D = st(z, v1, T, s)
        R1 = p['R1ref'] * arr(p['Ea1'], T) * (1 + p['rgrowth'])
        tau = p['tau1ref'] * arr(p['Eatau'], T)
        q = I * V + I * I * R0 + v1 * v1 / R1 - I * T * np.interp(np.clip(z, 0, 1), DUDT_Z, DUDT_V)
        return -I / (3600 * Q), (R1 * I - v1) / tau, (q - (T - TaK) / (p['kenv'] * p['Rth'])) / p['Cth']

    alive = np.ones(N, bool)
    tte = np.full(N, np.nan)
    soff = np.full(N, np.nan)
    ev = np.zeros(N, int)
    Vsys = p['Vcut'] - p['dsys']
    for i in range(nmax):
        s = S[i]
        I, V, Vk, R0, D = st(z, v1, T, s)
        g1 = V - p['Vcut']
        g2 = Vk - Vsys
        dead = alive & ((g1 <= 0) | (g2 <= 0) | (D <= 0))
        if dead.any():
            tte[dead] = i * h / 3600
            soff[dead] = (z[dead] - ze) / (1 - ze)
            ev[dead] = np.where(g1[dead] <= 0, 1, 2)
            alive &= ~dead
        if not alive.any():
            break
        k1 = f(z, v1, T, s)
        k2 = f(z + .5 * h * k1[0], v1 + .5 * h * k1[1], T + .5 * h * k1[2], s)
        k3 = f(z + .5 * h * k2[0], v1 + .5 * h * k2[1], T + .5 * h * k2[2], s)
        k4 = f(z + h * k3[0], v1 + h * k3[1], T + h * k3[2], s)
        z = np.where(alive, z + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), z)
        v1 = np.where(alive, v1 + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]), v1)
        T = np.where(alive, T + h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]), T)
    return tte, soff, ev


# --------------------------------------------------------------------------------------------- validation configs
def device(Q, diag_in, **ov):
    p = default_params()
    p['Qn'] = Q
    p['R0ref'] *= 4.70 / Q
    p['R1ref'] *= 4.70 / Q
    p['Adisp'] = 2.455 * diag_in ** 2
    p.update(ov)
    return p


def bonn_cases():
    airplane = dict(scr=1, b=0.4, apl=0.0, fr=60, x=0.0, wifi=0, cell=None, via=None, dPb=0.0)
    return {
        'idle': airplane,
        'cam': dict(airplane, cam=1),
        'gps': dict(airplane, gnss=1),
        'cxl': dict(airplane, cam=1, gnss=1),
    }


def main(sections):
    p = P0
    ze = z_empty(p)
    Eu = Eu_Wh(p)
    if 'base' in sections:
        print('== base ==')
        print(' z_e=%.4f  U(1)=%.4f  mean OCV=%.4f  E_n=%.3f Wh  E_u=%.3f Wh' % (
            ze, ocv(1, TREF, p), quad(lambda z: ocv(z, TREF, p), 0, 1)[0], p['Qn'] * quad(lambda z: ocv(z, TREF, p), 0, 1)[0], Eu))
        for z in (0.1, 0.2, 0.3, 0.5, 0.8, 1.0):
            print('  U(%.1f)=%.3f' % (z, ocv(z, TREF, p)), end='')
        print()
        for Tc in (40, 25, 10, 0, -10, -20):
            T = Tc + K0
            print('  T=%3d C: R0(0.5)=%.1f mOhm R1=%.1f tau1=%.1f s' % (Tc, 1e3 * R0f(0.5, T, p), 1e3 * R1f(T, p), tau1f(T, p)))
        print(' primitive powers (W, 25 C, before throttling):')
        for n in PR:
            P, br = Pmean(n, p)
            print('  %-5s %.3f  %s' % (n, P, {k: round(v, 3) for k, v in br.items() if v}))
        print(' steady (throttled) GM 25C: %.3f  35C: %.3f  15C: %.3f' % (
            steady_power('GM', p, 25), steady_power('GM', p, 35), steady_power('GM', p, 15)))
    if 'scen' in sections:
        print('== scenarios (L2, fresh, 25 C, s0=100%) ==')
        rows = [('S1 SB-W', [(3600.0, PR['SB-W'])], 25), ('S2a SB-L', [(3600.0, PR['SB-L'])], 25),
                ('S2b SB-5', [(3600.0, PR['SB-5'])], 25), ('S3 RD', [(3600.0, PR['RD'])], 25),
                ('S4 SO', [(3600.0, PR['SO'])], 25), ('S5a VW', [(3600.0, PR['VW'])], 25),
                ('S5b VL', [(3600.0, PR['VL'])], 25), ('S5c V5', [(3600.0, PR['V5'])], 25),
                ('S6 NV', [(3600.0, PR['NV'])], 25), ('S7 GM 25', [(3600.0, PR['GM'])], 25),
                ('S7 GM 35', [(3600.0, PR['GM'])], 35), ('S7 GM 15', [(3600.0, PR['GM'])], 15),
                ('S8a VC', [(3600.0, PR['VC'])], 25), ('S8b VC-W', [(3600.0, PR['VC-W'])], 25),
                ('S16 AU', [(3600.0, PR['AU'])], 25), ('CM', [(3600.0, PR['CM'])], 25),
                ('S12 day', S12(), 25)]
        for name, sch, Ta in rows:
            r = run(sch, Ta, p=p)
            print('  %-10s TTE=%7.2f h  event=%-9s s_off=%5.2f%%  T_end=%.1f Ts_end=%.1f P_end=%.3f' % (
                name, r['tte_h'], r['event'], 100 * r['s_off'], r['T_end'], r['Ts_end'], r['P_end']))
        pnt = dict(p, burst=False, Tlim=1e9)
        r = run([(3600.0, PR['GM'])], 25, p=pnt)
        print('  S7 GM 25 no throttle: %.2f h' % r['tte_h'])
        for s0 in (0.5, 0.2):
            r = run(S12(), 25, s0=s0, p=p)
            print('  S12 from %d%%: %.2f h (%s)' % (100 * s0, r['tte_h'], r['event']))
        E_day = sum(d * Pmean_seg(u, p) for d, u in S12()) / 3600
        print('  S12 mean demand P = %.4f W (energy/day %.2f Wh)' % (E_day / 24, E_day))
    if 'bonn' in sections:
        print('== Bonn (V1) ==')
        obs = {'14': dict(idle=55.03, cam=8.43, gps=48.49, cxl=8.38), '13': dict(idle=50.82, cam=9.53, gps=29.76, cxl=8.33)}
        cases = bonn_cases()
        for awake in (0.10, None):
            for dev, Q in (('14', 3.32), ('13', 3.16)):
                pd = device(Q, 6.1)
                if awake is not None:
                    pd['Pawake'] = awake
                for c, u in cases.items():
                    r = run([(3600.0, u)], 25, p=pd)
                    Pm = Pmean_seg(u, pd)
                    print('  aw=%s %s Pro %-4s P=%.3f TTE=%6.2f obs=%6.2f err=%+.1f%% (%s)' % (
                        awake, dev, c, Pm, r['tte_h'], obs[dev][c], 100 * (r['tte_h'] / obs[dev][c] - 1), r['event']))
        pd = device(3.32, 6.1)
        Eu14 = Eu_Wh(pd)
        f = lambda aw: run([(3600.0, cases['idle'])], 25, p=dict(pd, Pawake=aw))['tte_h'] - 55.03
        aw = brentq(f, 0.03, 0.2, xtol=1e-5)
        print('  E_u(14 Pro)=%.3f Wh; calibrated Pawake for 55.03 h: %.4f W (floor %.4f W)' % (
            Eu14, aw, aw + pd['Psusp'] + pd['Pd0']))
        # additivity of the model
        tt = {c: run([(3600.0, u)], 25, p=pd)['tte_h'] for c, u in cases.items()}
        res = (1 / tt['cxl']) / (1 / tt['cam'] + 1 / tt['gps'] - 1 / tt['idle']) - 1
        print('  model additivity residual 14 Pro: %+.3f%%' % (100 * res))
        # slope test D1: mean OCV in display windows at constant P
        for lo, hi in ((0.99, 0.50), (0.50, 0.20), (0.20, 0.01)):
            zl, zh = ze + lo * (1 - ze), ze + hi * (1 - ze)
            m = quad(lambda z: ocv(z, TREF, pd), zh, zl)[0] / (zl - zh)
            print('  window %.2f-%.2f mean OCV %.4f V' % (lo, hi, m))
    if 'bench' in sections:
        print('== published band (V6) ==')
        def nb(L):
            return math.sqrt(L / 1000.0)  # b such that (b)^2 = luminance fraction of a 1000-nit panel
        so = PR['SO']
        cases = [
            ('14P local video (APL14P)', device(3.32, 6.1), dict(PR['VW'], lam=0, ddl=0, bt=1, aud=0), 23.0),
            ('14P streamed video (APL14P)', device(3.32, 6.1), dict(PR['VW'], bt=1, aud=0), 20.0),
            ('14P audio BT (APL14P)', device(3.32, 6.1), PR['AU'], 75.0),
            ('14P audio BT, DSP offload (X-var)', device(3.32, 6.1), dict(PR['AU'], awake=0), 75.0),
            ('14PM Wi-Fi web 150nit (NBC22)', device(4.323, 6.7), dict(so, b=nb(150)), 15.8),
            ('14PM local video 150nit (NBC22)', device(4.323, 6.7), dict(PR['VW'], lam=0, ddl=0, b=nb(150), wifi=0, via=None), 28.13),
            ('14PM idle min bright (NBC22)', device(4.323, 6.7), dict(scr=1, b=0.0, apl=0.65, fr=60, wifi=0, cell='lte', via=None), 57.9),
            ('14P cell web 150nit (LMAG22)', device(3.32, 6.1), dict(PR['SO-L'], b=nb(150)), 10.37),
            ('Pixel8 Wi-Fi web (NBC23)', device(4.575, 6.2), dict(so, b=nb(150)), 13.57),
            ('Pixel8 local video (NBC23)', device(4.575, 6.2), dict(PR['VW'], lam=0, ddl=0, b=nb(150), wifi=0, via=None), 22.72),
            ('Pixel8 idle min bright (NBC23)', device(4.575, 6.2), dict(scr=1, b=0.0, apl=0.65, fr=60, wifi=0, cell='lte', via=None), 24.88),
            ('S24 Wi-Fi web (NBC24)', device(4.0, 6.2), dict(so, b=nb(150)), 14.67),
            ('S24 local video (NBC24)', device(4.0, 6.2), dict(PR['VW'], lam=0, ddl=0, b=nb(150), wifi=0, via=None), 20.45),
            ('S24 idle min bright (NBC24)', device(4.0, 6.2), dict(scr=1, b=0.0, apl=0.65, fr=60, wifi=0, cell='lte', via=None), 39.97),
        ]
        for name, pd, u, obs in cases:
            r = run([(3600.0, u)], 25, p=pd)
            Pm = Pmean_seg(u, pd)
            print('  %-34s P=%.3f TTE=%6.2f obs=%6.2f err=%+.1f%%' % (name, Pm, r['tte_h'], obs, 100 * (r['tte_h'] / obs - 1)))
        print('  GM capped at b=1: %.3f W' % steady_power('GM', dict(p)) , ' (b=1 demand %.3f)' % Pmean('GM', p, b=1.0)[0])
        pd = device(4.492, 6.1)
        print('  Pixel 8a S12 day: %.2f h (GPX >= 31 h)' % run(S12(), 25, p=pd)['tte_h'])
    if 'cold' in sections:
        print('== cold hike S9 and stranded map ==')
        print('  BU burst power %.3f W, TRK %.3f W' % (Pmean('BU', p)[0], Pmean('TRK', p)[0]))
        for Ta in (25, -10, -20):
            for q in (0.0, 0.2):
                pp = dict(p, qloss=q, rgrowth=2.5 * q)
                r = run(S9(Ta), Ta, T0_C=20.0, p=pp)
                print('  S9 Tamb=%d SOH=%.1f: TTE=%.2f h  %s  s_off=%.1f%%' % (Ta, 1 - q, r['tte_h'], r['event'], 100 * r['s_off']))
        print('  stranded map (S4 social, s_off %, TTE h):')
        for Ta in (25, 10, 0, -10, -20):
            row = []
            for q in (0.0, 0.1, 0.2):
                pp = dict(p, qloss=q, rgrowth=2.5 * q)
                r = run([(3600.0, PR['SO'])], Ta, p=pp)
                row.append('%5.1f%% (%5.2f %s)' % (100 * r['s_off'], r['tte_h'], r['event'][:3]))
            print('   %4d C: ' % Ta + ' | '.join(row))
        pp = dict(p, qloss=0.2, rgrowth=0.5)
        for mul in (1.0, 0.5, 0.0):
            r = run([(3600.0, PR['SO'])], -10, p=dict(pp, dPbmul=mul))
            print('  governor dPb x%.1f aged -10C social: TTE=%.2f s_off=%.1f%%' % (mul, r['tte_h'], 100 * r['s_off']))
        r = run([(3600.0, PR['NV'])], -10, p=p)
        r2 = run([(3600.0, PR['NV'])], 25, p=p)
        print('  NV -10C fresh: %.2f h vs 25C %.2f h' % (r['tte_h'], r2['tte_h']))
    if 'aged' in sections:
        print('== aged S12 day / S4 ==')
        for Ta in (25, 0):
            for q in (0.0, 0.2):
                pp = dict(p, qloss=q, rgrowth=2.5 * q)
                r = run(S12(), Ta, p=pp)
                print('  S12 Ta=%d SOH=%.1f: %.2f h (%s, s_off %.1f%%)' % (Ta, 1 - q, r['tte_h'], r['event'], 100 * r['s_off']))
        for rs in (-85, -95, -105, -110, -120):
            r = run([(3600.0, dict(PR['SO-L'], rsrp=rs))], 25, p=p)
            r2 = run([(3600.0, dict(PR['VC'], rsrp=rs))], 25, p=p)
            print('  S11 RSRP=%d: SO-L %.2f h (P %.3f)  VC %.2f h (P %.3f)' % (
                rs, r['tte_h'], Pmean('SO-L', p, rsrp=rs)[0], r2['tte_h'], Pmean('VC', p, rsrp=rs)[0]))
    if 'ctmc' in sections or 'mc' in sections:
        print('== CTMC ==')
        for prof in ('typical', 'heavy'):
            L = gen_product(prof)
            Pv = state_powers(p)
            m = moments(L, Pv, Eu)
            pia = m['pi'].reshape(-1, 2).sum(1)
            Pa = Pv.reshape(-1, 2)
            share = {a: float(m['pi'].reshape(-1, 2)[i] @ Pa[i]) / m['Pbar'] for i, a in enumerate(ACTS)}
            scr_on = 1 - pia[0] - pia[-1]
            print(' %s: Pbar=%.4f W  s2=%.4f W2h  fluid=%.2f h  dynkin=%.2f h  sd=%.2f h  CV=%.3f  screen-on=%.3f (%.1f h/day)' % (
                prof, m['Pbar'], m['s2'], m['mu_fluid'], m['mu_dynkin'], m['sd'], m['sd'] / m['mu_dynkin'], scr_on, 24 * scr_on))
            print('   pi:', {a: round(float(x), 4) for a, x in zip(ACTS, pia)})
            print('   energy share:', {a: round(v, 3) for a, v in share.items()})
            # hub closed-form on Wi-Fi-only chain (activity only) and variance attribution
            lam, pj, dmin = PROFILES[prof]
            La = gen_activity(prof)
            PaW = np.array([steady_power(MAPW[a], p) for a in ACTS])
            mm = moments(La, PaW, Eu)
            dh = {a: dmin[a] / 60 for a in ACTS[1:]}
            EL = 1 / lam + sum(pj[a] * dh[a] for a in ACTS[1:])
            contrib = {a: 2 * pj[a] * dh[a] ** 2 * (PaW[i + 1] - mm['Pbar']) ** 2 / EL for i, a in enumerate(ACTS[1:])}
            s2h = sum(contrib.values())
            print('   hub (Wi-Fi only): s2 dev-matrix=%.6f closed=%.6f ; variance shares %s' % (
                mm['s2'], s2h, {a: round(v / s2h, 3) for a, v in contrib.items()}))
            if 'mc' in sections:
                t0 = mc_L0(L, Pv, Eu, N=20000)
                q = np.percentile(t0, [5, 50, 95])
                print('   MC L0 (20000): mean=%.2f sd=%.2f q5/50/95=%s P(TTE<24h)=%.3f' % (
                    t0.mean(), t0.std(), np.round(q, 2), np.mean(t0 < 24)))
                states = [M[a] for a in ACTS for M in (MAPW, MAPL)]
                tt, so, ev = mc_L2(L, states, N=400, h=10.0)
                q2 = np.nanpercentile(tt, [5, 50, 95])
                print('   MC L2 (400, h=10 s): mean=%.2f sd=%.2f q5/50/95=%s  ratio L2/L0 mean=%.3f  s_off mean=%.2f%% brownout frac=%.2f' % (
                    np.nanmean(tt), np.nanstd(tt), np.round(q2, 2), np.nanmean(tt) / t0.mean(), 100 * np.nanmean(so),
                    np.mean(ev == 2)))
    if 'recs' in sections:
        print('== recommendations on S12 (L2) ==')
        base = run(S12(), 25, p=p)['tte_h']
        print('  base S12: %.2f h' % base)

        def mod(fn):
            return [(d, fn(dict(u))) for d, u in S12()]
        levers = {
            'no game (GM->SB-W)': mod(lambda u: dict(PR['SB-W']) if u.get('x') == 0.80 else u),
            'dark theme APL x0.25 (UI states)': mod(lambda u: dict(u, aplmul=0.25) if u.get('scr') and u.get('apl') in (0.65, 0.55) else u),
            'brightness -0.1': mod(lambda u: dict(u, b=u['b'] - 0.1) if u.get('scr') else u),
            '60 Hz cap': mod(lambda u: dict(u, fr=60) if u.get('scr') else u),
            'voice instead of video call': mod(lambda u: dict(u, scr=0, awake=1, camf=0, x=0.05, cls='rt') if u.get('camf') else u),
            'half background rate': mod(lambda u: dict(u, lam=u['lam'] / 2) if (not u.get('scr') and u.get('lam')) else u),
            'standby on 5G NSA (no Wi-Fi)': mod(lambda u: dict(PR['SB-5']) if u is not None and u.get('scr') == 0 and u.get('awake', 0) == 0 else u),
            'weak signal -110 dBm all day': mod(lambda u: dict(u, rsrp=-110)),
            'GNSS off when not navigating': None,
        }
        for k, sch in levers.items():
            if sch is None:
                continue
            r = run(sch, 25, p=p)
            print('  %-34s %.2f h (%+.1f%%)' % (k, r['tte_h'], 100 * (r['tte_h'] / base - 1)))
        # energy shares of S12 by component
        tot = {}
        for d, u in S12():
            Pns, dem, _, br = power_parts(u, p)
            br = dict(br)
            br['soc'] = throttle(dem, 25 + K0 + 5.5 * (Pns + dem), p) if dem > 0 else 0.0
            for kk, v in br.items():
                if kk == 'soc_dem':
                    continue
                tot[kk] = tot.get(kk, 0.0) + v * d
        S = sum(tot.values())
        print('  S12 energy shares:', {k: round(v / S, 3) for k, v in sorted(tot.items(), key=lambda x: -x[1]) if v > 0})
    if 'aging' in sections:
        aging_proto()
    if 'recs2' in sections:
        recs_smooth()
    if 'sens' in sections:
        print('== elasticities (central +-5 %) ==')
        names = ['Qn', 'PSoCmax', 'pc', 'Pd0', 'Pawake', 'R0ref', 'Ea0', 'R1ref', 'dPbmul', 'Rskin', 'Rth', 'Cth']
        for label, sch, Ta, q in (('S4 25C', [(3600.0, PR['SO'])], 25, 0.0), ('S4 aged -10C', [(3600.0, PR['SO'])], -10, 0.2),
                                  ('S7 25C', [(3600.0, PR['GM'])], 25, 0.0)):
            out = []
            for nm in names:
                vals = []
                for f in (0.95, 1.05):
                    pp = dict(p, qloss=q, rgrowth=2.5 * q)
                    pp[nm] = pp[nm] * f
                    vals.append(run(sch, Ta, p=pp)['tte_h'])
                out.append('%s %+.3f' % (nm, (math.log(vals[1]) - math.log(vals[0])) / (math.log(1.05) - math.log(0.95))))
            print('  %-13s ' % label + ', '.join(out))


def Pmean_seg(u, p):
    Pns, dem, _, _ = power_parts(u, p)
    return Pns + dem


# --------------------------------------------------------------------------------------------- aging (L3 slow layer)
def a_deg(T, E):
    return math.exp(-E / RG * (1.0 / T - 1.0 / TREF))


U08 = None


def gz(z, kHV=1.0):
    """Calendar SOC stress: KEIL16 graphite step x SCH14 linear-voltage factor above z = 0.8 (E16)."""
    global U08
    if U08 is None:
        U08 = ocv(0.8, TREF, P0)
    s = lambda y: 1 / (1 + math.exp(-y))
    g = 0.7 + 0.3 * s((z - 0.2) / 0.02) + 0.75 * s((z - 0.6) / 0.03)
    U = ocv(z, TREF, P0)
    hv = max(1.0, (7.543 * U - 23.75) / (7.543 * U08 - 23.75))
    return g * hv ** kHV


AG = dict(kc=0.027 / math.sqrt(365.0), Ecal=58e3, EH=33e3, EL=55.5e3, betaL=2.64, cref=1.0, wL=0.3, gR=2.5,
          dT_ch=8.0, kHV=1.0)


def day(strategy, Tamb_C, dsoc, c_cc=1.0, c_cv=0.3, z_cc=0.8, dT_use=1.32, heat_h=0.0, heat_C=45.0, a=AG,
        cal_mode=False, Tch_amb=None):
    """Hourly-resolution day: (dt_h, z, T_K, dn_tot, dn_ch, c_ch) segments. Day starts 07:00 unplugged."""
    seg = []
    ztop = 0.8 if strategy == 'limit80' else 1.0
    zlow = max(ztop - dsoc, 0.0)
    nd = 16
    for i in range(nd):  # 07:00-23:00 discharge, linear in SOC
        T = (Tamb_C + (0 if cal_mode else dT_use)) + K0
        if heat_h and i < heat_h:
            T = heat_C + K0
        seg.append((1.0, ztop - dsoc * (i + .5) / nd, T, (dsoc / nd) / 2, 0.0, 0.0))
    z = zlow
    t = 0.0
    Tch = (Tamb_C if Tch_amb is None else Tch_amb) + (0 if cal_mode else a['dT_ch']) + K0
    target = 0.8 if strategy in ('limit80', 'obc') else 1.0
    while z < target - 1e-9 and t < 8.0:
        c = c_cc if z < z_cc - 1e-9 else c_cv
        zstop = min(target, z_cc) if z < z_cc - 1e-9 else target
        dt = min(0.25, (zstop - z) / c, 8.0 - t)
        seg.append((dt, z + c * dt / 2, Tch, c * dt / 2, c * dt, c))
        z += c * dt
        t += dt
    hold_end = 8.0 - (1.0 if strategy == 'obc' else 0.0)
    if t < hold_end:
        seg.append((hold_end - t, z, Tamb_C + K0, 0.0, 0.0, 0.0))
        t = hold_end
    if strategy == 'obc':
        while z < 1.0 - 1e-9 and t < 8.0:
            dt = min(0.25, (1.0 - z) / c_cv, 8.0 - t)
            seg.append((dt, z + c_cv * dt / 2, Tch, c_cv * dt / 2, c_cv * dt, c_cv))
            z += c_cv * dt
            t += dt
        if t < 8.0:
            seg.append((8.0 - t, z, Tamb_C + K0, 0.0, 0.0, 0.0))
    return seg


def age(strategy, Tamb_C, days, dsoc_fn, a=AG, kH=None, kL=None, cal_mode=False, **kw):
    QC = QH = QL = 0.0  # squared losses (E13'-E15')
    hist = []
    for d in range(days):
        q = math.sqrt(QC) + math.sqrt(QH) + math.sqrt(QL)
        for dt, z, T, dn, dnch, c in day(strategy, Tamb_C, dsoc_fn(q), a=a, cal_mode=cal_mode, **kw):
            QC += (a['kc'] * a_deg(T, a['Ecal']) * gz(z, a['kHV'])) ** 2 * dt / 24.0
            QH += (kH * a_deg(T, a['EH'])) ** 2 * dn
            if dnch > 0:
                QL += (kL * a_R(T, a['EL']) * math.exp(a['betaL'] * (c - a['cref']))) ** 2 * dnch
        hist.append((d + 1, math.sqrt(QC), math.sqrt(QH), math.sqrt(QL)))
    return hist


def aging_proto():
    print('== aging (L3 slow layer, hourly day profile) ==')
    a = AG
    print('  k_c=%.4e d^-1/2' % a['kc'])
    # calibration protocol C0: 25 C cell, 1 FCE/day (100->0 %), CC 1C to 80 %, 0.3C to 100 %, hold
    one = lambda q: 1.0
    h = age('full', 25.0, 1000, one, kH=1.0, kL=1.0, cal_mode=True)
    qc = h[-1][1]
    nH2 = h[-1][2] ** 2
    nL2 = h[-1][3] ** 2
    need = 0.20 - qc
    kH = (1 - a['wL']) * need / math.sqrt(nH2)
    kL = a['wL'] * need / math.sqrt(nL2)
    print('  C0: q_cal(1000 d)=%.4f  mean g=%.3f ; k_H=%.4e  k_L=%.4e per sqrt(FCE)' % (
        qc, qc / (a['kc'] * math.sqrt(1000)), kH, kL))
    chk = age('full', 25.0, 1000, one, kH=kH, kL=kL, cal_mode=True)[-1]
    print('  check SOH(1000 FCE)=%.4f' % (1 - chk[1] - chk[2] - chk[3]))
    Eday = 0.440 * 24 / 18.19  # typical-day energy fraction of nominal (S13 typical Pbar = 0.440 W)
    dfn = lambda q: min(Eday / (1 - q), 0.95)
    print('  g(0.1,0.5,0.7,0.8,0.9,1.0)=', [round(gz(z), 3) for z in (0.1, 0.5, 0.7, 0.8, 0.9, 1.0)])
    rows = (('A1 full', 'full', 25, {}), ('A2 limit80', 'limit80', 25, {}), ('A3 obc', 'obc', 25, {}),
            ('A4 hot 35C', 'full', 35, {}), ('A5 car 2h@45C', 'full', 25, dict(heat_h=2.0)),
            ('A6a fast 1.5C', 'full', 25, dict(c_cc=1.5)), ('A6b charge at 0C amb, 1C', 'full', 25, dict(Tch_amb=0.0)),
            ('A6c charge at 0C amb, 0.5C (JEITA)', 'full', 25, dict(Tch_amb=0.0, c_cc=0.5, c_cv=0.3)))
    for lab, strat, Ta, extra in rows:
        hh = age(strat, Ta, 3 * 365, dfn, kH=kH, kL=kL, **extra)
        sv = [1 - sum(hh[y * 365 - 1][1:]) for y in (1, 2, 3)]
        yrs = next((d for d, *qq in hh if 1 - sum(qq) <= 0.8), None)
        print('  %-34s SOH 1/2/3 yr = %.3f / %.3f / %.3f  (cal/H/L at 2 yr: %.3f %.3f %.3f) years to 0.8: %s' % (
            lab, *sv, *hh[2 * 365 - 1][1:], None if yrs is None else round(yrs / 365, 2)))
    dfn1 = lambda q: 0.95
    hh = age('full', 25, 3 * 365, dfn1, kH=kH, kL=kL)
    yrs = next((d for d, *qq in hh if 1 - sum(qq) <= 0.8), None)
    print('  A7 1 FCE/day (0.95 DoD)  SOH 1/2/3 yr = %.3f / %.3f / %.3f ; years to 0.8: %.2f' % (
        *tuple(1 - sum(hh[y * 365 - 1][1:]) for y in (1, 2, 3)), yrs / 365))


# --------------------------------------------------------------------------------------------- smooth lever metrics
UI = ('RD', 'SO', 'SO-L', 'NV')


def lever_maps():
    """Each lever maps (primitive name, u) -> u'.  Applied to S12 and to the CTMC state powers."""
    def scr(fn):
        return lambda n, u: fn(u) if u.get('scr') else u
    return {
        'L1 skip gaming (GM -> standby)': None,
        'L2 dark theme (APL x0.25 on UI states)': lambda n, u: dict(u, aplmul=0.25) if n in UI else u,
        'L3 brightness -0.1 (screen-on states)': scr(lambda u: dict(u, b=max(u['b'] - 0.1, 0.0))),
        'L4 refresh capped at 60 Hz': scr(lambda u: dict(u, fr=60)),
        'L5 video call on Wi-Fi instead of LTE': lambda n, u: dict(PR['VC-W']) if n == 'VC' else u,
        'L6 voice call instead of video call': lambda n, u: dict(u, scr=0, awake=1, camf=0, x=0.05) if n in ('VC', 'VC-W') else u,
        'L7 halve background wake-ups': lambda n, u: dict(u, lam=u['lam'] / 2) if n.startswith('SB') else u,
        'L8 standby on 5G NSA (no Wi-Fi)': lambda n, u: dict(PR['SB-5'], lam=u['lam']) if n.startswith('SB') else u,
        'L9 weak signal -110 dBm all day': lambda n, u: dict(u, rsrp=-110.0),
        'L10 brightness +0.2 (outdoor)': scr(lambda u: dict(u, b=min(u['b'] + 0.2, 1.0))),
        'O1 tail coalescing x5 (OS)': lambda n, u: dict(u, lam=u['lam'] / 5) if n.startswith('SB') else u,
        'O2 game frame cap: x_c 0.8 -> 0.64 (OS)': lambda n, u: dict(u, x=0.64) if n.startswith('GM') else u,
    }


def sp_u(u, p=P0, Tamb_C=25.0):
    Pns, dem, _, _ = power_parts(u, p)
    f = lambda Ts: Ts - (Tamb_C + K0) - p['kenv'] * p['Rskin'] * (Pns + throttle(dem, Ts, p))
    Ts = brentq(f, Tamb_C + K0 - 1, Tamb_C + K0 + 200)
    return Pns + throttle(dem, Ts, p)


def recs_smooth():
    print('== levers: smooth metrics (S12 mean-power TTE_inf = E_u/Pbar ; CTMC typical Dynkin mean) ==')
    Eu = Eu_Wh(P0)
    sched = [(60 * m, n) for m, n in [(30, 'RD'), (45, 'NV')] + [(26, 'SB-W'), (4, 'SO')] * 7 + [(30, 'VW'), (30, 'SO')] +
             [(26, 'SB-W'), (4, 'SO')] * 8 + [(45, 'NV'), (30, 'VC'), (60, 'GM'), (120, 'VW'), (60, 'SO'), (60, 'RD'), (480, 'SB-W')]]

    def s12_pbar(fn):
        E = 0.0
        for d, n in sched:
            u = PR['SB-W'] if (fn is None and n == 'GM') else (PR[n] if fn is None or fn == 'base' else fn(n, PR[n]))
            E += d * sp_u(u)
        return E / 86400.0

    def ctmc_mu(fn, prof='typical'):
        L = gen_product(prof)
        Pv = []
        for a in ACTS:
            for M in (MAPW, MAPL):
                n = M[a]
                if fn is None:
                    u = PR['SB-W'] if a == 'game' else PR[n]
                else:
                    u = PR[n] if fn == 'base' else fn(n, PR[n])
                Pv.append(sp_u(u))
        return moments(L, np.array(Pv), Eu)['mu_dynkin']
    P0b, mu0 = s12_pbar('base'), ctmc_mu('base')
    print('  base: S12 Pbar=%.4f W (E_u/Pbar=%.2f h)  CTMC typical Dynkin=%.2f h' % (P0b, Eu / P0b, mu0))
    for k, fn in lever_maps().items():
        Pb, mu = s12_pbar(fn), ctmc_mu(fn)
        print('  %-42s S12 %+6.1f%%   CTMC %+6.1f%%' % (k, 100 * (P0b / Pb - 1), 100 * (mu / mu0 - 1)))
    # O2 performance cost
    print('  O2 perf: f ratio 0.8 (x_c 0.64 vs 0.80); steady GM power %.3f -> %.3f W' % (
        sp_u(PR['GM']), sp_u(dict(PR['GM'], x=0.64))))
    # Wi-Fi share in CTMC 50 -> 80 %
    for mw, mc in ((2.5, 2.5), (4.0, 1.0)):
        L = gen_product('typical', mean_w=mw, mean_c=mc)
        Pv = state_powers(P0)
        print('  Wi-Fi share %.0f %%: CTMC Dynkin %.2f h' % (100 * mw / (mw + mc), moments(L, Pv, Eu)['mu_dynkin']))

if __name__ == '__main__':
    secs = sys.argv[1:] or ['base', 'scen', 'bonn', 'bench', 'cold', 'aged', 'ctmc', 'mc', 'recs', 'recs2', 'aging', 'sens']
    main(secs)
