"""Scratch prototype of the physics-first electro-thermal ECM (docs/design/proposal_physics_first.md).

Reference implementation used to compute the [PROTO] numbers in the proposal.  It is NOT the
production model (that lives in src/phonebatt).  State x = (z, v1, T); constant-power load built
from the component power model; terminal events: z<=0, V<=V_cut, V_pk<=V_sys.
Run:  python3 -I docs/design/scripts/pf_proto_model.py
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RG = 8.314
TREF = 298.15

P = dict(
    Qn=4.70,            # Ah
    b0=0.60117, b1=36.426, a0=3.6193, a1=0.44311, a2=-0.48605, a3=0.63127,
    dV_hv=0.242, p_hv=3.86, ocv_off=0.0,
    R0ref=0.060, A_R=0.80, z0_R=0.139, Ea0=18e3,
    R1ref=0.010, tau1ref=30.0, Ea1=30e3, Eatau=21e3,
    Vcut=3.40, Vsys=3.30,
    Cth=150.0, Rth=3.0, Rskin=5.5, Tthr=316.15, wthr=1.0, PSoCfloor=0.5,
    entropic=1.0,
    # device
    Psusp=0.020, Pawake=0.10,
    Pd0=0.15, Pcref=0.25, APLref=0.70, gam_b=2.0, kRR=0.0028,
    PSoCmax=5.0, v0frac=0.446,
    # radios: idle, tail, Ttail, Epro, e_b (W/Mbps), Prt (real-time)
    wifi=dict(idle=0.004, tail=0.12, Ttail=0.24, Epro=0.0098, eb=0.015, rt=0.30),
    lte=dict(idle=0.012, tail=0.15, Ttail=10.0, Epro=0.31, eb=0.030, rt=1.20),
    nr=dict(idle=0.020, tail=0.25, Ttail=10.3, Epro=1.0, eb=0.030, rt=1.40),
    kul=5.0, gsig=0.25, rsrp_ref=-85.0, msig_max=3.5,
    Pgnss=0.25, Pcam=1.4, Pcamf=0.30, Paudio=0.05, Pbt=0.04,
    qloss=0.0, rgrowth=0.0,
)

DUDT_Z = np.array([0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00])
DUDT_V = np.array([-0.4485, -0.706, -0.534, -0.479, -0.372, -0.307, -0.248, -0.248, -0.205, -0.153, -0.097, -0.046]) * 1e-3


def ocv(z, T, p=P):
    z = min(max(z, 0.0), 1.0)
    u = (-p['b0'] * math.exp(-p['b1'] * z) + p['a0'] + p['a1'] * z + p['a2'] * z**2 + p['a3'] * z**3
         + p['dV_hv'] * z ** p['p_hv'] + p['ocv_off'])
    return u + (T - TREF) * float(np.interp(z, DUDT_Z, DUDT_V))


def arr(Ea, T):
    return math.exp(Ea / RG * (1.0 / T - 1.0 / TREF))


def R0f(z, T, p=P):
    z = max(z, 0.0)
    Rinf = p['R0ref'] / (1 + p['A_R'] * math.exp(-0.5 / p['z0_R']))
    return Rinf * (1 + p['A_R'] * math.exp(-z / p['z0_R'])) * arr(p['Ea0'], T) * (1 + p['rgrowth'])


def R1f(T, p=P):
    return p['R1ref'] * arr(p['Ea1'], T) * (1 + p['rgrowth'])


def tau1f(T, p=P):
    return p['tau1ref'] * arr(p['Eatau'], T)


def z_empty(p=P):
    return brentq(lambda z: ocv(z, TREF, p) - p['Vcut'], 1e-6, 0.2)


def msig(rsrp, p=P):
    return min(p['msig_max'], 1 + p['gsig'] * max(0.0, (p['rsrp_ref'] - rsrp) / 10.0))


def radio_traffic(r, lam, ddl, dul, kul, m):
    tailfrac = 1 - math.exp(-lam * r['Ttail'])
    return m * ((r['tail'] - r['idle']) * tailfrac + lam * math.exp(-lam * r['Ttail']) * r['Epro']
                + r['eb'] * (ddl + kul * dul))


def power(u, Ts, p=P):
    """Return (P_mean, P_burst, breakdown dict) at the battery terminals [W]."""
    br = {}
    br['susp'] = p['Psusp']
    br['awake'] = p['Pawake'] * u.get('awake', u.get('scr', 0))
    scr = u.get('scr', 0)
    b = u.get('b', 0.4)
    apl = u.get('apl', 0.7)
    br['disp'] = scr * (p['Pd0'] + p['Pcref'] * (b / 0.5) ** p['gam_b'] * apl / p['APLref']
                        + p['kRR'] * max(u.get('frr', 60) - 60, 0))
    x = u.get('x', 0.0)
    dem = p['PSoCmax'] * x * (p['v0frac'] + (1 - p['v0frac']) * x) ** 2
    if dem > p['PSoCfloor']:
        sig = 1.0 / (1.0 + math.exp(-(p['Tthr'] - Ts) / p['wthr']))
        br['soc'] = p['PSoCfloor'] + (dem - p['PSoCfloor']) * sig
    else:
        br['soc'] = dem
    # radios
    net = 0.0
    if u.get('wifi', 1):
        net += p['wifi']['idle']
    cell = u.get('cell', 'lte')
    m = msig(u.get('rsrp', -85.0), p)
    if cell:
        net += p[cell]['idle']
    via = u.get('via', 'wifi')
    cls = u.get('cls', 'burst')
    r = p['wifi'] if via == 'wifi' else p[cell] if via == 'cell' else None
    mm = 1.0 if via == 'wifi' else m
    if r is not None:
        if cls == 'rt':
            net += mm * r['rt']
        else:
            net += radio_traffic(r, u.get('lam', 0.0), u.get('ddl', 0.0), u.get('dul', 0.0), p['kul'], mm)
    br['net'] = net
    br['gnss'] = p['Pgnss'] * u.get('gnss', 0.0)
    br['cam'] = p['Pcam'] * u.get('cam', 0.0) + p['Pcamf'] * u.get('camf', 0.0)
    br['aud'] = p['Paudio'] * u.get('audio', 0.0) + p['Pbt'] * u.get('bt', 0.0)
    br['sens'] = u.get('sens', 0.0)
    Pm = sum(br.values())
    return Pm, Pm + u.get('burst', 0.0), br


def make_rhs(u, Tamb, p=P):
    Q = p['Qn'] * (1 - p['qloss'])

    def state(x):
        z, v1, T = x
        U = ocv(z, T, p)
        R0 = R0f(z, T, p)
        Ts = Tamb + p['Rskin'] / p['Rth'] * (T - Tamb)
        Pm, Pk, _ = power(u, Ts, p)
        Ue = U - v1
        D = Ue * Ue - 4 * R0 * Pm
        I = 2 * Pm / (Ue + math.sqrt(D)) if D > 0 else Ue / (2 * R0)
        V = Ue - I * R0
        Dk = Ue * Ue - 4 * R0 * Pk
        Vk = 0.5 * (Ue + math.sqrt(Dk)) if Dk > 0 else 0.0
        return U, R0, I, V, Vk, Pm, D

    def rhs(t, x):
        z, v1, T = x
        U, R0, I, V, Vk, Pm, D = state(x)
        R1 = R1f(T, p)
        tau = tau1f(T, p)
        dz = -I / (3600.0 * Q)
        dv = (R1 * I - v1) / tau
        dudt = float(np.interp(min(max(z, 0), 1), DUDT_Z, DUDT_V))
        q = I * V + I * I * R0 + v1 * v1 / R1 - p['entropic'] * I * T * dudt
        dT = (q - (T - Tamb) / p['Rth']) / p['Cth']
        return [dz, dv, dT]

    def ev_z(t, x):
        return x[0]
    ev_z.terminal = True; ev_z.direction = -1

    def ev_v(t, x):
        return state(x)[3] - p['Vcut']
    ev_v.terminal = True; ev_v.direction = -1

    def ev_b(t, x):
        return state(x)[4] - p['Vsys']
    ev_b.terminal = True; ev_b.direction = -1

    return rhs, [ev_z, ev_v, ev_b], state


def run(schedule, Tamb, s0=1.0, T0=None, p=P, tmax_h=2000.0, cyc=True):
    """schedule: list of (duration_s, u). Repeats cyclically until an event."""
    ze = z_empty(p)
    z0 = ze + s0 * (1 - ze)
    x = [z0, 0.0, Tamb if T0 is None else T0]
    t = 0.0
    names = ['z<=0', 'V<=Vcut', 'Vpk<=Vsys']
    k = 0
    traj = []
    while t < tmax_h * 3600:
        dur, u = schedule[k % len(schedule)]
        Ta = u.get('Tamb', Tamb)
        rhs, evs, state = make_rhs(u, Ta, p)
        sol = solve_ivp(rhs, (t, t + dur), x, method='LSODA', events=evs, rtol=1e-7,
                        atol=[1e-9, 1e-7, 1e-5], max_step=60.0)
        traj.append(sol)
        for i, te in enumerate(sol.t_events):
            if len(te):
                xe = sol.y_events[i][0]
                return te[0] / 3600, names[i], xe, (xe[0] - ze) / (1 - ze), traj
        x = sol.y[:, -1]
        t = sol.t[-1]
        k += 1
        if not cyc and k >= len(schedule):
            break
    return np.inf, 'none', x, None, traj


# ------------------------------------------------------------------ primitives
U = {}
U['standby'] = dict(scr=0, awake=0, wifi=1, cell='lte', via='cell', lam=12 / 3600, ddl=0.0, burst=1.0)
U['standby5g'] = dict(scr=0, awake=0, wifi=0, cell='nr', via='cell', lam=60 / 3600, burst=1.0)
U['read'] = dict(scr=1, b=0.4, apl=0.7, frr=60, x=0.05, wifi=1, cell='lte', via='wifi', lam=1 / 30, ddl=0.05, burst=3.0)
U['social'] = dict(scr=1, b=0.4, apl=0.7, frr=120, x=0.25, wifi=1, cell='lte', via='wifi', lam=0.5, ddl=2.0, burst=3.0)
U['social_lte'] = dict(U['social'], via='cell', wifi=0)
U['video_wifi'] = dict(scr=1, b=0.4, apl=0.35, frr=60, x=0.10, wifi=1, cell='lte', via='wifi', lam=0.25, ddl=3.0, audio=1, burst=2.0)
U['video_5g'] = dict(U['video_wifi'], wifi=0, cell='nr', via='cell')
U['video_lte'] = dict(U['video_wifi'], wifi=0, cell='lte', via='cell')
U['nav'] = dict(scr=1, b=0.8, apl=0.6, frr=60, x=0.20, wifi=0, cell='lte', via='cell', lam=0.2, ddl=0.1, gnss=1.0, sens=0.047, burst=3.0)
U['game'] = dict(scr=1, b=0.6, apl=0.5, frr=120, x=0.80, wifi=1, cell='lte', via='wifi', lam=20.0, ddl=0.2, dul=0.05, audio=1, burst=4.0)
U['vcall'] = dict(scr=1, b=0.5, apl=0.5, frr=60, x=0.30, wifi=0, cell='lte', via='cell', cls='rt', camf=1, audio=1, burst=3.0)
U['camera'] = dict(scr=1, b=0.5, apl=0.5, frr=60, x=0.05, wifi=1, cell='lte', via='wifi', lam=0.0, cam=1, burst=4.0)
U['pocket'] = dict(U['standby'])


def report(name, sched, Tamb=298.15, s0=1.0, T0=None, p=P):
    tte, ev, xe, sdisp, traj = run(sched, Tamb, s0, T0, p)
    Pm = power(sched[0][1], Tamb)[0]
    print(f"{name:34s} P0={Pm:6.3f} W  TTE={tte:8.2f} h  event={ev:10s} z_end={xe[0]:.3f} s_disp_end={sdisp*100 if sdisp is not None else float('nan'):5.1f}%  T_end={xe[2]-273.15:5.1f} C")
    return tte, ev, xe, sdisp


if __name__ == '__main__':
    print('z_empty', z_empty())
    print('U(1)', ocv(1, TREF), 'mean OCV', np.mean([ocv(z, TREF) for z in np.linspace(0, 1, 2001)]))
    for k, u in U.items():
        Pm, Pk, br = power(u, 30.0)
        print(f"{k:12s} P={Pm:6.3f}  " + ' '.join(f"{a}={b:.3f}" for a, b in br.items()))
    print()
    K = 273.15
    report('S1 standby LTE+WiFi', [(3600, U['standby'])])
    report('S2 standby 5G heavy push', [(3600, U['standby5g'])])
    report('S3 reading WiFi', [(3600, U['read'])])
    report('S4 social WiFi 120Hz', [(3600, U['social'])])
    report('S5a video WiFi', [(3600, U['video_wifi'])])
    report('S5b video 5G', [(3600, U['video_5g'])])
    report('S5c video LTE', [(3600, U['video_lte'])])
    report('S6 navigation', [(3600, U['nav'])])
    report('S7 gaming 25C', [(3600, U['game'])])
    report('S7 gaming 35C', [(3600, U['game'])], Tamb=35 + K)
    report('S8 video call LTE', [(3600, U['vcall'])])
    report('camera', [(3600, U['camera'])])
    # S9 cold outdoor: nav 10 min + photo 1 min + pocket? keep in hand at -10C
    cold = [(600, U['nav']), (60, dict(U['camera'], b=0.9)), (1140, U['nav'])]
    report('S9 cold -10C nav+photo T0=20C', cold, Tamb=-10 + K, T0=20 + K)
    report('S9 same at 25C', cold, Tamb=25 + K)
    for TT in (0, -10, -20):
        report(f'social at {TT}C', [(3600, U['social'])], Tamb=TT + K)
    # aged
    pa = dict(P, qloss=0.2, rgrowth=0.5)
    report('S10 aged social 25C', [(3600, U['social'])], p=pa)
    report('S10 aged social 0C', [(3600, U['social'])], Tamb=0 + K, p=pa)
    report('S10 aged social -10C', [(3600, U['social'])], Tamb=-10 + K, p=pa)
    report('S10 aged nav -10C', cold, Tamb=-10 + K, T0=20 + K, p=pa)
    # weak signal
    report('S11 social LTE -85', [(3600, U['social_lte'])])
    report('S11 social LTE -110', [(3600, dict(U['social_lte'], rsrp=-110))])
    report('S11 standby LTE -110', [(3600, dict(U['standby'], rsrp=-110))])
