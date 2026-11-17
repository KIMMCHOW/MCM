"""Checks for the [PROTO] numbers in docs/design/proposal_physics_first.md.

Usage:  python3 -I docs/design/scripts/pf_proto_checks.py {day|ctmc|quad|bonn|lee|elast|cold|rec|all}
  day   S12 typical-day TTE (from 100/50/20 %, aged, 0 C)
  ctmc  S12-R: renewal-CLT mean/sd of TTE vs energy-bucket Monte Carlo
  quad  quasi-steady quadrature formula F2 vs full ODE
  bonn  blind predictions for the Bonn iPhone 14 Pro / 13 Pro runs
  lee   qualitative 2015-class high-resistance pack in the cold (V7)
  elast local elasticities (S4, S9-aged, S7)
  cold  aged-cold elasticities, peak-power governor, stranded-charge map
  rec   behaviour ranking on the S12 day (mean-demand based)
"""
import math, sys, time
import numpy as np
from scipy.optimize import brentq
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pf_proto_model as M
from pf_proto_model import P, U, run, power, ocv, R0f, R1f, z_empty, TREF

K = 273.15
what = sys.argv[1] if len(sys.argv) > 1 else 'all'


def tte(sched, Tamb=TREF, s0=1.0, T0=None, p=P):
    return run(sched, Tamb, s0, T0, p)[:4]


# (a) day in the life -------------------------------------------------------
def day_schedule():
    sb = U['standby']; so = U['social']; rd = U['read']; vw = U['video_wifi']
    nav = U['nav']; vc = dict(U['vcall']); gm = U['game']
    work = [(26 * 60, sb), (4 * 60, so)] * 7
    aft = [(26 * 60, sb), (4 * 60, so)] * 8
    s = ([(30 * 60, rd), (45 * 60, nav)] + work + [(30 * 60, vw), (30 * 60, so)] + aft
         + [(45 * 60, nav), (30 * 60, vc), (60 * 60, gm), (30 * 60, vw), (90 * 60, vw), (60 * 60, so), (60 * 60, rd)]
         + [(8 * 3600, sb)])
    return s


if what in ('all', 'day'):
    s = day_schedule()
    tot = sum(d for d, _ in s)
    E = sum(d * power(u, 300.0)[0] for d, u in s)
    print('day schedule length h', tot / 3600, 'mean P W', E / tot)
    for s0 in (1.0, 0.5, 0.2):
        r = tte(s, s0=s0)
        print(f'day-in-life s0={s0}: TTE={r[0]:.2f} h event={r[1]} sdisp_end={r[3]}')
    pa = dict(P, qloss=0.2, rgrowth=0.5)
    r = tte(s, p=pa); print('day aged 25C', r[0], r[1])
    r = tte(s, Tamb=0 + K, p=pa); print('day aged 0C', r[0], r[1])
    r = tte(s, Tamb=0 + K); print('day fresh 0C', r[0], r[1])

# (b) CTMC analytic vs MC (energy bucket) -----------------------------------
if what in ('all', 'ctmc'):
    names = ['standby', 'read', 'social', 'video_wifi', 'nav', 'game', 'vcall']
    hold = np.array([25, 6, 8, 20, 30, 25, 10]) * 60.0  # s
    Pst = np.array([power(U[n], 300.0)[0] for n in names])
    # jump matrix
    J = np.zeros((7, 7))
    J[0, 1:] = [0.35, 0.35, 0.12, 0.08, 0.05, 0.05]
    for i in range(1, 7):
        J[i, 0] = 0.7
        others = [j for j in range(1, 7) if j != i]
        w = J[0, others] / J[0, others].sum()
        J[i, others] = 0.3 * w
    Qg = (J - np.eye(7)) / hold[:, None]
    # stationary
    A = np.vstack([Qg.T, np.ones(7)])
    pi = np.linalg.lstsq(A, np.r_[np.zeros(7), 1.0], rcond=None)[0]
    Pbar = pi @ Pst
    # Poisson eq Q g = -(P - Pbar), pi g = 0
    Ab = np.vstack([Qg, pi])
    g = np.linalg.lstsq(Ab, np.r_[-(Pst - Pbar), 0.0], rcond=None)[0]
    sig2 = 2 * np.sum(pi * (Pst - Pbar) * g)
    ze = z_empty()
    zz = np.linspace(ze, 1, 4001)
    Eu = 3600 * P['Qn'] * np.trapezoid([ocv(z, TREF) for z in zz], zz)  # J (OCV energy, upper bound)
    Eu_eff = Eu * 0.985
    mu = Eu_eff / Pbar; sd = math.sqrt(sig2 * Eu_eff / Pbar**3)
    print('pi', np.round(pi, 3), 'Pbar', Pbar, 'sigma2_inf W^2 s', sig2, 'E_ocv Wh', Eu / 3600)
    print(f'analytic TTE mean {mu/3600:.2f} h sd {sd/3600:.2f} h CV {sd/mu:.3f}')
    rng = np.random.default_rng(1)
    out = []
    for k in range(4000):
        st = 0; e = 0.0; t = 0.0
        st = rng.choice(7, p=pi)
        while True:
            dt = rng.exponential(hold[st])
            if e + Pst[st] * dt >= Eu_eff:
                t += (Eu_eff - e) / Pst[st]; break
            e += Pst[st] * dt; t += dt
            st = rng.choice(7, p=J[st])
        out.append(t)
    out = np.array(out)
    print(f'MC energy-bucket TTE mean {out.mean()/3600:.2f} h sd {out.std()/3600:.2f} h; 5-95% {np.percentile(out,5)/3600:.2f}-{np.percentile(out,95)/3600:.2f}')
    print('state powers', dict(zip(names, np.round(Pst, 3))))

# (c) quasi-steady quadrature vs ODE ----------------------------------------
if what in ('all', 'quad'):
    def tte_quad(Pw, Tamb=TREF, p=P, Rextra=0.0):
        Q = p['Qn'] * (1 - p['qloss'])
        def Vqs(z):
            Rs = R0f(z, Tamb, p) + R1f(Tamb, p)
            Uq = ocv(z, Tamb, p)
            D = Uq * Uq - 4 * Rs * Pw
            return 0.5 * (Uq + math.sqrt(D)) if D > 0 else 0.0
        ze = brentq(lambda z: Vqs(z) - p['Vcut'], 1e-4, 0.9)
        zz = np.linspace(ze, 1.0, 4001)
        z0 = z_empty(p)
        return 3600 * Q / Pw * np.trapezoid([Vqs(z) for z in zz], zz) / 3600, ze
    for name in ('read', 'social', 'video_wifi', 'nav', 'vcall'):
        Pw = power(U[name], 300.0)[0]
        uu = dict(U[name], burst=0.0)
        r = tte([(3600, uu)])
        q, ze = tte_quad(Pw)
        # isothermal ODE comparison: large Cth? just report
        print(f'{name:10s} P={Pw:.3f} ODE(no burst)={r[0]:.3f} h  quad(isothermal 25C)={q:.3f} h  rel={(q-r[0])/r[0]*100:+.2f}%  energy-bucket={18.19/Pw:.2f} h')

# (d) Bonn predictions ------------------------------------------------------
if what in ('all', 'bonn'):
    for lab, Q, obs in (('14 Pro', 3.32, dict(idle=55.0, cam=8.43, gps=48.5, cxl=8.38)), ('13 Pro', 3.16, dict(idle=50.8, cam=9.53, gps=29.8, cxl=8.33))):
        pb = dict(P, Qn=Q, R0ref=0.060 * 4.7 / Q, R1ref=0.010 * 4.7 / Q)
        base = dict(scr=1, b=0.0, apl=0.0, frr=60, x=0.0, wifi=0, cell=None, via=None, burst=0.0)
        sc = dict(idle=base, cam=dict(base, cam=1), gps=dict(base, gnss=1), cxl=dict(base, cam=1, gnss=1))
        res = {}
        for k, u in sc.items():
            r = tte([(3600, u)], p=pb)
            res[k] = r[0]
            print(f'Bonn {lab} {k:5s} P={power(u,300.0,pb)[0]:.3f} W pred TTE={r[0]:.2f} h obs={obs[k]} h err={(r[0]-obs[k])/obs[k]*100:+.1f}%')
        # display-hold mapping not needed for TTE

# (e) LEE20-style older phone cold check -------------------------------------
if what in ('all', 'lee'):
    pl = dict(P, Qn=2.7, R0ref=0.12, R1ref=0.03)
    u = dict(scr=1, b=0.5, apl=0.7, frr=60, x=0.75, wifi=1, cell='lte', via='wifi', lam=0.5, ddl=1.0, burst=3.0)
    for T in (25, 0):
        r = tte([(3600, u)], Tamb=T + K, p=pl)
        print('LEE-like', T, 'C', 'P', power(u, 300.0, pl)[0], 'TTE min', r[0] * 60, r[1], 'sdisp_end', r[3])

# (f) elasticities ----------------------------------------------------------
if what in ('all', 'elast'):
    cases = {
        'social25': ([(3600, U['social'])], TREF, None, {}),
        'nav-10aged': ([(600, U['nav']), (60, dict(U['camera'], b=0.9)), (1140, U['nav'])], -10 + K, 20 + K, dict(qloss=0.2, rgrowth=0.5)),
        'game25': ([(3600, U['game'])], TREF, None, {}),
    }
    keys = ['Qn', 'R0ref', 'R1ref', 'Ea0', 'Vcut', 'Vsys', 'Pd0', 'Pcref', 'gam_b', 'PSoCmax', 'Pawake', 'Rth', 'Rskin', 'Cth', 'ocv_off']
    for cname, (sched, Ta, T0, mod) in cases.items():
        base_p = dict(P, **mod)
        t0 = tte(sched, Ta, 1.0, T0, base_p)[0]
        line = []
        for k in keys:
            h = 0.05
            pp = dict(base_p)
            if k == 'ocv_off':
                pp[k] = 0.010; tt = tte(sched, Ta, 1.0, T0, pp)[0]; line.append(f'{k}:{(tt-t0)/t0*100:+.2f}%/10mV'); continue
            pp[k] = base_p[k] * (1 + h)
            tt = tte(sched, Ta, 1.0, T0, pp)[0]
            line.append(f'{k}:{(tt - t0) / t0 / h:+.3f}')
        # burst
        sched2 = [(d, dict(u, burst=u.get('burst', 0) * 1.05)) for d, u in sched]
        tt = tte(sched2, Ta, 1.0, T0, base_p)[0]; line.append(f'burst:{(tt-t0)/t0/0.05:+.3f}')
        print(cname, f'TTE0={t0:.2f}', ' '.join(line))

# (g) cold: aged-cold elasticities, peak governor, stranded-charge map -------
if what in ('all', 'cold'):
    sch = None
    pa=dict(P,qloss=0.2,rgrowth=0.5)
    sch=[(3600,U['social'])]
    t0=tte(sch,-10+K,p=pa)[0]
    out=[]
    for k in ['Qn','R0ref','R1ref','Ea0','Ea1','Pd0','PSoCmax','Rth','Cth']:
        pp=dict(pa); pp[k]=pa[k]*1.05; out.append(f"{k}:{(tte(sch,-10+K,p=pp)[0]-t0)/t0/0.05:+.3f}")
    for k,dv in [('Vcut',0.01),('Vsys',0.01),('ocv_off',0.01)]:
        pp=dict(pa); pp[k]=pa[k]+dv; out.append(f"{k}:{(tte(sch,-10+K,p=pp)[0]-t0)/t0*100:+.2f}%/10mV")
    s2=[(3600,dict(U['social'],burst=3.15))]; out.append(f"burst:{(tte(s2,-10+K,p=pa)[0]-t0)/t0/0.05:+.3f}")
    print('aged social -10C TTE0',t0,' '.join(out))
    # peak governor: cap burst so that P+dPb <= 0.9 Pcut,sys -> emulate by burst=1.5 W
    for b in (3.0,2.0,1.5,1.0,0.0):
        r=tte([(3600,dict(U['social'],burst=b))],-10+K,p=pa)
        print(f'aged social -10C burst={b}: TTE={r[0]:.2f} h event={r[1]} s_strand={r[3]*100:.1f}%')
    for T in (25,10,0,-10,-20):
        for soh,rg in ((0.0,0.0),(0.1,0.25),(0.2,0.5)):
            r=tte([(3600,U['social'])],T+K,p=dict(P,qloss=soh,rgrowth=rg))
            print(f'map T={T:4d} qloss={soh:.1f} TTE={r[0]:6.2f} ev={r[1]:9s} s_strand={r[3]*100:5.1f}%')

# (h) behaviour ranking on the S12 day (mean demand) ---------------------------
if what in ('all', 'rec'):
    S = day_schedule()

    def pbar(sched, p=P, Ts=300.0):
        tot = sum(d for d, _ in sched)
        return sum(d * power(u, Ts, p)[0] for d, u in sched) / tot

    P0 = pbar(S)

    def on_screen(fn):
        return [(d, fn(dict(u)) if u.get('scr', 0) else u) for d, u in S]

    tests = {
        'brightness -0.1': on_screen(lambda u: dict(u, b=max(u['b'] - 0.1, 0.05))),
        'brightness cap 0.5': on_screen(lambda u: dict(u, b=min(u['b'], 0.5))),
        'dark mode (APL x0.25 on UI)': on_screen(lambda u: dict(u, apl=u['apl'] * 0.25 if u['apl'] >= 0.6 else u['apl'])),
        '60 Hz cap': on_screen(lambda u: dict(u, frr=60)),
        'Wi-Fi instead of LTE for nav & call': on_screen(lambda u: dict(u, via='wifi', wifi=1) if u.get('via') == 'cell' else u),
        'background sync halved': [(d, dict(u, lam=u['lam'] / 2) if not u.get('scr', 0) else u) for d, u in S],
        'standby on 5G NSA': [(d, U['standby5g'] if not u.get('scr', 0) else u) for d, u in S],
        'GNSS off': [(d, dict(u, gnss=0) if u.get('gnss') else u) for d, u in S],
        'skip 1 h gaming': [(d, U['standby'] if u is U['game'] else u) for d, u in S],
        'weak signal -110 dBm': [(d, dict(u, rsrp=-110)) for d, u in S],
        'voice instead of video call': [(d, dict(u, camf=0, scr=0, awake=1, x=0.05) if u.get('cls') == 'rt' else u) for d, u in S],
    }
    print('baseline Pbar', P0)
    for k, s2 in tests.items():
        p2 = pbar(s2)
        print(f'{k:40s} Pbar={p2:.4f} W  dTTE={(P0 / p2 - 1) * 100:+.1f}%')
    comp = {}
    for d, u in S:
        for k, v in power(u, 300.0)[2].items():
            comp[k] = comp.get(k, 0) + d * v
    tot = sum(comp.values())
    print('component energy shares %', {k: round(v / tot * 100, 1) for k, v in comp.items()})
